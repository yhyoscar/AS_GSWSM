import numpy as np
from netCDF4 import Dataset as netcdf
from os import system
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from namelist import *
from mod_commons import *

nstepfoutvar = foutvarfreq
if foutvarfreq<0:   nstepfoutvar = int(-foutvarfreq*tunit/timestep)

nstepfoutfile = foutfilefreq
if foutfilefreq<0:   nstepfoutfile = int(-foutfilefreq*tunit/timestep)

nstepfig  = figfreq
if figfreq<0:    nstepfig  = int(-figfreq*tunit/timestep)

class NCfile:
    def __init__(self):
        self.time = np.array([])
        self.hs  = np.zeros([1,nlat,nlon])   # topography (m)
        self.hm  = np.zeros([1])             # mean thickness (m)
        self.h   = np.zeros([1,nlat,nlon])   # free surface height (m)
        self.u   = np.zeros([1,nlat,nlon])   # zonal wind (m/s)
        self.v   = np.zeros([1,nlat,nlon])   # meridional wind (m/s)
        self.vor = np.zeros([1,nlat,nlon])
        self.div = np.zeros([1,nlat,nlon])
        self.str = np.zeros([1,nlat,nlon])
        self.pot = np.zeros([1,nlat,nlon])
        self.sp_vor = np.zeros([1,2,M+1,N+1])  # vorticity in spetral space (real, image)
        self.sp_div = np.zeros([1,2,M+1,N+1])  # vorticity in spetral space (real, image)
        
    def addatm(self,atm,coord,sec):
        atm.gps.stog(coord)
        atm.gp.stog(coord)
        atm.svor_to_sstr()  
        atm.sdiv_to_spot()
        atm.vor.stog(coord)
        atm.div.stog(coord)
        atm.str.stog(coord)
        atm.pot.stog(coord)
        atm.U.g, atm.V.g = fvordivstouvg(atm.vor.s, atm.div.s,coord)

        if len(self.time)==0:
            self.time = np.array([sec])
            self.hs[0,:,:] = atm.gps.g / grav
            self.hm[0]     = atm.gpm / grav
            self.h[0,:,:]  = (atm.gp.g + atm.gps.g + atm.gpm) / grav
            self.u[0,:,:]  = atm.U.g /coord.clat2d
            self.v[0,:,:]  = atm.V.g /coord.clat2d
            self.vor[0,:,:] = atm.vor.g + 0
            self.div[0,:,:] = atm.div.g + 0
            self.str[0,:,:] = atm.str.g + 0
            self.pot[0,:,:] = atm.pot.g + 0
            self.sp_vor[0,0,:,:] = atm.vor.s.real + 0
            self.sp_vor[0,1,:,:] = atm.vor.s.imag + 0
            self.sp_div[0,0,:,:] = atm.div.s.real + 0
            self.sp_div[0,1,:,:] = atm.div.s.imag + 0
        else:
            self.time = np.append(self.time, sec)   
            self.hs = np.append(self.hs, np.reshape(atm.gps.g/grav, (1,nlat,nlon)),axis=0)
            self.hm = np.append(self.hm, np.array([atm.gpm/grav]),axis=0)
            self.h  = np.append(self.h, np.reshape((atm.gp.g+atm.gps.g+atm.gpm)/grav, (1,nlat,nlon)),axis=0)
            self.u  = np.append(self.u, np.reshape(atm.U.g /coord.clat2d, (1,nlat,nlon)),axis=0)
            self.v  = np.append(self.v, np.reshape(atm.V.g /coord.clat2d, (1,nlat,nlon)),axis=0)
            self.vor = np.append(self.vor, np.reshape(atm.vor.g, (1,nlat,nlon)),axis=0)
            self.div = np.append(self.div, np.reshape(atm.div.g, (1,nlat,nlon)),axis=0)
            self.str = np.append(self.str, np.reshape(atm.str.g, (1,nlat,nlon)),axis=0)
            self.pot = np.append(self.pot, np.reshape(atm.pot.g, (1,nlat,nlon)),axis=0)
            self.sp_vor = np.append(self.sp_vor, np.reshape(np.array([atm.vor.s.real, atm.vor.s.imag]), (1,2,M+1,N+1)),axis=0)
            self.sp_div = np.append(self.sp_div, np.reshape(np.array([atm.div.s.real, atm.div.s.imag]), (1,2,M+1,N+1)),axis=0)

    def writefile(self,fnout,coord):
        print fnout
        system('rm -f '+fnout)
        fidout = netcdf(fnout, 'w')

        fidout.createDimension('time', len(self.time))
        fidout.createDimension('lat', nlat)
        fidout.createDimension('lon', nlon)
        fidout.createDimension('n', N+1)
        fidout.createDimension('m', M+1)
        fidout.createDimension('complex', 2)

        ftime = fidout.createVariable('time', 'f', ('time',))
        ftime[:] = self.time + 0
        ftime.units = 'second'
        ftime.long_name = 'time'
        
        flat = fidout.createVariable('lat', 'f', ('lat',))
        flat[:] = coord.lat + 0
        flat.units = 'degrees_north'
        flat.long_name = 'latitude'

        fwt = fidout.createVariable('wt', 'f', ('lat',))
        fwt[:] = coord.wt + 0
        fwt.units = 'degrees_north'
        fwt.long_name = 'Gaussian weights'

        flon = fidout.createVariable('lon', 'f', ('lon',))
        flon[:] = coord.lon + 0
        flon.units = 'degrees_east'
        flon.long_name = 'longitude'

        fm = fidout.createVariable('m', 'f', ('m',))
        fm[:] = np.linspace(0,M,M+1)
        fm.units = '#'
        fm.long_name = 'zonal wavenumber'

        fn = fidout.createVariable('n', 'f', ('n',))
        fn[:] = np.linspace(0,N,N+1)
        fn.units = '#'
        fn.long_name = 'total wavenumber'

        fn = fidout.createVariable('complex', 'f', ('complex',))
        fn[:] = np.array([0,1])
        fn.units = '#'
        fn.long_name = 'real and image parts of complex numbers'

        if 'hm' in foutvars:
            fdata = fidout.createVariable('hm', 'f', ('time',))
            fdata[:] = self.hm + 0
            fdata.units = 'm'
            fdata.long_name = 'mean of free surface height'
        if 'hs' in foutvars:
            fdata = fidout.createVariable('hs', 'f', ('time','lat','lon',))
            fdata[:] = self.hs + 0
            fdata.units = 'm'
            fdata.long_name = 'surface height'
        if 'h' in foutvars:
            fdata = fidout.createVariable('h', 'f', ('time','lat','lon',))
            fdata[:] = self.h + 0
            fdata.units = 'm'
            fdata.long_name = 'free surface height'
        if 'u' in foutvars:
            fdata = fidout.createVariable('u', 'f', ('time','lat','lon',))
            fdata[:] = self.u + 0
            fdata.units = 'm/s'
            fdata.long_name = 'zonal wind'
        if 'v' in foutvars:
            fdata = fidout.createVariable('v', 'f', ('time','lat','lon',))
            fdata[:] = self.v + 0
            fdata.units = 'm/s'
            fdata.long_name = 'meridional wind'
        if 'vor' in foutvars:
            fdata = fidout.createVariable('vor', 'f', ('time','lat','lon',))
            fdata[:] = self.vor + 0
            fdata.units = '1/s'
            fdata.long_name = 'relative vorticity'
        if 'div' in foutvars:
            fdata = fidout.createVariable('div', 'f', ('time','lat','lon',))
            fdata[:] = self.div + 0
            fdata.units = '1/s'
            fdata.long_name = 'divergence'
        if 'str' in foutvars:
            fdata = fidout.createVariable('str', 'f', ('time','lat','lon',))
            fdata[:] = self.str + 0
            fdata.units = 'm2/s'
            fdata.long_name = 'stream function'
        if 'pot' in foutvars:
            fdata = fidout.createVariable('pot', 'f', ('time','lat','lon',))
            fdata[:] = self.pot + 0
            fdata.units = 'm2/s'
            fdata.long_name = 'potential function'
        if 'sp_vor' in foutvars:
            fdata = fidout.createVariable('sp_vor', 'f', ('time','complex','m','n',))
            fdata[:] = self.sp_vor + 0
            fdata.units = '#'
            fdata.long_name = 'spherical harmonic coefficients: relative vorticity'
        if 'sp_div' in foutvars:
            fdata = fidout.createVariable('sp_div', 'f', ('time','complex','m','n',))
            fdata[:] = self.sp_div + 0
            fdata.units = '#'
            fdata.long_name = 'spherical harmonic coefficients: divergence'

        self.__init__()


def fplotatm(ffig,tstr,atm,coord,varstrs=figvars):
    print ffig
    nvar = len(varstrs);  ncol = 2;  nrow = np.ceil(1.0*nvar/ncol)
    plt.figure(1, figsize=[6*ncol,4*nrow])

    atm.gps.stog(coord)
    atm.gp.stog(coord)
    atm.svor_to_sstr()  
    atm.sdiv_to_spot()
    atm.vor.stog(coord)
    atm.div.stog(coord)
    atm.str.stog(coord)
    atm.pot.stog(coord)
    atm.U.g, atm.V.g = fvordivstouvg(atm.vor.s, atm.div.s,coord)

    for ivar in range(nvar):
        if varstrs[ivar] == 'hs':  data = atm.gps.g/grav;  x=coord.lon; y=coord.lat; clev=np.linspace(0,5000,41); 
        if varstrs[ivar] == 'h':   data = (atm.gp.g+atm.gps.g+atm.gpm)/grav;  x=coord.lon; y=coord.lat; clev=np.linspace(7700,8300,61); 
        if varstrs[ivar] == 'u':   data = atm.U.g /coord.clat2d; x=coord.lon; y=coord.lat; clev = np.linspace(-20,20,41);
        if varstrs[ivar] == 'v':   data = atm.V.g /coord.clat2d; x=coord.lon; y=coord.lat; clev = np.linspace(-20,20,41)
        if varstrs[ivar] == 'vor': data = atm.vor.g; x=coord.lon; y=coord.lat; clev = np.linspace(-1,1,41)*1.0e-6
        if varstrs[ivar] == 'div': data = atm.div.g; x=coord.lon; y=coord.lat; clev = np.linspace(-1,1,41)*1.0e-6
        if varstrs[ivar] == 'str': data = atm.str.g; x=coord.lon; y=coord.lat; clev = np.linspace(-1,1,41)*1.0e6
        if varstrs[ivar] == 'pot': data = atm.pot.g; x=coord.lon; y=coord.lat; clev = np.linspace(-1,1,41)*1.0e6
        if varstrs[ivar] == 'SKE': 
            data = np.zeros([3,N+1]); cstrs=['r','g','b']; x=np.arange(N+1); y = []
            ns = np.linspace(0,N,N+1)
            data[0,:]   = np.sum(atm.vor.s.real **2 + atm.vor.s.imag **2, axis=0)
            data[0,1:] *= rad*rad/4/ns[1:]/(ns[1:]+1) 
            data[1,:]   = np.sum(atm.div.s.real **2 + atm.div.s.imag **2, axis=0)
            data[1,1:] *= rad*rad/4/ns[1:]/(ns[1:]+1) 
            data[2,:] = data[0,:] + data[1,:]
        
        plt.subplot(nrow,ncol,ivar+1)
        if len(y) == 0:
            for i in range(data.shape[0]):
                plt.plot(x,data[i],cstrs[i])

            y1 = 10.0; x1 = 10;  x2 = 200
            plt.plot([x1,x2], [y1, np.exp(-5.0/3* np.log(x2/x1))*y1],'k-')
            plt.plot([x1,x2], [y1, np.exp(-3.0* np.log(x2/x1))*y1],'k-')
            plt.xscale('log')
            plt.yscale('log')
        else:
            plt.contourf(x,y,data,30)
            plt.colorbar(orientation='horizontal')
            plt.contour(x,y,data,[0, 0],colors='k')
            if figcoast[0]:
                fplotcoast(coord.lon, color=figcoast[1])

        plt.title(varstrs[ivar]+': '+tstr)

    plt.tight_layout()
    plt.savefig(ffig)
    plt.close(1)


def fplotcoast(x, color=[0,0,0], lw=1):
    fn = '../topo/coast.txt'
    coast = np.loadtxt(fn)
    lat = coast[:,0]
    lon = coast[:,1]
    if x[0] < 0 :
        plt.plot(lon,lat,color=color,linewidth=lw)
    else:
        lon2 = np.mod(lon,360)
        lon2[abs(lon2[1:]-lon2[:-1]) > 350] = np.nan
        plt.plot(lon2,lat,color=color,linewidth=lw)
    return 


