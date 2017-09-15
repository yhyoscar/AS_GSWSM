import numpy as np
from netCDF4 import Dataset as netcdf
from namelist import *
from mod_commons import *

#   initialization
# ----------------------------------------------------------------------------
def initatm(atm, coord, flag=runcase):
    
    # RH4
    if flag == 1:
        w = K = 3.924e-6;   R = 4;  h0 = 8000.0
        for ilat in range(nlat):
            atm.str.g[ilat,:] = rad*rad*coord.sinlat[ilat] * (-w + K*(coord.coslat[ilat]**R)*np.cos(R*coord.lon*np.pi/180.0))

            Alat = 0.5*w*(2*omg+w)*(coord.coslat[ilat]**2) + 0.25*K*K*(coord.coslat[ilat]**(2*R))*((R+1)*(coord.coslat[ilat]**2) + (2*R*R-R-2) - 2*R*R/(coord.coslat[ilat]**2))
            Blat = 2*(omg+w)*K*(coord.coslat[ilat]**R)*((R*R+2*R+2) - (R+1)*(R+1)*(coord.coslat[ilat]**2))/(R+1)/(R+2)
            Clat = 0.25*K*K*(coord.coslat[ilat]**(2*R))*((R+1)*(coord.coslat[ilat]**2) - (R+2))
            atm.gp.g[ilat,:] = grav*h0 + rad*rad*(Alat + Blat*np.cos(R*coord.lon*np.pi/180) + Clat*np.cos(2*R*coord.lon*np.pi/180))

        atm.gpm = np.sum(np.mean(atm.gp.g, axis=1)*coord.wt)/2.0
        atm.gp.g -= atm.gpm
        atm.gp.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()

    # Horta and Simmons
    if flag == 2:
        for ilat in range(nlat):
            for ilon in range(nlon):
                latp = np.arcsin(max(min(coord.coslat[ilat]*np.cos(coord.lon[ilon]*np.pi/180.0),1),-1))*180.0/np.pi
                lonp = np.arccos(max(min(-coord.sinlat[ilat]/np.cos(latp*np.pi/180.0), 1), -1))*180.0/np.pi
                atm.vor.g[ilat,ilon] = omg*(coord.sinlat[ilat]/22.0 + 3*((9.0/8.0)**4)*(np.cos(latp*np.pi/180.0)**8) * np.sin(latp*np.pi/180.0)*np.cos(8*lonp*np.pi/180.0))
        atm.vor.gtos(coord)

    if flag==3:
        atm.vor.s[0,0:2] = 1.0e-4
#        atm.vor.s[1,1] = 1.0e-4

    # advection of cosine bell over the pole
    if flag==4:
        h0 = 1000.0;  R=rad/3.0;  u0=2*np.pi*rad/(12.0*86400)
        lonc = 270.0*np.pi/180.0;   latc = 0.0;   alpha = 0.05 + np.pi/2
        
        for i in range(nlat):
            for j in range(nlon):
                r = min(rad*np.arccos(np.sin(latc)*coord.sinlat[i] + np.cos(latc)*coord.coslat[i]*np.cos(coord.lon[j]*np.pi/180-lonc)), R)
                atm.gp.g[i,j] = 0.5*h0*(1+np.cos(np.pi*r/R))
                atm.str.g[i,j] = -rad*u0*(coord.sinlat[i]*np.cos(alpha) - np.cos(coord.lon[j]*np.pi/180.0)*coord.coslat[i]*np.sin(alpha))

        atm.gp.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()

    # steady state nonlinear zonal geostrophic flow
    if flag==5:
        alpha = np.pi/2;  u0 = 2*np.pi*rad/(12.0*86400);  gh0 = 2.94e4
        for i in range(nlat):
            atm.str.g[i,:] = -rad*u0*(coord.sinlat[i]*np.cos(alpha) - coord.lon*coord.coslat[i]*np.sin(alpha))
            temp = (-np.cos(coord.lon*np.pi/180)*coord.coslat[i]*np.sin(alpha) + coord.sinlat[i]*np.cos(alpha))
            coord.f2d.g[i,:] = 2*omg*temp
            atm.gp.g[i,:] = gh0 - (rad*omg*u0+u0*u0/2)*temp*temp

        coord.f2d.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()
        atm.gp.gtos(coord)

    # flow over mountain
    if flag==6:
        alpha = 0.0;  u0 = 20;  gh0 = grav*5400
        R = np.pi/9;  lonc = 90.0 *np.pi/180.0;  latc = 30 *np.pi/180.0;
        hs0 = 2000.0
        for i in range(nlat):
            atm.str.g[i,:] = -rad*u0*(coord.sinlat[i]*np.cos(alpha) - coord.lon*coord.coslat[i]*np.sin(alpha))
            temp = (-np.cos(coord.lon*np.pi/180)*coord.coslat[i]*np.sin(alpha) + coord.sinlat[i]*np.cos(alpha))
            atm.gp.g[i,:] = gh0 - (rad*omg*u0+u0*u0/2)*temp*temp
            for j in range(nlon):
                r = np.sqrt(min(R*R, (coord.lon[j]*np.pi/180-lonc)**2 + (coord.lat[i]*np.pi/190-latc)**2))
                atm.gps.g[i,j] = grav*hs0 * (1.0-r/R)

        atm.gps.gtos(coord)
        atm.gpm = np.sum(np.mean(atm.gp.g, axis=1)*coord.wt)/2.0
        atm.gp.g -= atm.gpm
        atm.gp.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()

    # baroclinic wave
    if flag==9:
        alpha = 0.0;  u0 = 20;  gh0 = grav*5400;  up=1.0
        R = np.pi/9;  lonc = 90.0 *np.pi/180.0;  latc = 30 *np.pi/180.0;
        hs0 = 2000.0
        for i in range(nlat):
            atm.str.g[i,:] = -rad*u0*(coord.sinlat[i]*np.cos(alpha) - coord.lon*coord.coslat[i]*np.sin(alpha))
            temp = (-np.cos(coord.lon*np.pi/180)*coord.coslat[i]*np.sin(alpha) + coord.sinlat[i]*np.cos(alpha))
            atm.gp.g[i,:] = gh0 - (rad*omg*u0+u0*u0/2)*temp*temp

        atm.gps.gtos(coord)
        atm.gpm = np.sum(np.mean(atm.gp.g, axis=1)*coord.wt)/2.0
        atm.gp.g -= atm.gpm
        atm.gp.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()

        atm.U.g, atm.V.g = fvordivstouvg(atm.vor.s,atm.div.s,coord)
        for i in range(nlat):
            for j in range(nlon):
                r = np.sqrt(min(R*R, (coord.lon[j]*np.pi/180-lonc)**2 + (coord.lat[i]*np.pi/190-latc)**2))
                atm.U.g[i,j] += up * (1.0-r/R) * coord.coslat[i]
        atm.vor.s, atm.div.s = fuvgtovordivs(atm.U.g,atm.V.g,coord,flag=0)

        
    # equatorial waves test
    if flag==7:
        gh0 = grav*5400;  ghp = grav*10.0;  up=1.0
        R = np.pi/18;  lonc = 90.0 *np.pi/180.0;  latc = 0 *np.pi/180.0;

        atm.gpm  = gh0
        for i in range(nlat):
            for j in range(nlon):
                r = np.sqrt(min(R*R, (coord.lon[j]*np.pi/180-lonc)**2 + (coord.lat[i]*np.pi/190-latc)**2))
                #atm.gp.g[i,j] = ghp * (1.0-r/R)
                atm.U.g[i,j] = up * (1.0-r/R) * coord.coslat[i]

        atm.vor.s, atm.div.s = fuvgtovordivs(atm.U.g,atm.V.g,coord,flag=0)

        
    # realistic topography
    if flag==8:
        fin = '../topo/topo_10min.nc'
        topo = fncread(fin,'htopo')
        lon0 = fncread(fin,'lon')
        lat0 = fncread(fin,'lat')
        atm.gps.g = grav*my_interp_2d(lon0,lat0,topo,coord.lon,coord.lat)
        for i in range(1):
            atm.gps.g = fsmooth2d(atm.gps.g)
        atm.gps.gtos(coord)
        #atm.gps.s *= 1-(coord.n2d/1.4/trun[1])**2
        
        alpha = 0.0;  u0 = 20;  gh0 = grav*20000
        for i in range(nlat):
            atm.str.g[i,:] = -rad*u0*(coord.sinlat[i]*np.cos(alpha) - coord.lon*coord.coslat[i]*np.sin(alpha))
            temp = (-np.cos(coord.lon*np.pi/180)*coord.coslat[i]*np.sin(alpha) + coord.sinlat[i]*np.cos(alpha))
            atm.gp.g[i,:] = gh0 - (rad*omg*u0+u0*u0/2)*temp*temp

        atm.gp.g -= atm.gps.g
        atm.gpm   = np.sum(np.mean(atm.gp.g, axis=1)*coord.wt)/2.0
        atm.gp.g -= atm.gpm
        atm.gp.gtos(coord)
        atm.str.gtos(coord)
        atm.sstr_to_svor()

         
#        atm.gp.g = gh0 + ghp*np.random.rand(nlat, nlon) - atm.gps.g
#        atm.gpm = np.sum(np.mean(atm.gp.g, axis=1)*coord.wt)/2.0
#        atm.gp.g -= atm.gpm
#        atm.gp.gtos(coord)
 
#        atm.U.g = -fdiff_lat(atm.gp.g, coord.lat)/rad/coord.f2d.g * coord.clat2d
#        atm.V.g =  fdiff_lon(atm.gp.g, coord.lon)/rad/coord.f2d.g
#        atm.U.g = np.min([atm.U.g*0 + 50, np.max([atm.U.g*0-50, atm.U.g],axis=0)], axis=0)
#        atm.V.g = np.min([atm.V.g*0 + 50, np.max([atm.V.g*0-50, atm.V.g],axis=0)], axis=0)

#        atm.vor.s, atm.div.s = fuvgtovordivs(atm.U.g, atm.V.g,coord)
        
        
def my_interp_2d(lon1,lat1,data1,lon2,lat2):
    nlon1 = np.size(lon1)
    nlat1 = np.size(lat1)
    nlon2 = np.size(lon2)
    nlat2 = np.size(lat2)

    temp = np.zeros([nlat1,nlon2])
    for ilat in range(nlat1):
        temp[ilat,:] = np.interp(lon2,lon1,data1[ilat,:])

    data2 = np.zeros([nlat2,nlon2])
    for ilon in range(nlon2):
        data2[:,ilon] = np.interp(lat2,lat1,temp[:,ilon])
    return data2        

        
def fncread(fn,var):
    fid = netcdf(fn, 'r')
    data = fid.variables[var][:]
    fid.close()
    return data


def fdiff_lat(x,lat0):
    nlat = np.size(lat0)
    lat = lat0*np.pi/180.0
    dx = x * 0.0
    for ilat in range(nlat-2):
        dx[ilat+1,:] = (x[ilat+2,:] - x[ilat,:]) / (lat[ilat+2] - lat[ilat])
    dx[0,:]      = (x[1,:] - x[0,:]) / (lat[1] - lat[0])
    dx[nlat-1,:] = (x[nlat-1,:] - x[nlat-2,:]) / (lat[nlat-1] - lat[nlat-2])
    return dx


def fdiff_lon(x,lon0):
    nlon = np.size(lon0)
    lon = lon0*np.pi/180.0
    dx = x * 0.0
    for ilon in range(nlon-2):
        dx[:,ilon+1] = (x[:,ilon+2] - x[:,ilon]) / (lon[ilon+2] - lon[ilon])
    dx[:,0]      = (x[:,1] - x[:,nlon-1]) / (lon[1] - lon[nlon-1])
    dx[:,nlon-1] = (x[:,0] - x[:,nlon-2]) / (lon[0] - lon[nlon-2])
    return dx        
    
def fsmooth2d(x):
    y = x+0;   wt = 0.5
    y[1:-1,1:-1] = wt*x[1:-1,1:-1] + (1-wt)*0.25*(x[0:-2,1:-1] + x[2:,1:-1] + x[1:-1,0:-2] + x[1:-1,2:])
    y[1:-1,0]  = wt* x[1:-1,0] + (1-wt)*0.25*(x[1:-1,-1] + x[1:-1,1] + x[0:-2,0] + x[2:,0])
    y[1:-1,-1] = wt*x[1:-1,-1] + (1-wt)*0.25*(x[1:-1,-2] + x[1:-1,0] + x[0:-2,-1] + x[2:,-1])
    y[0,1:-1]  = 0.5*x[0,1:-1] + 0.25*(x[0,2:] + x[0,0:-2])
    y[0,0]     = 0.5*x[0,0]  + 0.25*(x[0,2] + x[0,-1])
    y[0,-1]    = 0.5*x[0,-1] + 0.25*(x[0,0] + x[0,-2])
    y[-1,1:-1]  = 0.5*x[-1,1:-1] + 0.25*(x[-1,2:] + x[-1,0:-2])
    y[-1,0]     = 0.5*x[-1,0]  + 0.25*(x[-1,2] + x[-1,-1])
    y[-1,-1]    = 0.5*x[-1,-1] + 0.25*(x[-1,0] + x[-1,-2])    
    return y


    
