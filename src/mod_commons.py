#   use numpy.polynomial.legendre.leggauss instead of scipy.special.legendre, 
# to be able to generate high-resolution Gaussian grids

import numpy as np

#   we have to calculate the Normalized Associated Legendre Functions and 
# their derivation by ourselves, because the default function in python: 
# scipy.special.lpmn returns only the Associated Legendre Function, which is 
# not normalized; when n,m become large, it will be numerically unstable 
# (return inf or nan). 
# from scipy.special import lpmn  

# scipy.fftpack.fft is faster than numpy.fft

from scipy.fftpack import fft, ifft    

from namelist import *


#   Return the Normalized Associated Legendre Functions
# ----------------------------------------------------------------------------
def nalf(M,N,lat):
    slat = np.sin(lat*np.pi/180.0)
    clat = np.cos(lat*np.pi/180.0)
    lp = np.zeros([M+1, N+2])
    hp = np.zeros([M+1, N+1])
    
    lp[0,0] = 1.0/np.sqrt(2.0)
    for m in range(1,M+1):
        lp[m,m] = -clat*np.sqrt((2.0*m+1)/2.0/m) * lp[m-1,m-1]
        lp[m-1,m] = np.sqrt(2.0*m+1)*slat * lp[m-1,m-1]
        hp[m,m] = -m*slat*lp[m,m]
    
    for m in range(M+1):
        for n in range(m+1,N+1):
            fac  = np.sqrt((1.0*n*n-m*m)/(4.0*n*n-1.0))
            fac1 = np.sqrt(((n+1.0)*(n+1.0)-m*m)/(4.0*(n+1.0)*(n+1.0)-1.0))
            lp[m,n+1] = (slat*lp[m,n] - fac*lp[m,n-1]) / fac1
            hp[m,n] = (2.0*n+1)*fac*lp[m,n-1] - n*slat*lp[m,n]

    return lp[:,:-1], hp


#   ucos, vcos in grid points -> vorticity, divergence in spectral space
# ----------------------------------------------------------------------------
def fuvgtovordivs(ug,vg,coord,flag=0):
    (vor, div) = ( Field(), Field() )
    ufft = fft(ug, axis=1)/nlon # ucos, [nlat,nlon]
    vfft = fft(vg, axis=1)/nlon # vcos, [nlat,nlon]
    temp = coord.wt/rad/(coord.coslat**2)
    for m in range(M+1):
        div.s[m, m:N+1]  = np.dot( coord.lp[m, m:N+1,:], ufft[:,m]*1j*m*temp )
        div.s[m, m:N+1] -= np.dot( coord.hp[m, m:N+1,:], vfft[:,m]*temp )
    if flag==0:
        for m in range(M+1):
            vor.s[m, m:N+1]  = np.dot( coord.lp[m, m:N+1,:], vfft[:,m]*1j*m*temp )
            vor.s[m, m:N+1] += np.dot( coord.hp[m, m:N+1,:], ufft[:,m]*temp )
    return vor.s, div.s


#   vorticity, divergence in spectral space -> ucos, vcos in grid points
# ----------------------------------------------------------------------------
def fvordivstouvg(vors,divs,coord):
    (u,v) = ( Field(), Field() )  # ucos, vcos
    ufft  = np.zeros([nlat,nlon]) + 0j
    vfft  = np.zeros([nlat,nlon]) + 0j

    n0   = np.arange(1,N+1)
    temp = rad/n0/(n0+1)
    ufft[:,0] =  np.dot(vors[0,1:N+1]*temp, coord.hp[0,1:N+1,:])
    vfft[:,0] = -np.dot(divs[0,1:N+1]*temp, coord.hp[0,1:N+1,:])
    for m in range(1,M+1):
        ns   = np.arange(m,N+1,1.0)
        temp = rad/ns/(ns+1)
        ufft[:, m]  = -np.dot(divs[m,m:N+1]*1j*m*temp, coord.lp[m,m:N+1,:])
        ufft[:, m] +=  np.dot(vors[m,m:N+1]*temp,      coord.hp[m,m:N+1,:])
        ufft[:,-m]  =  np.conj(ufft[:,m])
        vfft[:, m]  = -np.dot(vors[m,m:N+1]*1j*m*temp, coord.lp[m,m:N+1,:])
        vfft[:, m] -=  np.dot(divs[m,m:N+1]*temp,      coord.hp[m,m:N+1,:])
        vfft[:,-m]  =  np.conj(vfft[:,m])
    u.g = ifft(ufft*nlon, axis=1).real  # ucos
    v.g = ifft(vfft*nlon, axis=1).real  # vcos
    return u.g, v.g


#   vorticity, divergence in spectral space -> stream, potential function 
# ----------------------------------------------------------------------------
def fvordivstostrpots(ssx):
    ssy = ssx * 0
    for n in range(1,N+1):
        mmax = min(M, n) + 1
        ssy[0:mmax,n] = -rad*rad* ssx[0:mmax,n] /n/(n+1.0)
    return ssy


#   stream, potential function in spectral space -> vorticity, divergence
# ----------------------------------------------------------------------------
def fstrpotstovordivs(ssx):
    ssy = ssx * 0
    for n in range(N+1):
        mmax = min(M, n) + 1
        ssy[0:mmax,n] = -n*(n+1.0)*ssx[0:mmax,n]/rad/rad
    return ssy



#   Coordinate and Legendre polynomial
# ----------------------------------------------------------------------------
class Coordinate:
    def __init__(self):
        self.lon  = np.linspace(0.0, 360.0-360.0/nlon, nlon)

        x = np.polynomial.legendre.leggauss(nlat)  # a tuple contains: sin(lat) and weights
        self.lat    = np.arcsin(x[0])*180.0/np.pi
        self.sinlat = x[0] + 0
        self.coslat = np.cos(self.lat*np.pi/180.0)
        self.wt     = x[1] + 0
        self.clat2d = np.zeros([nlat,nlon])
        self.f2d    = Field()
        for ilat in range(nlat):
            self.clat2d[ilat,:] = self.coslat[ilat]
            self.f2d.g[ilat,:] = 2*omg*self.sinlat[ilat]
        self.f2d.s[0,1] = omg/np.sqrt(0.375)
        
        self.m2d  = np.zeros([M+1,N+1])
        self.n2d  = np.zeros([M+1,N+1])
        self.lp   = np.zeros([M+1,N+1,nlat])
        self.hp   = np.zeros([M+1,N+1,nlat])

        for m in range(M+1):
            if trun[0] == 'T': nmax = N
            if trun[0] == 'R': nmax = min(N,m+M)
            self.m2d[m,m:] = m
            self.n2d[m,m:nmax+1] = np.arange(m,nmax+1)
        
        for ilat in range(nlat):
            self.lp[:,:,ilat], self.hp[:,:,ilat] = nalf(M,N,self.lat[ilat])

        
#   2d Fields: grid points, spherical harmonic coefficients
# ----------------------------------------------------------------------------
class Field:
    def __init__(self):
        self.g   = np.zeros([nlat,nlon])        # Grid point field [lat, lon]
        self.s   = np.zeros([M+1,N+1]) + 0j     # Spherical harmonic field [m, n]

    def gtos(self,coord):
        temp  = fft(self.g, axis=1)/nlon  # [nlat,nlon] ------- fastest 
        for m in range(M+1):
            self.s[m, m:N+1] = np.dot( coord.lp[m, m:N+1,:], temp[:,m]*coord.wt )

    def stog(self,coord):
        temp = np.zeros([nlat,nlon]) + 0j
        temp[:,0] = np.dot(self.s[0,0:N+1], coord.lp[0,0:N+1,:])
        for m in range(1,M+1):
            temp[:, m] = np.dot(self.s[m,m:N+1], coord.lp[m,m:N+1,:])
            temp[:,-m] = np.conj(temp[:,m])
        self.g = ifft(temp*nlon, axis=1).real


#   shallow water model fields
# ----------------------------------------------------------------------------
class Atmosphere:
    def __init__(self):
        self.gp   = Field()       # perturbation geopotential (g*h')
        self.gps  = Field()       # surface geopotential (topography) (g*hs)
        self.gpm  = 0.0           # mean geopotential (g*hm)
        self.U   = Field()        # u*coslat
        self.V   = Field()        # v*coslat
        self.vor = Field()        # relative vorticity
        self.div = Field()        # divergence
        self.str = Field()        # stream function
        self.pot = Field()        # potential function

    def sstr_to_svor(self):
        self.vor.s = fstrpotstovordivs(self.str.s)

    def spot_to_sdiv(self):
        self.div.s = fstrpotstovordivs(self.pot.s)

    def svor_to_sstr(self):
        self.str.s = fvordivstostrpots(self.vor.s)

    def sdiv_to_spot(self):
        self.pot.s = fvordivstostrpots(self.div.s)


#===============================================================================

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    coord = Coordinate()
    #for n in range(N+1):
    #    print n, np.sum(coord.fac[n,:])
    #atm = Atmosphere()

    #print coord.f2d.g.shape, coord.fac.shape, coord.lp.shape, coord.hp.shape



