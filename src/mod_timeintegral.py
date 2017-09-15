import numpy as np
from namelist import *
from mod_commons import *


def timeint(atms, flags, coord, dtime, intflag=0):

    # atms: 2 time layers
    # flags: 1 --- current time step; 0 --- t-dt back time

    ic = flags.index(1); ib = flags.index(0)

    vorsb = atms[ib].vor.s + 0
    divsb = atms[ib].div.s + 0
    gpsb  = atms[ib].gp.s + 0
    
    atms[ic].gp.stog(coord)  
    atms[ic].vor.stog(coord)
    atms[ic].div.stog(coord)
    atms[ic].U.g, atms[ic].V.g = fvordivstouvg(atms[ic].vor.s, atms[ic].div.s, coord)

    A = atms[ic].U.g * (atms[ic].vor.g + coord.f2d.g)
    B = atms[ic].V.g * (atms[ic].vor.g + coord.f2d.g)
    C = atms[ic].U.g * atms[ic].gp.g
    D = atms[ic].V.g * atms[ic].gp.g
    E = Field()
    E.g = (atms[ic].U.g **2 + atms[ic].V.g **2)/2/(coord.clat2d**2)
    E.gtos(coord)

    ns = coord.n2d*(coord.n2d+1)/rad/rad
    Pnm, Onm = fuvgtovordivs(A, B, coord)
    Pnm += (E.s + atms[ic].gps.s)*ns 
    Onm *= -1

    temp, Qnm = fuvgtovordivs(C, D, coord, flag=1)
    Qnm *= -1

    if intflag==0:   # semi-implicit
        Rnm = atms[ib].div.s + 2*dtime*Pnm + atms[ib].gp.s*dtime*ns
        Snm = atms[ib].gp.s  + 2*dtime*Qnm - atms[ib].div.s*atms[ic].gpm*dtime

        atms[ib].vor.s += 2*dtime*Onm

        Gnm = 1 + ns*atms[ic].gpm*(dtime**2)
        atms[ib].div.s = (Rnm + Snm*ns*dtime) / Gnm
        atms[ib].gp.s  = (Snm - Rnm*atms[ic].gpm*dtime) / Gnm

    if intflag==1:  # Euler-forward
        atms[ib].vor.s = atms[ic].vor.s + 2*dtime*Onm
        atms[ib].div.s = atms[ic].div.s + 2*dtime*(Pnm + ns*atms[ic].gp.s)
        atms[ib].gp.s  = atms[ic].gp.s + 2*dtime*(Qnm - atms[ic].gpm*atms[ic].gp.s)
        

    if intflag==2:  # Central
        atms[ib].vor.s += 2*dtime*Onm
        atms[ib].div.s += 2*dtime*(Pnm + ns*atms[ic].gp.s)
        atms[ib].gp.s  += 2*dtime*(Qnm - atms[ic].gpm*atms[ic].gp.s)

    #---------------------------------------
    # horizontal diffusion    
    ns = coord.n2d * (coord.n2d+1)
    if (hdiff[0]):
        atms[ib].vor.s /= 1+2*dtime*K2*ns/(rad**2)
        atms[ib].div.s /= 1+2*dtime*K2*ns/(rad**2)
        atms[ib].gp.s  -= 2*dtime*K2*ns*atms[ib].gps.s
        atms[ib].gp.s  /= 1+2*dtime*K2*ns

    if (hdiff[1]):
        atms[ib].vor.s /= 1+2*dtime*K4*(ns*ns - 4.0)/(rad**4)
        atms[ib].div.s /= 1+2*dtime*K4*(ns*ns - 4.0)/(rad**4)
        atms[ib].gp.s  -= 2*dtime*K4* ns*ns*atms[ib].gps.s /(rad**4) 
        atms[ib].gp.s  /= 1+2*dtime*K4*ns*ns/(rad**4)

#        atms[ib].vor.s /= 1+2*dtime*K4*ns*ns
#        atms[ib].div.s /= 1+2*dtime*K4*ns*ns
#        atms[ib].gp.s  -= 2*dtime*K4*ns*ns*atms[ib].gps.s
#        atms[ib].gp.s  /= 1+2*dtime*K4*ns*ns

    #---------------------------------------
    if runcase==4:
        atms[ib].vor.s = atms[ic].vor.s + 0
        atms[ib].div.s = atms[ic].div.s + 0

    flags[0] = 1-flags[0] 
    flags[1] = 1-flags[1]

    # time filter
    if timefilter > 0:
        ic = flags.index(1); ib = flags.index(0)
        atms[ib].vor.s = ftimefilter(vorsb, atms[ib].vor.s, atms[ic].vor.s, alpha = timefilter)
        atms[ib].div.s = ftimefilter(divsb, atms[ib].div.s, atms[ic].div.s, alpha = timefilter)
        atms[ib].gp.s  = ftimefilter(gpsb,  atms[ib].gp.s,  atms[ic].gp.s,  alpha = timefilter)

    return

def ftimefilter(xb, xc0, xf, alpha = 0.01):
    xc = xc0 + alpha*(xb - 2*xc0 + xf)
    return xc


