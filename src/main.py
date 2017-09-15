#===============================================================================
import numpy as np
import copy
from time import time as timer

from os import system
from namelist import *
from mod_initial import *
from mod_commons import *
from mod_timeintegral import *
from mod_io import *
#===============================================================================
casedir = '../cases/'+casename+'/'
system('mkdir '+casedir)

coord = Coordinate()            # coordinate
atms  = ( Atmosphere(), Atmosphere() )   # two wind fields for time integration
flags = [1,0]              # 1 current time,  0 previous time step
initatm(atms[flags.index(1)], coord)  # initialize atmosphere
initatm(atms[flags.index(0)], coord)  # initialize atmosphere

nstep =  timeend
if timeend < 0:  nstep = int(-timeend*tunit/timestep)

fncout = NCfile()
time = [0,0,0,0]  # day, hour, min, second

for istep in range(nstep+1):

    t0 = timer()

    sec = int(istep*timestep)
    time[0] = sec//86400
    time[1] = np.mod(sec,86400)//3600
    time[2] = np.mod(sec,3600)//60
    time[3] = np.mod(sec,60)

    tstr = 'd'+format(time[0],'04')+'_'+format(time[1],'02')+'-'+format(time[2],'02')+'-'+format(time[3],'02')

    #------------------------------
    if np.mod(istep,nstepfoutvar) == 0:
        fncout.addatm(atms[flags.index(1)],coord,sec)

    if np.mod(istep,nstepfoutfile) == 0:
        fnout = casedir+casename+'_'+tstr+'.nc'
        fncout.writefile(fnout,coord)
    
    if np.mod(istep,nstepfig) == 0:
        ffig = casedir+'fig_'+casename+'_'+tstr+'.'+figformat
        fplotatm(ffig,tstr,atms[flags.index(1)],coord)
    #------------------------------

    # time integration
    if istep == 0:
        timeint(atms, flags, coord, timestep/2, intflag=0)
    else:
        timeint(atms, flags, coord, timestep, intflag=0)

    t1 = timer()
    print 'step: ',istep,', ',tstr,', cputime=',t1-t0

#================================= End =========================================
