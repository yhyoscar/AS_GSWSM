#   Case name
# -----------------------------------------------------------------------------
casename = 'realtopo_T106_nohd'


#   Grid information
# -----------------------------------------------------------------------------
trun = ['T',106]     # Truncation method and number: 'T' -- triangular; 'R' -- Rhomboidal


#   Time integration information
# -----------------------------------------------------------------------------
tunit    = 86400     
timeend  = -1000       # Ending time: >=0 -- time steps; <0 -- per tunit second
runcase  = 8         # Run case
hdiff    = [False, False]   # K2, K4
timefilter = 0.01


#   Output NetCDF files 
# -----------------------------------------------------------------------------
# Output frequency: 0 -- no ouput; 
#                   + -- time steps; 
#                   - -- per tunit second
# Output variables: str -- stream function (m2/s); 
#                   pot -- potential function (m2/s);
#                   vor -- relative vorticity (1/s);
#                   div -- divergence (1/s);
#                   u,v -- horizonta wind (m/s)
#                   hm  -- mean free surface height (m);
#                   hs  -- topography (m);
#                   h   -- free surface height (m);
foutvarfreq  = -1
foutfilefreq = -30
foutvars = ['hm','hs','h','u','v','sp_vor','sp_div']


#   Output figures
# -----------------------------------------------------------------------------
figfreq   = -10                              # Frequency of plots, similar with Fout
figformat = 'png'                         # Figure format
figvars   = ['hs','h','u','v','SKE']    # Plot variables
figcoast  = [True, 'm']                       # plot coast



#===============================================================================
#                                  constants
#===============================================================================



# Parameters
# ----------------------------------------------------------------------------
rad  = 6371220.0             # Earth's radius    (m)
omg  = 7.292e-5              # Angular speed of Earth's rotation (1/s)
grav = 9.80616               # Gravity 


#   Size of grid, timestep, diffusion coefficient K4, K2
# ----------------------------------------------------------------------------
ntruns = [  15,   21,   30,   42,   63,   85,  106, 170, 213,  340,   680, 1360]
nlons  = [  48,   64,   96,  128,  192,  256,  320, 512, 640, 1024,  2048, 4096]
nlats  = [[ 24,   32,   48,   64,   96,  128,  160, 256, 320,  512,  1024, 2048], 
          [ 38,   54,   76,  108,  158,  214,  266, 426, 534,  854,  1700, 6826]]          
tsteps = [7200, 7200, 4800, 3600, 2400, 1800, 1200, 900, 600,  450,   300,  300]
K4s    = [   0,    0,   0,0.5e16,2.5e15,1.e15,5.0e14,1.5e14,7.e13, 1.5e13, 1.e12, 0]
K2 = 2.5e5


#   Size of grid, timestep, diffusion coefficient K4, K2
# ----------------------------------------------------------------------------
M = trun[1]
if trun[0] == 'T':  N = trun[1]
if trun[0] == 'R':  N = 2*trun[1]

itrun  = ntruns.index(trun[1])
nlon   = nlons[itrun]
if trun[0] == 'T':  nlat = nlats[0][itrun]
if trun[0] == 'R':  nlat = nlats[1][itrun]

timestep = tsteps[itrun]
K4 = K4s[itrun]



