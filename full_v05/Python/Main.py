# -*- coding: utf-8 -*-
"""
@author: Levi K C Fisher

Comments:
This is the main file of the Oscillator Add-on used to run FEL simulations. The
parameters in use have been set for the UH Manoa Accelerator group. To modify
simulation parameters the associated variable can either be: changed in the
script below, or set while calling this script from the command-line. Variables
changed frequently have been marked to the edge of the script.
"""

# Built-in Modules used
import os

# External Modules used
import numpy as np

# Custom Modules used
import cavity_utils as cav
from ginger3d_utils import run_g3d



# General constants
c = 299792458



##################################################
#               Oscillator Params
##################################################
# Cavity params
CavityLength = 2.0469                   # length of full cavity in meters
MirrorRadiusCurv = 1.3                  # lens (mirror) radius of curvature
Pwr_C = 0.93                            # cavity reflectivity


# Simulation run params #######################################################

# number of passes through oscillator
npass = 400

# cavity desynchronization
cav_d = 0
#cav_d = [0, 0.005, 0.01]



##################################################
#               Beam Params
##################################################
# Energy distribution
energy = 40                             # energy in MeV
eSpread = 0.5e-2                        # fractional energy spread
emitN = 8e-6                            # normalized transverse emittance

# Twiss params
betaX_twiss = 1.4
betaY_twiss = 0.24
alphaX_twiss= 0.4714
alphaY_twiss= 0

# Peak of the current in Amps
currentMax = 30


# Pulse manipulation params ###################################################

# Electron-bunch length in seconds (aka the pulse full-width half-max) 
bunchDuration = 2e-12

# Electron-bunch RMS                            # fwhm = rms * 2*(2*ln(2))^0.5
bunchRMS = bunchDuration / np.sqrt( np.log(256) )

# Electron time-of-flight Jitter
jitterRMS = 0

# Electron desync (ramping)
d_ramp = np.zeros(npass)



##################################################
#               FEL Params
##################################################
# Undulator params
unduPeriod = 0.023                      # undulator period/wavelength in meters
unduK = 1.2/np.sqrt(2)                  # rms of undulator parameter K
numPeriod = 47                          # number of periods in oscillator

unduL = unduPeriod * numPeriod          # length of undulator in meters

# Radiation params
radWavelength = 3.2281e-6               # radiation wavelength in meters
fcar = c/radWavelength                  # radiation carrier frequency in Hz
slipp = radWavelength*numPeriod         # slippage length in meters



##################################################
#        Ginger3D Simulation Params
##################################################
# Time slice specifications
ebunch_slices = 192
rad_slices = 256
#ebunch_slices = 96
#rad_slices = 192
intrct_slices = 15

# 2D spatial grid
r_grid = 0.08e-3
n_grid = 129



##################################################
#               Create Dictionary
##################################################
# Define dictionary to pass as argument for file-generating functions
Input_Dict = {
    'energy'  : energy,
    'eSpread' : eSpread,
    'emitx'   : emitN,
    'emity'   : emitN,
    'currentMax' : currentMax,
    'betaX'  : betaX_twiss,
    'betaY'  : betaY_twiss,
    'alphaX' : alphaX_twiss,
    'alphaY' : alphaY_twiss,
    'RMSwidth' : bunchRMS,
    
    'unduPeriod' : unduPeriod,
    'unduLength' : unduL,
    'RMSunduK'   : unduK,
    'radWavlen'  : radWavelength,
    
    'eslices'   : ebunch_slices,
    'rslices'   : rad_slices,
    'intslices' : intrct_slices,
    
    'rgrid' : r_grid,
    'ngrid' : n_grid
    }



##################################################
#           Initialize run_g3d Object
##################################################
# Define object from dictionary
runG3Dsim = run_g3d(Input_Dict)

# Set user directories ########################################################
my_dir   = '/Users/levikf/full_v05/Python/'
xg3d_path = '/Users/amirweinberg/Documents/ginger/3d/xg3D-static-1.1.0-a5'

runG3Dsim.home_dir     = my_dir
runG3Dsim.file_indir   = '/Users/levikf/full_v05/in_dir/'
runG3Dsim.file_outdir  = '/Users/levikf/full_v05/out_dir/'
runG3Dsim.g3d_exe_path = xg3d_path



##################################################
#              Cavity: First Pass
##################################################
# Run Ginger-3D by calling input_EE with 1
ar1, ai1, simParams = runG3Dsim.input_EE(1)

# Extract simulation parameters
delt, i0, dx1, rad_s, rad_e, ebm_s, ebm_e = simParams

# Define desynchronization values, in cavity shift length and time
desyn = cav_d
desyn_shift = desyn * slipp
dshifts_s = desyn_shift/c


## Define axes ##
# time-frame of radiation field exiting undulator
radT_ax = np.arange(rad_e - delt/2, rad_s - delt/2, -delt)

# time-frame of field entering undulator (post-cavity)
myT_ax = radT_ax + dshifts_s

# frequency
radF_ax = np.arange(len(radT_ax)) / (radT_ax[0] - radT_ax[-1])
radF_ax += fcar - np.median(radF_ax)    # shift central frequency to carrier


# Implement radiation-field into cavity
#arn1, ain1, waist_S = \
arn1, ain1 = \
                cav.inser_lens(ar1, ai1, radWavelength, CavityLength, unduL, \
                               MirrorRadiusCurv/2, dx1, desyn_shift)

# Implement desynchronization
ars1, ais1 = cav.desynchronize(arn1, ain1, radT_ax, myT_ax)

# Cavity loss
arL = ars1 * Pwr_C**(0.5)
aiL = ais1 * Pwr_C**(0.5)

# Calculate radiation pulse entering undulator
pwr_f, pwr_t, phs = cav.find_pwr(arL, aiL, i0, dx1)

# Display peak power
print(f"Pulse at pass {1}", f"Pulse peak power={np.max(pwr_t)} [W]", sep='\n')



##################################################
#           Cavity: Second Pass and Beyond
##################################################
# Set up pass array
rem_pass = np.arange(npass - 1) + 2

# Run all other passes
for pnum in rem_pass:
    # Get shift for current pass
    my_shft = jitterRMS*np.random.randn() + d_ramp[pnum-1]*slipp/c
    new_time = radT_ax + my_shft
    
    # Shift radiation to current frame
    ars, ais = cav.desynchronize(arL, aiL, radT_ax, new_time)

    # Through Ginger3D program
    #arn, ain, simParams = runG3Dsim.input_EE(pnum, ar=ar1, ai=ai1)
    ar, ai, sPar = runG3Dsim.input_EE(pnum, ar=ars, ai=ais, EFFICIENT=False)
    
    # Reverse timing-offset shift
    arns, ains = cav.desynchronize(ar, ai, new_time, radT_ax)
    
    # Through cavity
#    arn, ain, waist_S = \
    arn, ain = \
                cav.inser_lens(arns, ains, radWavelength, CavityLength, unduL,\
                               MirrorRadiusCurv/2, dx1, desyn_shift)
                    
    # Shift field to undulator entrance
    arnp, ainp = cav.desynchronize(arn, ain, radT_ax, myT_ax)
    
    # Cavity loss
    arL = arnp * Pwr_C**(0.5)
    aiL = ainp * Pwr_C**(0.5)
    
    # Calculate field entering undulator
    pwr_f, pwr_t, phs = cav.find_pwr(arL, aiL, i0, dx1)
    
    # Display peak
    print(f"Pulse at pass {pnum}\n")
    print(f"Pulse peak power={np.max(pwr_t)} [W]")




print('Simulation Complete!')

