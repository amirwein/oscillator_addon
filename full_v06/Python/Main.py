#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Levi K C Fisher

"""


# Set user directories ########################################################
my_dir   = '/Users/levikf/full_v06/Python/'
xg3d_path = '/Users/amirweinberg/Documents/ginger/3d/xg3D-static-1.1.0-a5'
fin_dir  = '/Users/levikf/full_v06/in_dir/'
fout_dir = '/Users/levikf/full_v06/out_dir/'


# General constants
c = 299792458                           # speed of light
m_e = 0.511                             # eletron mass in MeV


# Oscillator Params
CavityLength = 2.0469                   # length of full cavity in meters
MirrorRadiusCurv = 1.3                  # lens (mirror) radius of curvature
Pwr_C = 0.93                            # cavity reflectivity



def runsim_FELosc_g3d(
        npass=400,
        cav_d=0,
        jitterRMS=0,
        d_ramp=None,
        inp_dict=None,
        CavityLength = CavityLength,
        MirrorRadiusCurv = MirrorRadiusCurv,
        Pwr_C = Pwr_C,
        EFF=False,
        ):
    # Imports
    import numpy as np

    import cavity_utils as cav
    from ginger3d_utils import run_g3d
    
    
    # Assumes no ramping if none given
    if d_ramp is None:
        d_ramp = np.zeros(npass)
        
        
    # Assumes 2ps run
    if inp_dict is None:
        inp_dict = {
            'energy' : 40,
            'eSpread' : 0.5e-2,
            'emitx' : 8e-6,
            'emity' : 8e-6,
            'currentMax' : 30,

            'betaX' : 1.4,
            'betaY' : 0.24,
            'alphaX' : 0.4714,
            'alphaY' : 0,
        
            'RMSwidth' : 2e-12/np.sqrt( np.log(256)),
        
            'unduPeriod' : 0.023,
            'unduLength' : 1.081,
            'RMSunduK' : 1.2/np.sqrt(2),
            'radWavlen' : 3.203e-6,
        
            'eslices' : 192,
            'rslices' : 256,
            'intslices' : 15,
        
            'rgrid' : 0.08e-3,
            'ngrid' : 129
            }
        
    
    # Use specified parameters to derive quantities (slippage)
    radWavelength = inp_dict['radWavlen']
    unduL = inp_dict['unduLength']
    numPeriod = unduL/inp_dict['unduPeriod']
    slipp = radWavelength*numPeriod
    
    desyn = cav_d                       # desyn-units: 2 * del cav_length/slipp
    desyn_shift = desyn * slipp         # desync-induced shift in meters
    dshifts_s = desyn_shift/c           # "" shift in seconds
    
    # Initialize object from dictionary and set directories
    runG3Dsim = run_g3d(inp_dict)
    runG3Dsim.home_dir     = my_dir
    runG3Dsim.file_indir   = fin_dir
    runG3Dsim.file_outdir  = fout_dir
    runG3Dsim.g3d_exe_path = xg3d_path
    
    
    # First Pass + Get Params
    ar1, ai1, simParams = runG3Dsim.input_EE(1)
    delt, i0, dx1, rad_s, rad_e, ebm_s, ebm_e = simParams
    
    #####     Setup axes     #####
    # time-frame of radiation field exiting undulator
    radT_ax = np.arange(rad_e - delt/2, rad_s - delt/2, -delt)
    
    # time-frame of field entering undulator (post-cavity)
    myT_ax = radT_ax + dshifts_s
    
    
    # Implement radiation-field into cavity
    arn1, ain1 = cav.inser_lens(ar1, ai1, radWavelength, CavityLength, unduL, \
                                MirrorRadiusCurv/2, dx1, desyn_shift)
        
    # Implement desynchronization
    ars1, ais1 = cav.desynchronize(arn1, ain1, radT_ax, myT_ax)
    
    # Cavity loss
    arL = ars1 * Pwr_C**(0.5)
    aiL = ais1 * Pwr_C**(0.5)
    
    # Calculate radiation pulse entering undulator
    pwr_f, pwr_t, phs = cav.find_pwr(arL, aiL, i0, dx1)
    
    # Display peak power entering undulator
    print(f"Pulse at pass {1}", 
          f"Pulse peak power={np.max(pwr_t)} [W]", sep='\n')


    # Set up pass array
    rem_pass = np.arange(npass - 1) + 2
    
    # Run Second Pass +
    for pnum in rem_pass:
        # Get shift for current pass
        my_shft = jitterRMS*np.random.randn() + d_ramp[pnum-1]*slipp/c
        new_time = radT_ax + my_shft
    
        # Shift radiation to current (offset) frame
        ars, ais = cav.desynchronize(arL, aiL, radT_ax, new_time)


        # Through Ginger3D program with offset
        ar, ai, _ = runG3Dsim.input_EE(pnum, ar=ars, ai=ais, EFFICIENT=EFF)
        
        # Reverse timing-offset shift
        arns, ains = cav.desynchronize(ar, ai, new_time, radT_ax)

        # Through cavity
        arn, ain = cav.inser_lens(arns, ains, radWavelength, CavityLength, \
                                  unduL, MirrorRadiusCurv/2, dx1, desyn_shift)

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


if __name__ == "__main__":
    runsim_FELosc_g3d()

