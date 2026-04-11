# -*- coding: utf-8 -*-
"""
@author: Levi K C Fisher

Comments:
ginger3d_utils is a Python module used by the UH Manoa Accelerator group made
as part of an Oscillator Add-On package to be used in FEL simulations,
specifically for the Ginger-3D program.
The ginger3d_utils module defines a class structure run_g3d that is used to run
simulations with the Ginger-3D software.
The ginger3d_utils module also defines some functions that extract the data and
parameters used in a Ginger-3D run.
"""

# Built-in Modules used
import os
import glob
import shutil
import subprocess

# Ext. Modules used
import numpy as np
import h5py



# Constants
c = 299792458



# Run Ginger3D Simulations
class run_g3d:
    '''
    Defines the methods used to run simulations with the Ginger-3D software.
    Users should be wary that directory paths may not be the same as what has
    been defined below.
    '''
    ### User Directory Information ###
    # Modify the variables defined below or re-define the variable attribute
    # in your local run_g3d object before calling the file-generating methods.
    
    # Current working directory (should contain all modules)
    home_dir     = '/Users/levikf/full_v05/Python/'
    
    # Directory to generate Ginger3D "in" files
    file_indir   = ''
    
    # Directory to store Ginger3D "out" files
    file_outdir  = ''
    
    # Basename used for output files
    out_base = 'UHM_G3D'
    
    # Absolute path to the executable for Ginger3D
    g3d_exe_path = ''
    
    
    ### Other Ginger-3D Parameters ###
    # Change the number of simulated macroparticles here
    numParts = 8192
    
    
    
    def __init__(self, inp_dict):
        '''
        Initialize an object that runs the Ginger-3D program using parameters 
        defined in a dictionary.

        Parameters
        ----------
        inp_dict : Dictionary
            All parameters used to run Ginger-3D. A value must be defined for
            each of the following keys:
                'energy', 'eSpread', 'emitx', 'emity', 'currentMax', 'betaX',
                'betaY', 'alphaX', 'alphaY', 'RMSwidth', 'unduPeriod',
                'unduLength', 'RMSunduK', 'radWavlen', 'eslices', 'rslices',
                'intslices', 'rgrid', 'ngrid'

        '''
        # Electron-beam params
        self.energyMev  = inp_dict['energy']
        self.eSpread    = inp_dict['eSpread']
        self.currentMax = inp_dict['currentMax']
        self.bunchRMS   = inp_dict['RMSwidth']
        
        # Twiss params
        self.emitX  = inp_dict['emitx']
        self.emitY  = inp_dict['emity']
        self.betaX  = inp_dict['betaX']
        self.betaY  = inp_dict['betaY']
        self.alphaX = inp_dict['alphaX']
        self.alphaY = inp_dict['alphaY']
        
        # FEL params
        self.uPeriod = inp_dict['unduPeriod']
        self.uLength = inp_dict['unduLength']
        self.uK_RMS  = inp_dict['RMSunduK']
        
        self.rLambda = inp_dict['radWavlen']

        # Simulation params
        self.ebm_slices = inp_dict['eslices']
        self.rad_slices = inp_dict['rslices']
        self.int_slices = inp_dict['intslices']
        
        self.rgrid = inp_dict['rgrid']
        self.ngrid = inp_dict['ngrid']



    def make_in_file(self, filename):
        '''
        Function to make a Ginger3D "in" file. (1st pass)
        Specifies that Ginger3D will generate a fld file containing the complex
        E-field amplitudes exiting the undulator.

        Parameters
        ----------
        filename : Str or Path
            Path to (or name of) the file being generated.

        '''
        # Ensure correct path handle to new file
        fPath = os.path.abspath(filename)
        
        with open(fPath, 'w') as f:
            f.write("input file for parameters corresponding ")
            f.write("to UHM IR FEL benchmark$$$\n")
            f.write("&in_beams\n")
            
            f.write("!>>> e-beam parameters:\n")
            f.write(f"emitx_mks = {self.emitX:4.1e}\n")
            f.write(f"emity_mks = {self.emitY:4.1e}\n")
            f.write(f"jmg = {2:2d}\n")
            f.write(f"betaX_twiss = {self.betaX:4.4f}\n")
            f.write(f"betaY_twiss = {self.betaY:4.4f}\n")
            f.write(f"alphaX_twiss = {self.alphaX:4.4f}\n")
            f.write(f"alphaY_twiss = {self.alphaY:4.4f}\n")
            f.write(f"energy_MeV = {self.energyMev:4.4f}\n")
            f.write(f"sigmaE_norm = {self.eSpread:4.4f}\n")
            f.write("gam_load = 'gaussian'\n")
            f.write(f"current = {self.currentMax:4.4f}\n")
            f.write("ebeam_pulse_shape='gaussian'\n")
            f.write(f"ebeam_tbody={self.bunchRMS:1.4e}\n")
            f.write(f"ntestp = {self.numParts:4d}\n")
            f.write(f"nfold_sym = {8:4d}\n")
            f.write("lshot = .t\n")
            
            f.write("\n!>>> undulator parameters:\n")
            f.write(f"wavelw = {self.uPeriod:4.4f}\n")
            f.write("llinear = .t\n")
            f.write(f"idesign = {1:2d}\n")
            f.write(f"aw0 = {self.uK_RMS:4.4f}\n")
            f.write(f"rkxkw = {0.00:4.4f}\n")
            
            f.write("\n!>>> input laser characteristics:\n")
            f.write(f"wavels = {self.rLambda:1.4e}\n")
            f.write(f"plaser = {0:1.1e}\n")
            f.write("/END\n")
            
            f.write("&in_sim_diag\n")
            f.write("l_write_fld=.t\n")
            f.write(f"nside = {self.ebm_slices:4d}\n")
            f.write(f"nslice_sim = {self.rad_slices:4d}\n")
            f.write("lshort=.t\n")
            f.write(f"nsidep = {self.int_slices:4d}\n")
            f.write(" iseeds(1:2)=20241810	12418367\n")
            
            f.write("\n!>> simulation and diagnostic parameters:\n")
            f.write("lfred = .f\n")
            f.write("l_2Dxy=.t\n")
            f.write(f"zmxmeter = {self.uLength:4.4f}\n")
            f.write(f"dzstep = {self.uPeriod:4.4f}\n")
            f.write(f"!rmaxsim = {0.025:4.4f}\n")
            f.write(f"dxg = {self.rgrid:1.4e}\n")
            f.write(f"nxg = {self.ngrid:4d}\n")
            f.write("/END\n")



    def make_2in_file(self, filename):
        '''
        Function to make a Ginger3D "in" file. (2nd+ pass)
        Specifies that Ginger3D generates both a fld file as well a rst file,
        containing the electron data of the bunch.

        Parameters
        ----------
        filename : Str or Path
            Path to (or name of) the file being generated.

        '''
        # Ensure correct path handle to new file
        fPath = os.path.abspath(filename)
        
        with open(fPath, 'w') as f:
            f.write("input file for parameters corresponding ")
            f.write("to UHM IR FEL benchmark$$$\n")
            f.write("&in_beams\n")
            
            f.write("!>>> e-beam parameters:\n")
            f.write(f"emitx_mks = {self.emitX:4.1e}\n")
            f.write(f"emity_mks = {self.emitY:4.1e}\n")
            f.write(f"jmg = {2:2d}\n")
            f.write(f"betaX_twiss = {self.betaX:4.4f}\n")
            f.write(f"betaY_twiss = {self.betaY:4.4f}\n")
            f.write(f"alphaX_twiss = {self.alphaX:4.4f}\n")
            f.write(f"alphaY_twiss = {self.alphaY:4.4f}\n")
            f.write(f"energy_MeV = {self.energyMev:4.4f}\n")
            f.write(f"sigmaE_norm = {self.eSpread:4.4f}\n")
            f.write("gam_load = 'gaussian'\n")
            f.write(f"current = {self.currentMax:4.4f}\n")
            f.write("ebeam_pulse_shape='gaussian'\n")
            f.write(f"ebeam_tbody={self.bunchRMS:1.4e}\n")
            f.write(f"ntestp = {self.numParts:4d}\n")
            f.write(f"nfold_sym = {8:4d}\n")
            f.write("lshot = .t\n")
            
            f.write("\n!>>> undulator parameters:\n")
            f.write(f"wavelw = {self.uPeriod:4.4f}\n")
            f.write("llinear = .t\n")
            f.write(f"idesign = {1:2d}\n")
            f.write(f"aw0 = {self.uK_RMS:4.4f}\n")
            f.write(f"rkxkw = {0.00:4.4f}\n")
            
            f.write("\n!>>> input laser characteristics:\n")
            f.write(f"wavels = {self.rLambda:1.4e}\n")
            f.write(f"plaser = {0:1.1e}\n")
            f.write("/END\n")
            
            f.write("&in_sim_diag\n")
            f.write("l_write_fld=.t\n")
            f.write("l_write_rst=.t\n")
            f.write(f"nside = {self.ebm_slices:4d}\n")
            f.write(f"nslice_sim = {self.rad_slices:4d}\n")
            f.write("lshort=.t\n")
            f.write(f"nsidep = {self.int_slices:4d}\n")
            f.write(" iseeds(1:2)=20241810	12418367\n")
            
            f.write("\n!>> simulation and diagnostic parameters:\n")
            f.write("lfred = .f\n")
            f.write("l_2Dxy=.t\n")
            f.write(f"zmxmeter = {self.uLength:4.4f}\n")
            f.write(f"dzstep = {self.uPeriod:4.4f}\n")
            f.write(f"!rmaxsim = {0.025:4.4f}\n")
            f.write(f"dxg = {self.rgrid:1.4e}\n")
            f.write(f"nxg = {self.ngrid:4d}\n")
            f.write("/END\n")



    def input_EE(self, reflect, ar=0, ai=0, EFFICIENT=False):
        '''
        Function that calls the Ginger3D program to generate "out" files,
        returning the field amplitudes and some of the simulation parameters.

        Parameters
        ----------
        reflect : Int
            Number indicating which pass through undulator will be simulated.
        ar : Array (3D)
            The real components of the E-field entering the undulator. (if 2+)
        ai : Array (3D)
            The imaginary E-field components at the undulator entrance.
        EFFICIENT : Bool, optional
            Toggle for "efficient" storage, which will only save a fraction of
            the Ginger3D output files. The default is False.

        Returns
        -------
        E1 : Array (3D)
            The real components of the E-field exiting the undulator.
        E2 : Array (3D)
            The imaginary components of the E-field at undulator exit.
        Params : List of Floats     # size=(7,)
            List of parameters used in simulation, organized as follows:
                delt:   spacing of time-slices
                i0:     power coefficient
                dx:     grid spacing of x,y steps
                rad_s:  radiation window-start time
                rad_e:  radiation window-end time
                ebm_s:  electron beam-start time
                ebm_e:  electron beam-end time
        '''
        # Localize variables
        fin_dir  = self.file_indir
        fout_dir = self.file_outdir
        out_base = self.out_base
        exe_path = self.g3d_exe_path
        
        # Generate first pass (no input radiation field)
        if reflect == 1:
            # Empty input-dir
            for oldIn in glob.glob(f"{fin_dir}G3d*"):
                os.remove(oldIn)
            
            # Empty the output-directory
            for oldData in glob.glob(f"{fout_dir}*.hdf5"):
                os.remove(oldData)
            
            # Name Ginger3D input file and generate using function
            myinfile = f"{fin_dir}G3din{reflect}"
            self.make_in_file(myinfile)

            # Get command to run Ginger3D
            out_file = f"{out_base}{reflect}"
            myCmd = f"COLUMNS=1000 {exe_path} -i {myinfile} -r {out_file}"

        # Treat subsequent passes differently
        else:
            # Save fewer files
            if EFFICIENT and (reflect % 20):
                reflect = 2

            # Input radiation field
            fld_ifile = f"{fout_dir}tmp_fld.hdf5"
            E3 = np.concatenate((ar, ai), axis=0)
            E3 = np.asfortranarray(E3.T)
            with h5py.File(fld_ifile, "a") as f:
                del f['/radiation/2D_xy_field']
                f.create_dataset('/radiation/2D_xy_field', data=E3)

            # Name Ginger3D input file and generate using function
            myinfile = f"{fin_dir}G3d2in{reflect}"
            self.make_2in_file(myinfile)

            # Get command to run Ginger3D
            out_file = f"{out_base}2{reflect}"
            halfCmd = f"COLUMNS=400 {exe_path} -i {myinfile} -f {fld_ifile}"
            myCmd = halfCmd + f" -r {out_file}"

        # Call the Ginger3D command
        myCommand = myCmd
#        print('My Command : ', myCommand)
        result = subprocess.run(myCommand, capture_output=True, shell=True)

        # Group newly generated files together
        file_suffix = f"{out_file}.hdf5"
        n_files = glob.glob(f"*{file_suffix}")
        for nfile in n_files:
            
            '''# Open-close files for security
            with open(nfile, "r+") as f:
                pass'''
            
            # Move files to out-dir
#            shutil.move(nfile, fout_dir)
            os.rename(nfile, f"{fout_dir}{nfile}")
#            print(f"Moved {nfile} to {fout_dir}")

        # Find files in out-dir and create tmp_fld.hdf5 if DNE
        fld_out_file = f"{fout_dir}fld{file_suffix}"
        g3d_out_file = f"{fout_dir}{file_suffix}"
        fld_file = f"{fout_dir}tmp_fld.hdf5"
        shutil.copy2(fld_out_file, fld_file)
        
        # Read output values from Ginger3D out files
        with h5py.File(g3d_out_file, "r") as f:
            real_params = f["/base_param/real_params"]
            value1 = np.array(real_params["value1"])
            value2 = np.array(real_params["value2"])

        # Get relevant parameters
        delt = value1[28]
        i0 = value1[3]
        dx = value1[23]
        rad_s = value1[6]
        rad_e = value2[6]
        ebm_s = value1[14]
        ebm_e = value2[14]

        Params = [delt, i0, dx, rad_s, rad_e, ebm_s, ebm_e]

        # Read output E-field from Ginger3D fld file
        with h5py.File(fld_out_file, "r") as f:
            fund_field = np.array(f["/radiation/2D_xy_field"]).T
        
        # Separate real/imag components
        half_gb = fund_field.shape[0] // 2
        E1 = fund_field[:half_gb, :, :]
        E2 = fund_field[half_gb:, :, :]
        
        return E1, E2, Params



# Read Ginger3D data
def extract_g3d_par(dir_path, out_base, npass):
    '''
    Extract axes and parameters from Ginger3D output HDF5 files.

    Returns
    -------
    curr_dist : Array (1D)
        Current distribution of the input beam.
    t_ref : Array (1D)
        Reference time-axis for the electron beam.
    frq_ax : Array (1D)
        Frequency-axis for the radiation field.
    time_ax : Array (1D)
        Time-axis for the radiation field.
    runParams : List of Floats      # size=(11,)
        Simulation parameters for runs, in the following order:
            npass:      Number of passes through cavity
            curr_max:   Peak electron-beam current [A]
            ebm_rms:    RMS width of the electron bunch
            i0:         Power (intensity) coefficient [W/m^2]
            dx:         Grid spacing in x,y [m]
            rad_wavlen: Radiation wavelength [m]
            ebm_s:      Electron beam-start time
            ebm_e:      Electron beam-end time
            delt:       Time-slice spacing [s]
            rad_s:      Radiation window-start time
            rad_e:      Radiation window-end time
    '''
    # Ensure correct path handling
    param_get_file = os.path.join(dir_path, f"{out_base}1.hdf5")

    # Extract constants from output file
    with h5py.File(param_get_file, "r") as f:
        real_params = f["/base_param/real_params"]
        value1 = np.array(real_params["value1"])
        value2 = np.array(real_params["value2"])
        in_current = np.array(f["/input/input_env_vars"])

    # Extract relevant parameters
    i0 = value1[3]
    dx = value1[23]
    curr_max = value1[10]
    ebm_rms = value1[13]
    ebm_s = value1[14]
    ebm_e = value2[14]
    delt = value1[28]
    rad_s = value1[6]
    rad_e = value2[6]
    rad_wavlen = value1[1]

    runParams = [npass, curr_max, ebm_rms, i0, dx, rad_wavlen, \
                 ebm_s, ebm_e, delt, rad_s, rad_e]

    # Get time-axis associated with the current distribution
    curr_dist = in_current[1, :]
    t_ref = in_current[0, :]

    # Create time-axis for radiation field
    time_ax = np.arange(rad_e - delt/2, rad_s - delt/2, -delt)

    # Create radiation frequency-axis
    fcar = c / rad_wavlen
    f_ax = np.arange(len(time_ax)) / (time_ax[0] - time_ax[-1])
    frq_ax = f_ax - np.median(f_ax) + fcar

    return curr_dist, t_ref, frq_ax, time_ax, runParams



def extract_field_data(dir_path, out_base, npass):
    '''
    Extracts radiation data for a full run from Ginger3D output files.

    Ginger-3D file stores complex field data by concatenating the imaginary
    amplitudes to the real amplitudes in the dimension that we consider to have
    x-coordinates. We first extract the separate parts to get two equal-sized
    arrays of size (Ngrid,Ngrid,rad_time_slices) and use these parameters with
    npass to initialize two arrays of size (Ngrid,Ngrid,rad_time_slices, npass)
    which will store the full 3D complex radiation-field at each pass through
    the oscillator.

    Returns
    -------
    real_fld : Array (4D)
        Real-amplitudes of the 3D radiation field, stored for every pass. We
        use the first and second dimensions to index along the x and y axes
        respectively. The third dimension corresponds to the time-axis, while
        the last dimension identifies the oscillator pass number.
    imag_fld : Array (4D)
        Imaginary field amplitudes, complementing the data given by real_fld.
    '''
    # Ensure correct path handling
    first_fld_file = os.path.join(dir_path, f"fld{out_base}1.hdf5")

    # Get radiation field for first pass
    with h5py.File(first_fld_file, "r") as f:
        tmp_field = np.array(f["/radiation/2D_xy_field"])
    tmp_field = np.transpose(tmp_field)

    # Establish dimension sizes
    n_t_r_slices = tmp_field.shape[2]
    half_grid = tmp_field.shape[0] // 2

    # Initialize output arrays
    real_fld = np.zeros((half_grid, half_grid, n_t_r_slices, npass))
    imag_fld = np.zeros((half_grid, half_grid, n_t_r_slices, npass))

    # Separate the real and imaginary components of radiation field data
    fld_real = tmp_field[:half_grid, :, :]
    fld_imag = tmp_field[half_grid:, :, :]
    
    # Store first pass field data
    real_fld[:, :, :, 0] = fld_real
    imag_fld[:, :, :, 0] = fld_imag
    
    # Subsequent passes
    for pass_num in range(2, npass + 1):
        tmp_fld_file = os.path.join(dir_path, f"fld{out_base}2{pass_num}.hdf5")
        
        with h5py.File(tmp_fld_file, "r") as f:
            tmp_field = np.array(f["/radiation/2D_xy_field"]).T
        fld_real = tmp_field[:half_grid, :, :]
        fld_imag = tmp_field[half_grid:, :, :]
        
        real_fld[:, :, :, pass_num-1] = fld_real
        imag_fld[:, :, :, pass_num-1] = fld_imag
        
        if pass_num % 50:
            pass
        else:
            print(f"extracting field data at pass {pass_num} from {dir_path}")
        
    return real_fld, imag_fld



def extract_particle_data(dir_path, out_base, npass):
    '''
    Extracts electron-bunch data for a full run from Ginger3D output files.

    Ginger-3D file stores electron macroparticle data by concatenating the
    theta-data to the relativistic gamma-data in a new axis, separate from the
    number of macroparticles in the first dimension as well as the time slice
    in the last dimension. A single rst-file will store the data in an array of
    size (Nparts, 2, ebm_time_slices) so we separate the two datasets into two
    arrays of size (Nparts, ebm_time_slice) and initialize two arrays of size
    (Nparts, ebm_time_slice, npass) to store the data at each pass in the
    oscillator.

    Returns
    -------
    gam_dat : Array (3D)
        Relativistic gamma factor of the electron-bunch macroparticles saved at
        each pass. The first dimension specifies the particle, the second goes
        with bunch reference time-slice, while last describes the pass number.
    ang_dat : Array (3D)
        Longitudinal-phase data of the electron-bunch macroparticles at each
        oscillator pass.
    '''
    # Ensure correct path handling
    param_get_file = os.path.join(dir_path, f"{out_base}1.hdf5")

    # Get e-beam current distribution
    with h5py.File(param_get_file, "r") as f:
        in_current = np.array(f["/input/input_env_vars"])

    # Establish dimension sizes
    n_particles = 8192                          # modify this variable if wrong
    n_t_e_slices = in_current.shape[1]

    # Initialize output arrays
    gam_dat = np.zeros((n_particles, n_t_e_slices, npass))
    ang_dat = np.zeros((n_particles, n_t_e_slices, npass))    

    # Subsequent passes
    for pass_num in range(2, npass + 1):
        tmp_rst_file = os.path.join(dir_path, f"rst{out_base}2{pass_num}.hdf5")
        
        with h5py.File(tmp_rst_file, "r") as f:
            tmp_parts = np.array(f["/particles/gam-theta-data"]).T
        part_gamma = np.squeeze(tmp_parts[:, 0, :])
        part_theta = np.squeeze(tmp_parts[:, 1, :])
        
        gam_dat[:, :, pass_num-1] = part_gamma
        ang_dat[:, :, pass_num-1] = part_theta
        
        if pass_num % 50:
            pass
        else:
            print(f"extracting part-data from pass {pass_num} from {dir_path}")
        
    return gam_dat, ang_dat



