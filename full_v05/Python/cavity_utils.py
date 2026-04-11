# -*- coding: utf-8 -*-
"""
@author: Levi K C Fisher

Comments:
cavity_utils is a Python module used by the UH Manoa Accelerator group made as
part of an Oscillator Add-On package to be used in FEL simulations.
The cavity_utils module defines functions that are useful for calculating some
of the beam quantities from the radiation-field, manipulating the state of the
radiation-field, and even transmitting the radiation-field through a cavity
with variable length.
"""

# Ext. Modules used
import numpy as np
import scipy.fft as fft

# Ext. Functions used
from scipy.interpolate import interp1d



# Constants
c = 299792458                                       # speed of light [meter/s]
pi = np.pi
i = complex(1j)                                     # imaginary unit



# Calculation (Find) Functions: used to analyze radiation field

def find_pwr(ar, ai, i0, dx):
    '''
    Function to calculate the power-frequency spectrum, power-time pulse,
    and the overall complex-phase at each time slice of the Radiation field.

    Parameters
    ----------
    ar : Array (3D)
        Real component of E-field amplitudes.
    ai : Array (3D)
        Imaginary component of E-field amplitudes.
    i0 : Float
        Intensity coefficient [W/m**2].
    dx : Float
        Grid spacing [m].

    Returns
    -------
    power_z : Array (1D)
        Power-frequency spectrum of radiation field
    power_s : Array (1D)
        Power-time pulse of radiation field
    angl : Array (1D)
        Overall complex-phase of radiation field
    '''
    # Define Complex E-field at each gridpoint (slice intersection)
    E_field = ar + i * ai
    
    # Calculate power at each time slice from E-field intensity
    I = np.abs(E_field)**2
    power_s = np.sum(I, axis=(0,1)) *i0 * dx**2
    
    # Calculate overall complex phase at each time slice
    Ea = np.sum(E_field, axis=(0,1))
    angl = np.unwrap( np.angle(Ea) )

    # Calculate frequency spectrum that contributes to the fields power
    Ef = np.abs(  fft.fftshift( fft.fft(E_field, axis=2), axes=2 )  )**2
    power_z = np.sum(Ef, axis=(0,1))
    
    return power_z, power_s, angl


def find_mode2(Eint_2d, xgrid, ygrid):
    '''
    Function to calculate the second moment mode of the Radiation field 
    intensity, i.e. the spot-size of the 2-D transverse intensity profile.

    Parameters
    ----------
    Eint_2d : Array (2D)
        DESCRIPTION.
    xgrid : Array (1D)
        DESCRIPTION.
    ygrid : Array (1D)
        DESCRIPTION.

    Returns
    -------
    sig_x2 : Float
        Squared second moment of intensity in x.
    sig_y2 : Float
        Squared second moment of intensity in y.
    '''
    # Determine grid-size
    dx = abs(xgrid[1] - xgrid[0])
    dy = abs(ygrid[1] - ygrid[0])


#def find_fwhm_pwr(pwr_s,rad_s,rad_e,delt):



# Optical Utility Functions: used to manipulate radiation field

def propagate_z(E3_fft, KKmesh, z_drift):
    '''
    Function to propagate 3D E-field in the z-direction using Fourier optics.

    Parameters
    ----------
    E3_fft : Array (3D)
        Initial 3D E-field in frequency (kx,ky) space.
    KK : Array (2D)
        Meshgrid of kx, ky frequency coordinates expressed as kz.
    z_drift : Float
        Distance in the z-direction to propagate E-field.

    Returns
    -------
    E3n_fft : Array (3D)
        The kx, ky frequency-space E-field after propagation.
    '''
    # Create propagation kernel
    H_z = np.exp( i*z_drift * KKmesh )
    
    # Apply kernel and return field
    E3n_fft = E3_fft * H_z[:,:,None]
    
    return E3n_fft


def mirror_reflect(E3_xyz, XYmesh, k_wavnum, focal_len):
    '''
    Function thar applies a mirror reflection to a 3D E-field beam.

    Parameters
    ----------
    E3_xyz : Array (3D)
        Input E-field in spatial (x,y,z) space.
    XYmesh : Array (2D)
        Summed squares of X,Y meshgrids as single array.
    k_wavnum : Float
        Free-space wavenumber of radiation.
    focal_len : Float
        Focal length of the mirror (lens).

    Returns
    -------
    E3n_xyz : Array (3D)
        The spatial domain of the E-field after being reflected.
    '''
    # Create thin-lens transfer function
    k_f = k_wavnum/(2*focal_len)
    H_lens = np.exp( -i*k_f * XYmesh )
    
    # Apply lens and return field
    E3n_xyz = E3_xyz * H_lens[:,:,None]
    
    return E3n_xyz



# Oscillator (Cavity) Functions: used to incorporate a beam into a cavity

def desynchronize(arn1,ain1,time_axs,my_time):
    '''
    Function to interpolate/shift the E-field to a "desynchronized" time-axis.
    Returns the radiation field as observed in a new reference frame.

    Parameters
    ----------
    arn1 : Array (3D, size = [i, j, k])
        Real component of E-field amplitudes.
    ain1 : Array (3D, size = [i, j, k])
        Imaginary components of E-field.
    time_axs : Array of size (k,)
        Vector array of 'input' times the E-field is defined on.
    my_time : Array of size (m,)
        Vector array of 'output' times to interpolate to.
    
    Returns
    -------
    ar : Array (3D, size = [i, j, m])
        The real components of the E-field on new time-axis.
    ai : Array (3D, size = [i, j, m])
        The imaginary E-field components on new time-axis.
    '''
    # Apply linear interpolation onto new time axis
    r_interp = interp1d(
        time_axs,
        arn1,
        axis=2,
        kind='linear',
        bounds_error=False,
        fill_value= 0 )

    i_interp = interp1d(
        time_axs,
        ain1,
        axis=2,
        kind='linear',
        bounds_error=False,
        fill_value= 0 )
    
    # Return interpolated amplitudes
    ar = r_interp(my_time)
    ai = i_interp(my_time)
    
    return ar, ai


def inser_lens(ar,ai, lamb, Z, LL, fN, dxN, desyn_s):
    '''
    Function for propagating a radiation beam into a cavity, incorporating 
    natural diffraction and focusing lenses using Fourier optics. The length of
    the cavity can be adjusted via "desynchronization" with the convention that
    positive desynchronization corresponds to cavity shortening while negative
    corresponds to cavity lengthening. The desynchronization value d is defined
    as the ratio between the change in cavity length over the slippage length.

    Parameters
    ----------
    ar : Array (3D)
        Real component of E-field that exits undulator.
    ai : Array (3D)
        Imaginary component of E-field that exits undulator.
    lamb : Float
        Radiation wavelength in meters.
    Z : Float
        Original length of the full cavity in meters.
    LL : Float
        Length of the undulator in meters.
    fN : Float
        Focal length of the mirror lens in meters.
    dxN : Float
        2D grid spacing.
    desyn_s : Float
        The change in length of the cavity, given in units of meters.

    Returns
    -------
    arn : Array (3D)
        Real component of field amplitude as it would enter the undulator.
    ain : Array (3D)
        Imaginary component of E-field entering undulator.
    Waist_list : List of Floats     # size=(4,)
        Waist-size of radiation beam at different points in cavity, as follows:
            wx_in:  Exiting undulator
            wx_z1:  Before first focusing lens
            wx_z2:  Before second focusing lens
            wx_z1b: Before entering undulator
    '''
    ##### Setup #####
    # Define field parameters
    k0 = 2*pi / lamb                                # Free-space wavenumber
    
    # Calculate drift distances
    z1  = (Z-LL-desyn_s)/2                          # Undulator to far lens
    z2  = Z - desyn_s/2                             # Far lens to close lens
    z1b = (Z-LL)/2                                  # Close lens to undulator
    
    # Establish grid (and padding) parameters
    padd = 100                                      # Number/dim of grids added
    N1   = ar.shape[0]
    mid  = len(ar[0,0,:]) // 2                      # Central time-slice index
    dx = dy = dxN                                   # Grid step size in x,y
    
    p2 = padd // 2
    N  = N1 + padd                                  # Total number of x,y grids
    Lx = dx*N                                       # Physical size in x
    Ly = dy*N                                       # Physical size in y
    
    # Set up spatial grid arrays
    x    = ( np.arange(N) - (N-1)/2 ) * dx          # X,Y-spatial coordinates
    y    = ( np.arange(N) - (N-1)/2 ) * dy
    X, Y = np.meshgrid(x, y)                        # Meshgrid spatial coords
    
    # Set up frequency grid arrays
    kx = ( np.arange(N) - (N-1)/2 )* (2*pi / Lx)    # X-frequency coordinates
    ky = ( np.arange(N) - (N-1)/2 )* (2*pi / Ly)    # Y-frequency coordinates
    
    KX, KY = np.meshgrid(kx, ky)                    # Meshgrid frq coordinates
    
    # Combine meshgrids
    XY = X**2 + Y**2
    KK = np.sqrt(k0**2 - KX**2 - KY**2)
    
    # Pad the x,y axes of input field
    E_in0= (ar + i * ai)
    E_in = np.pad(E_in0, ( (p2,p2), (p2,p2), (0,0) ), \
                  mode='constant', constant_values=0)
    
        
    ##### Propagation #####
    # Fourier transform E-field to propagate in z direction (to back lens)
    E_in_fft = fft.fftshift( fft.fft2(E_in, axes=(0,1)), axes=(0,1) )
    E_z1_fft = propagate_z(E_in_fft, KK, z1)
    
    # Transform to spatial domain to apply lens (located at z=Z-desyn_s)
    E_z1 = fft.ifft2( fft.ifftshift(E_z1_fft, axes=(0,1)), axes=(0,1) )
    E_l1 = mirror_reflect(E_z1, XY, k0, fN)
    
    # Back to Fourier domain for propagation (backwards to the front lens)
    E_l1_fft = fft.fftshift( fft.fft2(E_l1, axes=(0,1)), axes=(0,1) )
    E_z2_fft = propagate_z(E_l1_fft, KK, z2)
    
    # Apply the front mirror lens in spatial domain (located at z=0)
    E_z2 = fft.ifft2( fft.ifftshift(E_z2_fft, axes=(0,1)), axes=(0,1) )
    E_l2 = mirror_reflect(E_z2, XY, k0, fN)
    
    # Transform to Fourier domain and finally propagate to undulator entrance
    E_l2_fft = fft.fftshift( fft.fft2(E_l2, axes=(0,1)), axes=(0,1) )
    E_z1b_fft= propagate_z(E_l2_fft, KK, z1b)
    
    
    ##### Complete #####
    # Return the spatial-domain field at undulator entrance (at z=(Z-LL)/2)
    E_z1b= fft.ifft2( fft.ifftshift(E_z1b_fft,axes=(0,1)), axes=(0,1) )
    arn = E_z1b.real
    arn = arn[p2:N-p2, p2:N-p2, : ]
    ain = E_z1b.imag
    ain = ain[p2:N-p2, p2:N-p2, : ]
    
    """
    # Also return the beam waist at different stages
    myE0 = E_in[:,:, mid]                           # Undulator exit
    wx_in, wy_in = find_waistxy(myE0, x, y)
    
    myE1 = E_z1[:,:, mid]                           # Before far lens
    wx_z1, wy_z1 = find_waistxy(myE1, x, y)
    
    myE2 = E_z2[:,:, mid]                           # Before close lens
    wx_z2, wy_z2 = find_waistxy(myE2, x, y)
    
    myE3 = E_z1b[:,:,mid]                           # Undulator entrance
    wx_z1b, wy_z1b = find_waistxy(myE3, x, y)
    
    
    # Group the waists into a list
    Waist_list = [wx_in, wx_z1, wx_z2, wx_z1b]
    
    return arn, ain, Waist_list"""
    return arn, ain
