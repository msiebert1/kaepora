############################################################################
#
## This module contains the functions for data fidelity:
## Smoothing functions, clip_data functions, and variance generation
#
## To include these functions, use
#  from datafidelity import *
#
############################################################################

# Import necessary python modules
import numpy as np
from math import *

############################################################################
#
## Smoothing functions
#
## Function gsmooth() is an inverse variance weighted Gaussian smoothing of spectra
## Optional imputs are smoothing velocity (vexp) and number of sigma (nsig)
## Syntax: new_y_array = gsmooth(x_array, y_array, var_y, vexp = 0.01, nsig = 5.0)

def gsmooth(x_array, y_array, var_y, vexp = 0.005, nsig = 5.0):
    
    # Check for non-zero variance points, and set to 1E-20
    for i in range(len(var_y)):
        if var_y[i] == 0:
            var_y[i] = 1E-20
    
    # Output y-array
    new_y = np.zeros(len(x_array), float)
    
    # Loop over y-array elements
    for i in range(len(x_array)):
        
        # Construct a Gaussian of sigma = vexp*x_array[i]
        gaussian = np.zeros(len(x_array), float)
        sigma = vexp*x_array[i]
        
        # Restrict range to +/- nsig sigma
        sigrange = np.nonzero(abs(x_array-x_array[i]) <= nsig*sigma)
        gaussian[sigrange] = (1/(sigma*sqrt(2*pi)))*np.exp(-0.5*((x_array[sigrange]-x_array[i])/sigma)**2)
        
        # Multiply Gaussian by 1 / variance
        W_lambda = gaussian / var_y
        
        # Perform a weighted sum to give smoothed y value at x_array[i]
        W0 = np.sum(W_lambda)
        W1 = np.sum(W_lambda*y_array)
        new_y[i] = W1/W0
    
    # Return smoothed y-array
    return new_y

############################################################################
#
## Function wsmooth() smooths data by convolving input data with a window of specified size
## Optional inputs are window length (window_len) and window type (window)
## Syntax: new_flux_array = wsmooth(flux_array, window_len=17, window='hanning')

def wsmooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
        
        output:
        the smoothed signal
        
        example:
        
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also:
        
        np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
        scipy.signal.lfilter
        
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2):-(window_len/2)] instead of just y.
        """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    
    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window needs to be 'flat', 'hanning', 'hamming', 'bartlett', or 'blackman'"
    
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    
    return y[(window_len/2):-(window_len/2)]

############################################################################
#
# Function to generate variance for files that are missing this data
#
## Function gvar() generates a variance spectrum.
## Required inputs are an array of wavelengths (wavelength) and an array of corresponding fluxes (flux)
## Optional inputs are velocity of smoothing (vexp) [default 0.005] and number of sigma (nsig) [default 5.0]
## genvar(wavelength, flux, float vexp = 0.005, float nsig = 5.0)
#

def genvar(wavelength, flux, vexp = 0.005, nsig = 5.0):
    
    # Create variance from sky spectrum (Will add additional code here)
    varflux = np.zeros(len(wavelength))+1.0 # Place holder
    
    # Smooth original flux
    new_flux = gsmooth(wavelength, flux, varflux, vexp, nsig)
    
    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)
    
    # Smooth noise to find the variance
    variance = gsmooth(wavelength, error, varflux, vexp, nsig)
    
    # Return generated variance
    return variance

############################################################################
#
# Function clip() scans data and clips any bad points
# Required input is an array of fluxes
# Optional inputs are the upper and lower limits for the ratio
# Syntax is clip(flux_array, upper = 1.7, lower = 0.5)

def clip(flux, upper = 1.1, lower = 0.9):
    
    #Clip any bad data and replace it with the smoothed value.  Fine tune the ratio limits to cut more (ratio closer to one) or less (ratio farther from one) data
    
    new_flux = np.zeros(len(flux))
    clipped_points = []
    
    smooth_flux = wsmooth(flux)
    ratio = flux/smooth_flux
    
    for i in range(len(ratio)):
        if ratio[i] > upper:
            new_flux[i] = smooth_flux[i]
            clipped_points.append(i)
        elif ratio[i] < lower:
            new_flux[i] = smooth_flux[i]
            clipped_points.append(i)
        else:
            new_flux[i] = flux[i]
    
    return new_flux, clipped_points


############################################################################
#
# Function telluric_flag() scans data for telluric lines and flags any
# of the indices that have telluric contamination
# Required inputs are an array of wavelengths and an array of fluxes
# Optional input is the limit for flagging
# Syntax is telluric_flag(wavelength_array,flux_array, limit = 0.99)

def telluric_flag(wavelength, flux, limit=0.9):

    telluric_lines = np.loadtxt('telluric_lines.txt')

    min = telluric_lines[:,0]
    max = telluric_lines[:,1]

    new_flux = wsmooth(flux,window_len=35)

    ratio = flux/new_flux
    telluric_clip = np.zeros(len(flux))

#Look at the flux/smoothed flux ratios for a given telluric absorption range as defined by the min and max arrays. If the ratio is less than the given condition, clip and replace with the smoothed flux value.

    for i in range(len(wavelength)):
        for j in range(len(min)):
            if wavelength[i] > min[j]:
                if wavelength[i] < max[j]:
                    if ratio[i] < limit:
                        telluric_clip.append(i)

#Return the indices of telluric absorption
    return telluric_clip

############################################################################
#
# Function update_variance() uses the returns from the telluric_flag and clip
# functions to update the inverse variance spectrum
# Required inputs are an array of wavelengths, an array of fluxes, and the
# inverse variance array
# Syntax is update_variance(wavelength_array,flux_array, variance_array)

def update_variance(wavelength, flux, variance):

#Determine the clipped indices
    telluric_clip=telluric_flag(wavelength, flux)
    clipped_points=clip(flux)

    for i in range(len(clipped_points)):
            variance[clipped_points[i]] = 0


    for j in range(len(telluric_clip)):
        variance[telluric_clip[i]] = 0


    #Return the updated inverse variance spectrum
    return variance


