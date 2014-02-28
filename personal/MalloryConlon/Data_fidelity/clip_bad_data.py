#This code uses a median filter to smooth out single-pixel deviations. Then, using sigma-clipping to remove large variations between the actual and smoothed image, we produce a cosmic ray cleaned image.

import numpy as np
import matplotlib.pyplot as plt


#Read in data file and put wavelength, flux and error into separate arrays.  Should make this to read in a list of spectra paths, then do the smoothing for that list.
SN=np.genfromtxt('sn2001l-20010201-ui.flm')

wavelength = SN[:,0]
flux = SN[:,1]
#error = SN[:,2]

def smooth(x,window_len=25,window='hanning'):
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



new_flux=smooth(flux)

ratio = flux/new_flux

for i in range(len(ratio)):
    if ratio[i] > 1.05:
        flux[i] = new_flux[i]
    if ratio[i] < 0.95:
        flux[i] = new_flux[i]

plt.plot(wavelength,flux)
#plt.plot(wavelength,new_flux)
plt.show()
