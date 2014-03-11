#This code uses a median filter to smooth out single-pixel deviations. Then, using sigma-clipping to remove large variations between the actual and smoothed image, we produce a smoothed image.  To change the amount of smoothing, change the window_len parameter.  This is similar to a polynomial smoother i.e., the spectrum becomes more smooth as we increase the window_len. Once the original data is smoothed, we clip any data that is different from the smoothed data by given factors that can be changed depending on desired smoothness.

import pyfits
import numpy as np
import matplotlib.pyplot as plt


#Read in data file and put wavelength, flux and error into separate arrays.  Should make this to read in a list of spectra paths, then do the smoothing for that list.
SN=np.genfromtxt('../../../data/bsnip/sn1998bp-19980921-ui.flm')

wavelength = SN[:,0]
flux = SN[:,1]

def smooth(x,window_len=55,window='hamming'):
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



new_flux = smooth(flux)
flux_update = []
clipped = []

ratio = flux/new_flux

#Clip any bad data and replace it with the smoothed value.  Fine tune the ratio limits to cut more (ratios closer to one) or less (ratios farther from one) data

for i in range(len(ratio)):
    if ratio[i] > 1.95:
        flux_update.append(new_flux[i])
        clipped.append(i)
    elif ratio[i] < 0.05:
        flux_update.append(new_flux[i])
        clipped.append(i)
    else:
        flux_update.append(flux[i])


#print clipped[i] #Uncomment to print clipped indices.  This can be used to fix the variance values after they are generated

#Plot old and new flux arrays vs wavelength to visually see changes

plt.plot(wavelength,flux,'k')
plt.plot(wavelength,flux_update,'r')
plt.show()




