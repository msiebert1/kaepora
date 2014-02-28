#This code uses a median filter to smooth out single-pixel deviations. Then, using sigma-clipping to remove large variations between the actual and smoothed image, we produce a smoothed image.  To change the amount of smoothing, change the window_len parameter.  This is similar to a polynomial smoother i.e. The spectrum becomes more smooth as we increase the window_len. Once the original data is smoothed, we clip any data that is different from the smoothed data by 5%.

import pyfits
import numpy as np
import matplotlib.pyplot as plt


#Read in data file and put wavelength, flux and error into separate arrays.  Should make this to read in a list of spectra paths, then do the smoothing for that list.
SN=np.genfromtxt('sn2008z-20080229-ir.flm')

wavelength = SN[:,0]
flux = SN[:,1]
#error = SN[:,2]

def smooth(x,window_len=21,window='hanning'):
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
flux_update=flux

ratio = flux/new_flux

#Clip any bad data and replace it with the smoothed value

for i in range(len(ratio)):
    if ratio[i] > 1.05:
        flux[i] = new_flux[i]
        #print wavelength[i] #Uncomment to print clipped wavelengths
    if ratio[i] < 0.95:
        flux_update[i] = new_flux[i]
        #print wavelength[i] #Uncomment to print clipped wavelengths

#Plot old and new flux arrays vs wavelength to visually see changes

plt.plot(wavelength,flux,'k')
plt.plot(wavelength,flux_update,'r')
plt.show()

#Generate the variance based on the smoothed flux and original flux.  Subtract the two to get the noise value and smooth again to get the variance.

sky=pyfits.open('../../personal/malloryconlon/data_fidelity/kecksky.fits')


sky_flux=sky[0].data[0]

#print sky_flux

wavelength = SN[:,0]
flux = SN[:,1]
#error = SN[:,2]


noise=flux-new_flux
smooth_noise = smooth(noise)
#print noise

plt.plot(smooth_noise)
plt.plot(sky_flux)
plt.show()
