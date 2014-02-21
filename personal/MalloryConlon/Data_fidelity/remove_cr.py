#This code uses a median filter to smooth out single-pixel deviations. Then, using sigma-clipping to remove large variations between the actual and smoothed image, we produce a cosmic ray cleaned image.

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt


#Read in data file and put wavelength, flux and error into separate arrays
SN=np.genfromtxt('sn1996ca-19970102.09-fast.flm')

wavelength = SN[:,0]
flux = SN[:,1]
error = SN[:,2]

#De-reshift the data
SN_info = np.genfromtxt('../cfasnIa_param.dat')
SN_name = SN_info[:,0]
z = SN_info[:,1]

for i in range(len(SN_name)):
    if SN_name[i]=='sn1996ca':
        wavelength /=1+z[i]

#Smooth out pixel-to-pixel variations and remove pixels that are deviant from the median by a factor of 3.  This removes cosmic rays fairly well and some sky lines, but it tends to miss galaxy lines.
flux_smooth = scipy.signal.medfilt(flux)
sigma = np.median(error)
bad = np.abs(flux - flux_smooth) / sigma > 3
flux_cr = flux.copy()
flux_cr[bad] = flux_smooth[bad]

plt.plot(wavelength,flux)
plt.plot(wavelength, flux_cr)
plt.show()