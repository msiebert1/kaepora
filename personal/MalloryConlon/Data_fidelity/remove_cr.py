import numpy as np
import scipy.signal
import matplotlib.pyplot as plt


SN=np.genfromtxt('sn1994ae-19941129.51-fast.flm')

wavelength = SN[:,0]
flux = SN[:,1]
error = SN[:,2]

flux_smooth = scipy.signal.medfilt(flux, 5)
sigma = np.median(error)
bad = np.abs(flux - flux_smooth) / sigma > 2.0
flux_cr = flux.copy()
flux_cr[bad] = flux_smooth[bad]

plt.plot(flux_cr)
plt.show()