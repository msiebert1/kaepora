import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy

SN = np.loadtxt("sn1994ae-19941129.51-fast.flm")

wavelength = SN[:, 0]
flux = SN[:, 1]
error = SN[:, 2]

print len(wavelength)

for i in range(len(wavelength)):
    if i < 3:
        l = 0
        r = i+4
    elif i > len(wavelength)-4:
        l = i-3
        r = len(wavelength)-1
    else:
        l = i-3
        r = i+4
    median = np.median(flux[l:r])
    if (flux[i]-median)/float(flux[i]) > 0.10:
        wave_inter = np.concatenate((wavelength[l:i], wavelength[i+1:r]))
        flux_inter = np.concatenate((flux[l:i], flux[i+1:r]))
        spline_rep = interpolate.splrep(wave_inter, flux_inter)
        new_flux_point = interpolate.splev(wavelength[i], spline_rep)
        flux[i] = new_flux_point
        error[i] = float("inf")

plt.plot(wavelength, flux, 'r')
plt.ylabel("Flux")
plt.xlabel("Wavelength")
plt.show()
