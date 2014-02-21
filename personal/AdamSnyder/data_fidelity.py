import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 

SN = np.loadtxt("../../data/bsnip/sn2002cc-20020420-ui.flm")

wavelength = SN[:, 0]
flux = SN[:, 1]
try:
    error = SN[:, 2]
    err_col = bool(1)
except IndexError:
    err_col = bool(0)
    
for i in range(len(wavelength)):
    if i < 4:
        l = 0
        r = i+5
    elif i > len(wavelength)-5:
        l = i-4
        r = len(wavelength)-1
    else:
        l = i-4
        r = i+5
    median = np.median(flux[l:r])
    if (flux[i]-median)/float(flux[i]) > 0.5:
        print wavelength[i]
        wave_inter = np.concatenate((wavelength[l:i], wavelength[i+1:r]))
        flux_inter = np.concatenate((flux[l:i], flux[i+1:r]))
        spline_rep = interpolate.splrep(wave_inter, flux_inter)
        new_flux_point = interpolate.splev(wavelength[i], spline_rep)
        flux[i] = new_flux_point
        if err_col == True:
            error[i] = float("inf")


plt.plot(wavelength, flux, 'r')
plt.ylabel("Flux")
plt.xlabel("Wavelength")
plt.show()
