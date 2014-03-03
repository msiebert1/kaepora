import numpy as np
import matplotlib.pyplot as plt
import scipy 
from math import *

import curve_smoothing        

def genvar(filepath, vexp = 0.005, nsig = 5.0):

    SN = np.loadtxt(filepath)

    # unpack data from file and initialize objects
    wavelength = SN[:, 0]
    flux = SN[:, 1]
    try:
        varflux = SN[:, 2]
    except IndexError:
        varflux = np.zeros(len(wavelength))+1.0 # Space holder for sky spectrum (Add code here)

    new_flux = curve_smoothing.smoothing(wavelength, flux, varflux, vexp, nsig)

    # Generate absolute value of the noise from original flux
    error = abs(flux - new_flux)

    # Initialize array for smoothed error
    variance = curve_smoothing.smoothing(wavelength, error, varflux, vexp, nsig)    

    # Plots to check the variance flux
    plt.plot(wavelength, flux, 'b', wavelength, new_flux, 'r', wavelength, variance, 'g')
    plt.show()

    return variance

variance = genvar("../../data/bsnip/sn2002cc-20020420-ui.flm")

