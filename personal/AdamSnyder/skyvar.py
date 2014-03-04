import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 
from math import *
import pyfits

sky = pyfits.open("kecksky.fits")
prihdr = sky[0].header
print prihdr

#plt.plot(wavelength, flux, 'b', wavelength, new_flux, 'r')
#plt.show()
