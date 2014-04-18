import datafidelity as df
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *

#SN = np.genfromtxt('../../data/spectra/bsnip/sn2004bz-20040613.453-ui.flm')
SN = np.genfromtxt('../data/spectra/cfa/sn1996ai/sn1996ai-19960620.17-fast.flm')

wavelength = SN[:, 0]
flux = SN[:, 1]
error = SN[:, 2]

ivar = df.genivar(wavelength, flux, error)

plt.plot(wavelength, flux)
plt.show()
