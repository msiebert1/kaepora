import datafidelity as df
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *

#SN = np.genfromtxt('../../data/spectra/bsnip/sn2004bz-20040613.453-ui.flm')
SN = np.genfromtxt('../personal/AdamSnyder/sn1996C-19960217.48-fast.flm')

wavelength = SN[:, 0]
flux = SN[:, 1]

ivar = df.genivar(wavelength, flux)

plt.plot(wavelength, ivar)
plt.show()
