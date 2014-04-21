import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table
import msgpack as msg
import msgpack_numpy as mn
from scipy.optimize import leastsq
import datafidelity as df

SN = np.genfromtxt('../../data/spectra/cfa/sn1994ae/sn1994ae-19941129.51-fast.flm')

wavelength = SN[:,0]
flux = SN[:,1]
error = SN[:,2]
ivar = df.genivar(wavelength, flux, error)

plt.plot(wavelength, ivar)
plt.show()
