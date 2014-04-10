import numpy as np
import pyfits
import math
import datafidelity as df
import matplotlib.pyplot as plt

SN = np.genfromtxt('sn1989m-19890710-oi1i2.flm')

wavelength = np.array(SN[:,0])
flux = np.array(SN[:,1])

ivar = df.genivar(wavelength, flux)
