import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 
from math import *

SN = np.loadtxt("../../data/cfa/sn1993ac/sn1993ac-19931016.49-mmt.flm")

wavelength = SN[:, 0]
flux = SN[:, 1]
varflux = SN[:, 2]
for i in range(len(varflux)):
    if varflux[i] == 0:
        varflux[i] = 1E-20
vexp = 0.005
new_flux = np.zeros(len(wavelength), float)

outflux = np.array([])

nsig = 5.0

for i in range(len(wavelength)):
    gaussian = np.zeros(len(wavelength), float)
    sigma = vexp*wavelength[i]

    sigrange = np.nonzero(abs(wavelength-wavelength[i]) <= nsig*sigma)

    gaussian[sigrange] = (1/(sigma*sqrt(2*pi)))*np.exp(-0.5*((wavelength[sigrange]-wavelength[i])/sigma)**2)

    W_lambda = gaussian / varflux
    W0 = np.sum(W_lambda)
    W1 = np.sum(W_lambda*flux)

    new_flux[i] = W1/W0

plt.plot(wavelength, flux, 'b', wavelength, new_flux, 'r')
plt.show()
