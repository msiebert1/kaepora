import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 
from math import *

SN = np.loadtxt("../../data/bsnip/sn2002cc-20020420-ui.flm")

wavelength = SN[:, 0]
flux = SN[:, 1]
varflux = np.zeros(len(wavelength))+1.0
vexp = 0.01
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

error = flux - new_flux

error = abs(error)
new_error = np.zeros(len(error))

for i in range(len(wavelength)):
    gaussian = np.zeros(len(wavelength), float)
    sigma = vexp*wavelength[i]

    sigrange = np.nonzero(abs(wavelength-wavelength[i]) <= nsig*sigma)

    gaussian[sigrange] = (1/(sigma*sqrt(2*pi)))*np.exp(-0.5*((wavelength[sigrange]-wavelength[i])/sigma)**2)

    W_lambda = gaussian / varflux
    W0 = np.sum(W_lambda)
    W1 = np.sum(W_lambda*error)

    new_error[i] = W1/W0

plt.plot(wavelength, flux, 'b', wavelength, new_flux, 'r', wavelength, new_error, 'g')
plt.show()
