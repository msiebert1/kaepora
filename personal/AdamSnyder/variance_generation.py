import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 

SN = np.loadtxt("../../data/bsnip/sn2002cc-20020420-ui.flm")

wavelength = SN[:, 0]
flux = SN[:, 1]
varflux = []
vexp = 2

outflux = np.array([])

nsig = 5

for i in range(len(wavelength)):
    gaussian = np.zeros(len(wavelength), float)
    sigma = vexp*wavelength[i]

    sigrange = [for point in wavelength if point >= wavelength[i]-nsig*sigma and point <= wavelength[i]+nig*sigma)

    gaussian = np.random.normal(wavelength[i], sigma)
