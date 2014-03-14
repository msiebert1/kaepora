#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from scipy import *
from datafidelity import *
import matplotlib.pyplot as plt
import pyfits
import glob
import os

cor = []


for dirs,subdirs,files in os.walk('../../data/cfa/'):
    for subdir in subdirs:
        list = glob.glob("*.flm")


#for i in range(len(list)):
#   SN = np.genfromtxt('../../../'+'list[i]')
#    wavelength = SN[:,0]
#   flux = SN[:,1]
#    var = SN[:,2]
#    var1 = genvar(wavelength,flux)
#    plt.plot(wavelength,var)
#    plt.plot(wavelength,var1)
#    plt.show()

#open = pyfits.open('bstarr.fits')
#sky = open[0].data

#cor = correlate(flux,sky)

#for i in range(len(cor)):
#   if cor[i]!=0:
#       print cor[i]



#plt.plot(sky_telluric)
#plt.show()





#plt.plot(wavelength,var)
#plt.plot(wavelength,var1)
#plt.show()