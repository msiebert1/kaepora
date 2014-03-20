#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from scipy import *
from datafidelity import *
import matplotlib.pyplot as plt
import pyfits
import glob
import os


root = '/users/malloryconlon/astr596/data/cfa/'

list = []
scales = []

for path, subdirs, files in os.walk(root):
    for name in files:
        list.append(os.path.join(path,name))
print list

for i in range(len(list)):
    try:
        SN = np.genfromtxt(list[i])
        wavelength = SN[:,0]
        flux = SN[:,1]
        var = SN[:,2]
        var1 = genvar(wavelength,flux)
        new = var/var1
        scale = np.average(new)
        scales.append(scale)
    except:
        print list[i]

print np.average(scales)


plt.hist(scales)

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