import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from datafidelity import *
from math import *
import matplotlib.pyplot as plt

list = []
scale = []
scales = []

root = '/users/malloryconlon/astr596/data/cfa/'

for path, subdirs, files in os.walk(root):
    for name in files:
        list.append(os.path.join(path,name))

for i in range(len(list)):
    try:
        SN = np.genfromtxt(list[i])
        wavelength = SN[:,0]
        flux = SN[:,1]
        var = SN[:,2]
        var1 = genvar(wavelength,flux)
        #plt.plot(wavelength,var)
        #plt.plot(wavelength,var1)
        #plt.show()
        scale.append(var/var1)
        avg = np.average(scale)
        scales.append(avg)
    except:
        print 'Zero Variance'

avg_scales = np.average(scales)
print np.median(scales)
print np.average(scales)
plt.hist(scales,100)
plt.show()

#for i in range(len(list)):
#    name = list[i]
#    SN = np.genfromtxt(name)
#   wavelength = SN[:,0]
#   flux = SN[:,1]
#   var = SN[:,2]
#   var1 = genvar(wavelength,flux)
#   plt.plot(wavelength,var)
#   plt.plot(wavelength,var1)
#   plt.show()

