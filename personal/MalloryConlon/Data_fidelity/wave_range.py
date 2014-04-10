import os
import numpy as np
from math import *
import matplotlib.pyplot as plt

root = '/users/malloryconlon/astr596/data/cfa/'

list = []
all_range = []

for path, subdirs, files in os.walk(root):
    for name in files:
        if name.endswith('.flm'):
            list.append(os.path.join(path,name))

for i in range(len(list)):
    try:
        SN = np.genfromtxt(list[i])
        wavelength = SN[:,0]
        start = wavelength[0]
        end = wavelength[len(wavelength)-1]
        range = end-start
        all_range.append(range)
    #print range
    except:
            print list[i]

print len(all_range)
maximum = max(all_range)
print maximum

print all_range.index(max(all_range))


#print all_range

SN = np.genfromtxt(list[0])
wavelength = SN[:,0]
flux = SN[:,1]
start = wavelength[0]
end = wavelength[len(wavelength)-1]

print start
print end
print list[0]

plt.plot(wavelength,flux)
plt.show()