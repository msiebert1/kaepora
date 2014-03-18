import os
import numpy as np
from math import *

root = '/users/malloryconlon/astr596/data/cfa/'

list = []
all_range = []

for path, subdirs, files in os.walk(root):
    for name in files:
        list.append(os.path.join(path,name))

for i in range(len(list)):
    try:
        SN = np.genfromtxt(list[i])
        wavelength = SN[:,0]
        start = wavelength[0]
        end = wavelength[len(wavelength-1)]
        range = end-start
        all_range.append(range)
        print range
    except:
        print 'hey!'

print all_range

