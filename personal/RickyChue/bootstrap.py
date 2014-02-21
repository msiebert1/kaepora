# Authors: Ricky, Lunan
# Input number of tries (argv[1])

import numpy as np
import re, sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

rootdir = 'sn1995al/'

filename = open('name.txt', 'r').read()

list = filename.split('\n') # Name of the files.
num = len(list)             # Number of samples.



spec = [0] * num

for i in range(num):
    spec[i] = np.loadtxt(rootdir + list[i], unpack = True, usecols = (0, 1,))


num_arr = np.arange(0, num, 1)


tries = (int)(sys.argv[1])                 # Number of tries.

sel_spec = [0] * tries      # An array for the bootstraped spectra. (tries = 100)

for i in range(tries):      # For our try, there are 100 spectra.
    sel_spec[i] = np.abs(np.round(np.random.uniform(-0.49999999999, num - 0.5000000001, num))).astype(int)
    
    # Create a simple composite spectrum 
    wave = [0] * len(sel_spec[i])
    flux = [0] * len(sel_spec[i])

    for m in range(len(sel_spec[i])):
        wave[m] = spec[sel_spec[i][m]][0]
        flux[m] = spec[sel_spec[i][m]][1]
    

file = [0] * tries
spect = [0] * tries

for i in range(tries):
    file[i] = np.loadtxt('random/' + str(i) + '.dat', unpack = True)
    spect[i] = file[i][1]

### Standard deviation of the spectrum (for scatter plot)
stv = np.std(spect, axis = 0)

### Write the wavelength and scatter in a file.
outfile = open('scatter.dat', 'w')

for i in range(len(stv)):
    outfile.write(str(file[0][0][i]) + '\t' + str(stv[i]) + '\n')

plt.scatter(file[0][0], stv)
plt.ylim((0, 1))
plt.plot()
plt.show()

outfile.close()
    