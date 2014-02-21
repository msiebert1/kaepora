# Authors: Ricky, Lunan
# Input number of tries (argv[1])

import numpy as np
import re, sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import erf

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
    sel_spec[i] = np.floor(np.random.uniform(0, num, num)).astype(int)
    
    
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

composite = np.mean(spect, axis = 0)    ### Composite spectrum (just taking mean)


### 16th and 84th percentile of the spectrum (for scatter plot)
percentile = erf(1/np.sqrt(2.))

low_pc = 0.5 - percentile / 2.
up_pc = 0.5 + percentile / 2.


### The 16th and 84th percentile index
low_ind = np.round(tries * low_pc).astype(int)
up_ind = np.round(tries * up_pc).astype(int)

### Sort the fluxes in each wavelength, and put the 16th and 84th percentile fluxes into two arrays
low_arr = np.sort(spect, axis = 0)[low_ind - 1]
up_arr = np.sort(spect, axis = 0)[up_ind - 1]


### Write the wavelength and scatter in a file.
outfile = open('scatter.dat', 'w')

### Writing the wavelength, composite, 16th and 84th percentile to a file.
for i in range(len(low_arr)):
    outfile.write(str(file[0][0][i]) + '\t' + str(composite[i]) + '\t' + str(low_arr[i]) + '\t' + str(up_arr[i]) + '\n')

plt.plot(file[0][0], low_arr, color = 'r')
plt.plot(file[0][0], composite, color = 'b')
plt.plot(file[0][0], up_arr, color = 'k')
plt.legend(['16th', 'comp', '84th'])
plt.title('# of spectra = ' + str(tries))
plt.ylim((0, 1.4))
plt.plot()
plt.show()

outfile.close()
