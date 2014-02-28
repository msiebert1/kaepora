# Authors: Ricky, Lunan
# Input number of tries (argv[1])

import numpy as np
import re, sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import erf

rootdir = '/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/MalloryConlon/Galaxy/'

filename = open(rootdir + 'MaxSpectra.dat', 'r').read()

list_long = filename.split('\r') # Name of the files.
new_num = len(list_long) - 1             # Number of samples.
num = new_num - 1

list = [0] * num    # Files to be input

for i in range(1, new_num):
    list[i - 1] = list_long[i].split()[1]

spec = [0] * num    # Input the wavelength and flux.


for i in range(num):
    spec[i] = np.loadtxt(list[i], unpack = True, usecols = (0, 1,))


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
    file[i] = np.loadtxt(list[i], unpack = True)
    spect[i] = file[i][1]


### To be completed by Sam's group.
composite = np.mean(spect, axis = 0)    ### Composite spectrum (just taking mean)
median = np.median(spect, axis = 0)     ### Median of the spectrum (for scaling)
scaled_comp = np.divide(composite, median)  ### Scaled spectrum
###



### 16th and 84th percentile of the spectrum (for scatter plot)
percentile = erf(1/np.sqrt(2.))

low_pc = 0.5 - percentile / 2.
up_pc = 0.5 + percentile / 2.


### The 16th and 84th percentile index
low_ind = np.round(tries * low_pc).astype(int)
up_ind = np.round(tries * up_pc).astype(int)

### Sort the fluxes in each wavelength, and put the 16th and 84th percentile fluxes into two arrays
low_arr = np.divide(np.sort(spect, axis = 0)[low_ind - 1], median)
up_arr = np.divide(np.sort(spect, axis = 0)[up_ind - 1], median)


### Write the wavelength and scatter in a file.
outfile = open('scatter.dat', 'w')

### Writing the wavelength, composite, 16th and 84th percentile to a file.
for i in range(len(low_arr)):
    outfile.write(str(file[0][0][i]) + '\t' + str(scaled_comp[i]) + '\t' + str(low_arr[i]) + '\t' + str(up_arr[i]) + '\n')

outfile.close()

"""
plt.plot(file[0][0], low_arr, color = 'r')
plt.plot(file[0][0], scaled_comp, color = 'b')
plt.plot(file[0][0], up_arr, color = 'k')
plt.legend(['16th', 'comp', '84th'])
plt.title('# of spectra = ' + str(tries))
#plt.ylim((0, 1.4))
plt.plot()
plt.show()
"""