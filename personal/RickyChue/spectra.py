### This program is to plot the spectra. Also, the average of the two spectra is obtained then the average is plotted.
### The link in line 9 has to be changed to fit your machine.

import numpy as np
import matplotlib.pyplot as plt
import os, glob, re

### Grab the file names. (You have to change the path to the directory where you put the spectra files.)
link = '/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596bak/data/'

os.chdir(link)

spectra = glob.glob('*.flm')

file0 = np.loadtxt(link + spectra[0], unpack = True)
file1 = np.loadtxt(link + spectra[1], unpack = True)


### Fluxes and wavelengths
lamda0 = file0[0]
lamda1 = file1[0]
flux0 = file0[1]
flux1 = file1[1]

min = np.minimum(np.min(flux0), np.min(flux1))
max = np.maximum(np.max(flux0), np.max(flux1))



### The graph for two spectra.
plt.scatter(lamda0, flux0, color = 'b', marker = '.')
plt.scatter(lamda1, flux1, color = 'r', marker = '.')
## Scaling the y axis
plt.ylim((min * 1.05, max * 1.05))
plt.xlabel(r'Wavelength, $\lambda$ (nm)')
plt.ylabel(r'Flux, $f_{\lambda}$ (Jy)')
plt.plot()
plt.savefig(link + 'two_spectra.png')
plt.close()



### The graph of the average for two spectra (as there are fewer data points for 'sn2011fe-visit3-hst.flm' than 'sn2011by-hst+lick.flm', the wavelength which 'sn2011fe-visit3-hst.flm' does not cover is deleted.
### The average of the two spectra is taken as the mean of the spectra.
lamda0_new = lamda0[np.searchsorted(lamda0, lamda1)]
flux0_new = flux0[np.searchsorted(lamda0, lamda1)]

### The mean of the two spectra
spectra_mean = np.mean([flux0_new, flux1], axis = 0)


### The graph for the mean spectra.
plt.scatter(lamda0_new, spectra_mean, color = 'b', marker = '.')
plt.ylim((np.min(spectra_mean) * 1.05, np.max(spectra_mean) * 1.05))
plt.xlabel(r'Wavelength, $\lambda$ (nm)')
plt.ylabel(r'Flux, $f_{\lambda}$ (Jy)')
plt.plot()
plt.savefig(link + 'average_spectra.png')
plt.close()

### Writing the average spectra in a file.
spectra_outfile = open(link + 'average_spectra.dat', 'w')

for i in range(flux0_new.size):
    spectra_outfile.write(str(lamda0_new[i]) + '\t' + str(spectra_mean[i]) + '\n')

spectra_outfile.close()









