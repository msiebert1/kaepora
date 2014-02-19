import numpy as np
import matplotlib.pyplot as plt
import os, glob, re

# set the address of data

link = '/Users/fripsun/astr596/personal/lsun/'
address = '/Users/fripsun/astr596/data/'

# take the data out
os.chdir(address)

data = glob.glob('*.flm')

# read the data

text0 = np.loadtxt(address + data[0], unpack = True)

text1 = np.loadtxt(address + data[1], unpack = True)

print text0

# set coordinate, wavelength and flux
wavelength0 = text0[0]
wavelength1 = text1[0]
flux0 = text0[1]
flux1 = text1[1]

# set max and min for scaled plotting

miny = np.minimum(np.min(flux0), np.min(flux1))
maxy = np.maximum(np.max(flux0), np.max(flux1))

# plot two spectra

m = plt.scatter(wavelength0, flux0, color = 'c', marker = '.' )

n = plt.scatter(wavelength1, flux1, color = 'm', marker = '.' )

plt.legend((m,n), ('sn2011by', 'sn2011fe'),scatterpoints=1,loc='upper right', ncol=3)


#scaled plotting
plt.ylim((miny * 1.01, maxy * 1.01))
plt.xlabel(r'Wavelength, $\lambda$ ($10^{-9}$m)')
plt.ylabel(r'Flux, $f_{\lambda}$ (Jy)')

plt.plot()
plt.show()
plt.savefig(link + 'original spectra.png')
plt.close()


#since one set of data has more points, we want to truncate the excess points. For the rest, take average with the other data set and plot.

minx = np.minimum(wavelength0.size, wavelength1.size)

spectra_outfile = open(link + 'average_spectra.dat', 'w')

for i in range(minx):

#meanflux = (flux0[i]+flux1)/2  #np.mean([flux0,flux1], axis=0)
    plt.scatter(wavelength0[i], (flux0[i]+flux1[i])/2  , color = 'b', marker = '.')
    plt.ylim((miny * 1.01, maxy * 1.01))
    plt.xlabel(r'Wavelength, $\lambda$ ($10^{-9}$m)')
    plt.ylabel(r'Flux, $f_{\lambda}$ (Jy)')
    plt.plot()
    spectra_outfile.write(str(wavelength0[i]) + '\t' + str((flux0[i]+flux1[i])/2) + '\n')



plt.savefig(link + 'averaged spectra.png')

#plt.show()
plt.close()
