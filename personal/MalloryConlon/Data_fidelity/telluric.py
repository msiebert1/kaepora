#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from clip_bad_data import smooth
import matplotlib.pyplot as plt

SN = SN=np.genfromtxt('sn1998bp-19980921-ui.flm')
wavelength = SN[:,0]
flux = SN[:,1]

telluric_lines = np.loadtxt('../../personal/malloryconlon/Data_fidelity/telluric_lines.txt')

min = telluric_lines[:,0]
max = telluric_lines[:,1]

new_flux = smooth(flux,window_len=35)

ratio = flux/new_flux
flux_update = []
clipped = []

#Look at the flux/smoothed flux ratios for a given telluric absorption range as defined by the min and max arrays. If the ratio is less than the given condition, clip and replace with the smoothed flux value.

for i in range(len(wavelength)):
    flux_update.append(flux[i])
    for j in range(len(min)):
        if wavelength[i] > min[j]:
            if wavelength[i] < max[j]:
                if ratio[i] < 0.99:
                    flux_update[i]=(new_flux[i])
                    clipped.append(i)

print len(flux_update)
print len(wavelength)
print clipped

plt.plot(wavelength,flux,'y')
plt.plot(wavelength,flux_update,'k')
plt.show()

#Replace the inverse variance of each clipped value with 0