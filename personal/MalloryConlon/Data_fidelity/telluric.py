#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from datafidelity import wsmooth
import matplotlib.pyplot as plt

SN = SN=np.genfromtxt('SN04dt_040914_b01_DUP_WF.dat')
wavelength = SN[:,0]
flux = SN[:,1]

telluric_lines = np.loadtxt('../../personal/malloryconlon/Data_fidelity/telluric_lines.txt')

min = telluric_lines[:,0]
max = telluric_lines[:,1]

new_flux = wsmooth(flux,window_len=35)

ratio = flux/new_flux
telluric_clip = []

#Look at the flux/smoothed flux ratios for a given telluric absorption range as defined by the min and max arrays. If the ratio is less than the given condition, clip and replace with the smoothed flux value.

for i in range(len(wavelength)):
    for j in range(len(min)):
        if wavelength[i] > min[j]:
            if wavelength[i] < max[j]:
                if ratio[i] < 0.9:
                    telluric_clip.append(i)

print telluric_clip