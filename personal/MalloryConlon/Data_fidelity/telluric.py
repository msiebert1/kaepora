#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from clip_bad_data.py import smooth

SN = SN=np.genfromtxt('sn2002cc-20020420-ui.flm')
wavelength = SN[:,0]
flux = SN[:,1]

telluric_lines = np.loadtxt('../../personal/malloryconlon/Data_fidelity/telluric_lines.txt')

min = telluric_lines[:,0]
max = telluric_lines[:,1]

new_flux = smooth(flux)

ratio = flux/new_flux
flux_update = []

#Look at the flux/smoothed flux ratios for a given telluric absorption range as defined by the min and max arrays. If the ratio is less than the given condition, clip and replace with the smoothed flux value.


#Replace the inverse variance of each clipped value with 0