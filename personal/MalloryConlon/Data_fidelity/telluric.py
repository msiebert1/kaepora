#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

import numpy as np
from datafidelity import wsmooth
import matplotlib.pyplot as plt

SN = SN=np.genfromtxt('../../../SN04dt_040914_b01_DUP_WF.dat')
wavelength = SN[:,0]
flux = SN[:,1]

telluric_lines = np.loadtxt('../../personal/malloryconlon/Data_fidelity/telluric_lines.txt')

mi = telluric_lines[:,0]
ma = telluric_lines[:,1]

new_flux = wsmooth(flux,window_len=35)
flux1 = []

ratio = flux/new_flux
telluric_clip = []


#Look at the flux/smoothed flux ratios for a given telluric absorption range as defined by the min and max arrays. If the ratio is less than the given condition, clip and replace with the smoothed flux value.

for i in range(len(wavelength)):
    for j in range(len(mi)):
        if wavelength[i] > mi[j]:
            if wavelength[i] < ma[j]:
                if ratio[i] < 0.99:
                    telluric_clip.append(i)

for k in range(len(flux)):
    flux1.append(flux[k])
for l in range(len(telluric_clip)):
    index=telluric_clip[l]
    flux1[l]=new_flux[l]



plt.plot(wavelength,flux,'k')
plt.plot(wavelength,flux1,'y')
plt.show()