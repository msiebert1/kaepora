import numpy as np
import matplotlib.pyplot as plt
import datafidelity as df

# Import Supernova data to work with
SN=np.genfromtxt('sn2002cc-20020420-ui.flm')

wavelength = SN[:,0]
flux = SN[:,1]

try:
    variance = SN[:,2]
except IndexError:
    variance = df.genvar(wavelength, flux)

# Clip flux file
new_flux1, clipped = df.clip(flux)

# Smooth curve
new_flux2 = df.gsmooth(wavelength, flux, error, vexp = 0.004)

variance = df.update_variance(wavelength, new_flux1, variance)

plt.plot(wavelength, new_flux2, 'b')
plt.plot(wavelength, variance)
plt.show()
