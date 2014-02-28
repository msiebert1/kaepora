import numpy as np
import pyfits
import matplotlib.pyplot as plt


SN=np.genfromtxt('sn1999ac-19990309.53-fast.flm')

sky=pyfits.open('../../../personal/malloryconlon/data_fidelity/kecksky.fits')


sky_flux=sky[0].data[0]
print sky_flux

wavelength = SN[:,0]
flux = SN[:,1]
error = SN[:,2]

print error

#noise=flux-sky_flux

#plt.plot(sky_flux)
plt.plot(error)
plt.show()