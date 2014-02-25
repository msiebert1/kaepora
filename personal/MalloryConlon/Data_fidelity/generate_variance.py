import numpy as np
import pyfits
import matplotlib.pyplot as plt


SN=np.genfromtxt('sn1989a-19890427-o1i1.flm')

sky=pyfits.open('../../personal/malloryconlon/data_fidelity/kecksky.fits')
sky.info()

print sky[0].header[0]

wavelength = SN[:,0]
flux = SN[:,1]



#noise=flux-sky_flux

#plt.plot(sky_flux)
#plt.show()