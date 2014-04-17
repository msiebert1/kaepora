import datafidelity as df
import numpy as np
import matplotlib.pyplot as plt

def clip(wave, flux, ivar):
    # Create an array of all ones
    var = np.ones(len(flux), float)
    
    # Create 2 smoothed fluxes, of varying vexp
    sflux = df.gsmooth(wave, flux, var, 0.002)
    sflux2 = df.gsmooth(wave, flux, var, 0.004)
    
    # Take the difference of the two fluxes and smooth
    err = abs(flux - sflux)
  
    serr = df.gsmooth(wave, err, var, 0.008)

    plt.plot(wave, sflux, 'r')

    # Find the wavelengths that need to be clipped (omitting 5800-6000 region)
    bad_wave = wave[np.where((err/serr > 4) & ((wave < 5800.0) | (wave > 6000.0)))]

    # Find indices for general clipping
    bad = np.array([], int)
    for i in range(len(bad_wave)):
        bad = np.append(bad, np.where(abs(wave - bad_wave[i]) < 8))

    # Set ivar to 0 for those points and return
    ivar[bad] = 0
    return ivar

#SN = np.genfromtxt('../../data/spectra/bsnip/sn2004bz-20040613.453-ui.flm')
SN = np.genfromtxt('../personal/AdamSnyder/sn1996C-19960217.48-fast.flm')

wavelength = SN[:, 0]
flux = SN[:, 1]

plt.plot(wavelength, flux)

ivar = df.genivar(wavelength, flux)

ivar_new = clip(wavelength, flux, ivar)

plt.plot(wavelength, flux)
plt.show()
