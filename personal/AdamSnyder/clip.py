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
    serr = df.gsmooth(wave, err, var, 0.0008)

    # Find the wavelengths that need to be clipped (omitting 5800-6000 region)
    bad_wave = wave[np.where((err/serr > 3.5) & ((wave < 5800.0) | (wave > 6000.0)))]

    # Find indices for wavelengths and surrounding wavelengths to clip
    bad2 = np.array([], int)
    for i in range(len(bad_wave)):
        bad2 = np.append(bad2, np.where(abs(wave-bad_wave[i]) < 8))

    # Set ivar to 0 for those points and return
    ivar[bad2] = 0
    return ivar

SN = np.genfromtxt('sn1995bd-19960125.23-fast.flm')

wavelength = SN[:, 0]
flux = SN[:, 1]
error = SN[:, 2]

ivar = df.genivar(wavelength, flux, error)

ivar_new = clip(wavelength, flux, ivar)

plt.plot(wavelength, ivar_new)
plt.show()
