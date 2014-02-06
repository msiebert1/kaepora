import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate

# Read spectra files from folder

spectra_files = glob.glob("../../data/cfa/*/*.flm")
spectra_arrays = []
bad_files = []

# Create list of spectra graphs

for i in range(20):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
    except ValueError:
        bad_files.append(spectra_files[i]) # These files are not txt files, and aren't read for some reason

num_spectra = len(spectra_arrays)

# De-redshift data (Add code to get redshift data here)

# for i in range(num_spectra):
#    spectra_arrays[i][:, 0] /= 1 + redshifts[i]

# Interpolate data, first find min and max

wave_min = max(spectra[0, 0] for spectra in spectra_arrays)
wave_max = min(spectra[-1, 0] for spectra in spectra_arrays)
wavelength = scipy.linspace(wave_min, wave_max, (wave_max-wave_min)*100)

new_spectra = [] # generate a new array of spectra graphs that share common x-values

for i in range(num_spectra):
    spline_rep = interpolate.splrep(spectra_arrays[i][:, 0], spectra_arrays[i][:, 1])
    new_flux = interpolate.splev(wavelength, spline_rep)
    new_flux /= np.median(new_flux)
    new_spectra.append([wavelength, new_flux])

# Generate composite spectra by averaging

Comp_wave = wavelength
Comp_flux = []

for i in range(len(wavelength)): # For each wavelength, sum associated flux for each spectra and average 
    flux = sum(spectra[1][i] for spectra in new_spectra)
    Comp_flux.append(flux/float(num_spectra))

# Plot composite spectra

plt.plot(Comp_wave, Comp_flux)
plt.show()                                 
                                         

