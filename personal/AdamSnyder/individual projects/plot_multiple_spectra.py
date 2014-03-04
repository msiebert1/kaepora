import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
import math
from scipy import interpolate

# Read spectra files from folder /cfa

spectra_files = glob.glob("../../data/cfa/*/*.flm")
spectra_arrays = []
bad_files = []

# Create list of spectra graphs

for i in range(20): # Set to read 20 SN spectra, can change
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
    except ValueError:
        bad_files.append(spectra_files[i]) # These files are not txt files, and aren't read for some reason

num_spectra = len(spectra_arrays)

# De-redshift data (Add code to get redshift data here)

# for i in range(num_spectra):
#    spectra_arrays[i][:, 0] /= 1 + redshifts[i]

# Interpolate data, first find min and max
# floor/ceiling to round to nearest integer wavelength

wave_min = math.floor(max(spectra[0, 0] for spectra in spectra_arrays))
wave_max = math.floor(min(spectra[-1, 0] for spectra in spectra_arrays))
wavelength = scipy.linspace(wave_min, wave_max, (wave_max-wave_min)*10)

new_spectra = [] # generate a new array of spectra graphs that share common x-values

for i in range(num_spectra):
    spline_rep = interpolate.splrep(spectra_arrays[i][:, 0], spectra_arrays[i][:, 1])
    new_flux = interpolate.splev(wavelength, spline_rep)
    new_flux /= np.median(new_flux)
    new_spectra.append([wavelength, new_flux])

# Generate composite spectra by averaging

Comp_spectra = np.mean(new_spectra, axis = 0)

res_flux = []

for i in range(num_spectra):
    res_flux.append(np.array(new_spectra[i][1]-Comp_spectra[1, :]))

res_flux = np.array(res_flux)
RMS = np.sqrt(np.mean(res_flux*res_flux, axis = 0))
flux_err_pos = Comp_spectra[1, :] + RMS
flux_err_neg = Comp_spectra[1, :] - RMS
scatter = RMS / Comp_spectra[1, :]

# Plot composite spectra (Jut need to figure out what to plot, and look up to make single figures)

plt.figure(1)
plt.subplot(211)
plt.plot(Comp_spectra[0, :], Comp_spectra[1, :], 'b' , Comp_spectra[0, :], flux_err_pos, 'r', Comp_spectra[0, :], flux_err_neg, 'g')
plt.title('Composite Spectra and (+/-) RMS')
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux')
plt.subplot(212)
plt.plot(wavelength, scatter)
plt.title('Residual')
plt.xlabel('Wavelength (A)')
plt.ylabel('Percantage')
plt.show()

                                         

