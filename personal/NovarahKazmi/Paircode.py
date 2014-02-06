import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from math import log
"""
! Create program which reads in at least 20 spectra (from multiple SNe)
>>> deredshifts them
! scales them
! and generates a composite spectrum.
>>> Measure the scatter for the composite spectrum. 
! Make a two-panel plot which shows the composite spectrum with the RMS spectrum around it on top and the residual RMS spectrum on bottom.
"""

spectra_files  = glob.glob("../astr596/data/cfa/*/*.flm")
spectra_arrays = []
bad_files      = []

for i in range(2):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
    except ValueError:
        bad_files.append(spectra_files[i])

num_spectra = len(spectra_arrays)

#Dereddining spectrum

#Want to skip the first line in cfasnIa_mjdspec.dat
"""
with open("../astr596/data/cfa/cfasnIa_mjdspec.dat") as redshift_files:
    next(redshift_files)
   
SN_name  = [] 
redshift = []
redshift_files[:,0] = SN_name
redshift_files[:,1] = redshift

# need to match SN_name with the spectra file being redshifted
for i in range(num_spectra):
	if spectra_arrays[i][:,0] == SN_name[i]:
		spectra_arrays[i][:, 0] /= 1 + redshifts[i]
	else continue:

print spectra_arrays
"""
# Interpolate Data
wave_min    = max(spectra[0, 0] for spectra in spectra_arrays)
wave_max    = min(spectra[-1, 0] for spectra in spectra_arrays)
wavelength  = scipy.linspace(wave_min,wave_max,(wave_max-wave_min)*100)

new_spectra = []

for i in range(num_spectra):
	spline_rep = interpolate.splrep(spectra_arrays[i][:,0],spectra_arrays[i][:,1])
	new_flux = interpolate.splev(wavelength, spline_rep)
	new_flux /= np.median(new_flux)
	new_spectra.append([wavelength, new_flux])

# Generate a composite spectrum
Comp_wave = wavelength
Comp_flux = []
sum_flux  = sum(spectra[1][:] for spectra in spectra_arrays)

for i in range(len(wavelength)):
	flux = sum(spectra[1][i] for spectra in new_spectra)
	Comp_flux.append(flux/float(num_spectra))

# Measure the scatter for the composite spectrum
"""
# Need to fill in information, general steps are written out
res_flux     = []
rms          = []
flux_err_pos = []
flux_err_neg = []

res_flux = new_spectra[] - Comp_flux[]

rms = sqrt(mean(res_flux[]^2))

flux_err_pos = Comp_flux + rms
flux_err_neg = Comp_flux - rms

scatter = rms[]/Comp_flux[]
"""


# Plotting both the Comoposite and Residual in a two-panel plot
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(Comp_wave, Comp_flux)
plt.title('Composite & RMS Spectrums')
#plt.xlabel('Wavelength (??)')
plt.ylabel('Scaled Flux')
ax2 = fig.add_subplot(212)
ax2.plot(Comp_wave, Comp_flux)
plt.title('Residual RMS Spectrum')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Scaled Flux')
plt.show()

