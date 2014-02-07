import numpy as np
import glob
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from math import log
# This code needs to be made more efficients. Too many for-loops.

# Read data from repository /cfa/, need ALL the Supernovae data
spectra_files  = glob.glob('../data/cfa/*/*.flm')

spectra_arrays = []
bad_files      = []

for i in range(2):
    try:
        spectra_arrays.append(np.loadtxt(spectra_files[i]))
    except ValueError:
        bad_files.append(spectra_files[i])

num_spectra = len(spectra_arrays)

# Dereddining spectrum. It is so embarrassed, it is blushing  
red = np.genfromtxt('../data/cfa/cfasnIa_param.dat',dtype = None)

# Fun Facts: red[row][column]
# The for-loop needs to be able to check all the spectra_array file names
# and compare that the red[i][0] (the SN name). 
# It probably isn't smart enough to do that yet. Such a pity. 
for i in range(num_spectra):
	if spectra_arrays[i][:,0] == red[i][0]:
		spectra_arrays[i][:, 0] /= 1 + red[i][1]

# Interpolate Data. Yes.
wave_min    = max(spectra[0, 0] for spectra in spectra_arrays)
wave_max    = min(spectra[-1, 0] for spectra in spectra_arrays)
wavelength  = scipy.linspace(wave_min,wave_max,(wave_max-wave_min)*100)

new_spectra = []

for i in range(num_spectra):
	spline_rep = interpolate.splrep(spectra_arrays[i][:,0],spectra_arrays[i][:,1])
	new_flux = interpolate.splev(wavelength, spline_rep)
	new_flux /= np.median(new_flux)
	new_spectra.append([wavelength, new_flux])

# Generate a composite spectrum.
Comp_wave = wavelength
Comp_flux = []
sum_flux  = sum(spectra[1][:] for spectra in spectra_arrays)

for i in range(len(wavelength)):
	flux = sum(spectra[1][i] for spectra in new_spectra)
	Comp_flux.append(flux/float(num_spectra))

# Measure the scatter for the composite spectrum
res_flux_array = []
RMS            = []
flux_err_pos   = []
flux_err_neg   = []

for i in range(num_spectra): 
	res_flux_array = new_spectra[i][1] - Comp_flux[i]

RMS = np.array(RMS)

for i in range(num_spectra):
	RMS = np.sqrt(np.mean(np.multiply(res_flux_array[i],res_flux_array[i])))

flux_err_pos = np.add(Comp_flux,RMS)
flux_err_neg = Comp_flux - RMS

scatter = np.divide(RMS,Comp_flux)

# Plotting both the Comoposite and Residual in a two-panel plot
fig = plt.figure()

# Labeling the plots. Pt1.
ax1 = fig.add_subplot(211)
plot1, = ax1.plot(Comp_wave, Comp_flux,'b')
# Plot RMS Value once it is calculated and the correct dimensions
plot2, = ax1.plot(Comp_wave, Comp_flux,'m')
plt.title('Composite & RMS Spectrums')
plt.ylabel('Scaled Flux')
plt.legend([plot1,plot2],('Composite Flux','RMS'),'upper right',numpoints=1)

# Labeling. Pt2.
ax2 = fig.add_subplot(212)
ax2.plot(Comp_wave, res_flux_array,'k')
plt.title('Residual RMS Spectrum')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Scaled Flux')

plt.show()

