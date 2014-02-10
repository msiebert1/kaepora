import numpy as np
import os
import scipy
from astropy.table import Table
import scipy.interpolate as inter
import glob
import matplotlib.pyplot as plt
import math

#Reads in spectra file names
spectra_files = glob.glob("../../../data/cfa/*/*.flm")

spectra_arrays=[]
bad_files = []
file_name =[]

num=20 #len(spectra_files)

for i in range(num):
    try:
		spectra_arrays.append(Table.read(spectra_files[i],format='ascii'))
		file_name.append(spectra_files[i])
    except ValueError:
        bad_files.append(spectra_files[i])


#deredshift data
parameters = Table.read('../../../data/cfa/cfasnIa_param.dat',format='ascii')
sn_name = parameters["col1"]
sn_z = parameters["col2"]
for i in range(len(file_name)):
	old_spectrum=spectra_arrays[i]
	z=0
	#file_name = spectra_files[i]
	for j in range(len(sn_name)):
		if sn_name[j] in file_name[i]:
			z=sn_z[j]
	lambda_obs=old_spectrum["col1"]
	lambda_emit= lambda_obs/(1+z)
	spectra_arrays[i]=Table([lambda_emit,old_spectrum["col2"]],names=('col1','col2'))

	
#scale spectra		
wave_min=0  #arbitrary minimum of wavelength range
wave_max=1000000   #arbitrary maximum of wavelength range

for i in range(len(spectra_arrays)):
	spectra = spectra_arrays[i]
	if (min(spectra["col1"]) > wave_min): #changes minimum wavelength if larger than previous
		wave_min=min(spectra["col1"])
	if (max(spectra["col1"]) < wave_max):  #changes maximum wavelength if smaller than previous
		wave_max=max(spectra["col1"])


wavelength = np.linspace(wave_min,wave_max,wave_max-wave_min)  #creates 100 equally spaced wavelength values between the smallest range

fitted_flux=[]
#generates composite spectrum
for i in range(len(spectra_arrays)):
	new_spectrum=spectra_arrays[i]	#declares new spectrum from list
	new_wave=new_spectrum["col1"]	#wavelengths
	new_flux=new_spectrum["col2"]	#fluxes
	lines=np.where((new_wave>wave_min) & (new_wave<wave_max))	#creates an array of wavelength values between minimum and maximum wavelengths from new spectrum
	sm1=inter.splrep(new_wave[lines],new_flux[lines])	#creates b-spline from new spectrum
	y1=inter.splev(wavelength,sm1)	#fits b-spline over wavelength range
	y1 /= np.median(y1)
	fitted_flux.append(y1)

avg_flux = np.mean(fitted_flux,axis=0)

"""
sum_flux = fitted_flux[0]
#sums fluxes to be averaged
for i in range(len(fitted_flux[0])):
	sum_flux[i] = 0
	
for i in range(len(fitted_flux)):
	sum_flux = sum_flux + fitted_flux[i]

avg_flux=sum_flux/len(spectra_arrays)	#averages fluxes
"""

avg_spectrum=Table([wavelength,avg_flux],names=('col1','col2'))	#puts together the average spectrum

rms_flux_max=[]
rms_flux_min=[]
delta=[]
scatter=[]
#RMS Spectrum, Residual

for i in range(len(fitted_flux)):
	delta.append(avg_flux-fitted_flux[i])

rms_flux = np.sqrt(np.mean(np.square(delta),axis=0))

print len(rms_flux)
print len(avg_flux)

for i in range(len(rms_flux)):
	scatter.append(rms_flux[i]/avg_flux[i]*100)
	
"""

for i in range(len(wavelength)):
	flux_sum = 0
	for j in range(len(fitted_flux)):
		flux_sum = flux_sum + (avg_flux[i]-fitted_flux[j])**2
	#flux_rms=(flux_sum)**.5
	flux_rms=(flux_sum/len(fitted_flux))**.5
	rms_max=avg_flux[i] + flux_rms
	rms_min=avg_flux[i] - flux_rms
	rms_flux_max.append(rms_max)
	rms_flux_min.append(rms_min)
	scatter.append(flux_rms/avg_flux[i]*100)


for i in range(len(wavelength)):
	flux_sum = 0
	for j in range(len(spectra_arrays)):
		spectrum=spectra_arrays[j]
		flux = spectrum["col2"]
		flux_sum = flux_sum + (flux[i]-avg_flux[i])**2
	#flux_rms=(flux_sum)**.5
	flux_rms=(flux_sum/len(spectra_arrays))**.5
	rms_max=avg_flux[i] + flux_rms
	rms_min=avg_flux[i] - flux_rms
	rms_flux_max.append(rms_max)
	rms_flux_min.append(rms_min)
	scatter.append(flux_rms/avg_flux[i]*100)	


for i in range(len(wavelength)):
	chi_square = 0
	a_flux = avg_flux[i]
	for j in range(len(spectra_arrays)):
		spectrum=spectra_arrays[j]
		flux = spectrum["col2"]
		chi_square = chi_square + ((flux[i]-a_flux)**2/flux[i])
	flux_rms= chi_square
	rms_max=a_flux + flux_rms
	rms_min=a_flux - flux_rms
	rms_flux_max.append(rms_max)
	rms_flux_min.append(rms_min)
	scatter.append(flux_rms/a_flux)
	"""
	
#rms_max_spectrum=Table([wavelength,rms_flux_max],names=('col1','col2'))
#rms_min_spectrum=Table([wavelength,rms_flux_min],names=('col1','col2'))

#RMS Residual
plt.figure(1)
plt.subplot(211)
plot1,=plt.plot(wavelength,avg_flux+rms_flux,label= 'rms+')
plot2,=plt.plot(wavelength,avg_flux-rms_flux, label= 'rms-')
plot3,=plt.plot(wavelength,avg_flux,label='comp')
legend=plt.legend(loc='upper right', shadow=True)
plt.xlim(wave_min,wave_max)
plt.ylabel('Flux')
plt.subplot(212)
plot1,=plt.plot(wavelength,scatter,label='rms residual')
plt.xlim(wave_min,wave_max)
plt.xlabel('Wavelength')
plt.ylabel('RMS Flux/ Average Flux')
plt.savefig('rmsplot.png')
plt.show()

