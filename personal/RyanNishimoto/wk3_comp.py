import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intp
import math
"""
Things to work on/Known issues:
-Clean by making more functions that are modular in design
-Make code more efficient by minimizing loops and improving code design 
(combine tasks under same loop, binary searches, etc)
-Indexes in loops don't line up when a file is bad and is added to 'junk'
	-Solution?  change the for loop 'in range' to reflect paths, not number of spectra to get
-Using linspace needs to be cleaner with the parameters
"""

#Functions
"""
this function reads in a file (from an array of paths)
and extracts the relevant data (Pathname, Wavelength and Flux),
appending the folder/file.flm path and Wavelength/flux into 
separate lists 'path' and 'data' respectively
"""
def get_data(f):
	try:
		data.append(np.loadtxt(f))	#adds spectra data into array 
		path.append(f[14:-4])		#adds SN path name into array
	except ValueError:
		junk.append(f)			#adds faulty files into array


"""
this function reads in, line by line, the file containing
SNe names and z values (among other stats) and appends
relevent information (right now name and z value) into
lists 'name' and 'z' respectively
"""
def get_stats(s):				#reads in particular row from stat file
	name.append(s[0])			#adds SN name into array
	z.append(s[1])				#adds z-value into array


"""
this function takes a single spectra 'path' and 'data', 
the lists 'name' and 'z', and de-redshifts the wavelength
of this particular data 
"""
def dez(path,data,name,z):
	for i in range(len(name)):				#Goes through list of names
		if name[i] in path:				#Current name is in path list
			print "matched SN",name[i],"to",path
			print "z =",z[i]
			print "...de-redshifting..."			
			for j in range(len(data)):		#Goes through each wavelength in single spectra
				data[j,0] /= 1+z[i]		#De-redshifts the wavelength 
			break	

	
#Reads in files and stats 
#Initializes variables

files = glob.glob ('../../data/cfa/*/*.flm')	#array of files
stats = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
name = []	#sn name from stats
z = []		#z from stats

data = []	#holds spectra data from files
junk = []	#holds junk files 

path = []	#holds pathname of files
wave = []	#holds wavelengths			
flux = []	#holds flux

num = 30		#measure num spectra
wave_min = 0	
wave_max = 10000

###########################
print "===================="
print "There are",len(files),"files to read"
print "There are",len(stats),"Supernovae with stats"
print "===================="
###########################

#loop for getting data and de-redshifting
for i in range(num):
	print"=== Data #",i,"==="

	#read in data
	get_data(files[i])	#gets pathname, wavelength, and flux from files
	get_stats(stats[i])	#gets name and z-value from stats file
	#de-redshift data

	print "working on path ",path[i]
	print "orig wavelengths", data[i][:,0]

	dez(path[i],data[i],name,z)

	print "new wavelengths",data[i][:,0]

	wave.append(data[i][:,0])
	flux.append(data[i][:,1])

print "\ncalculating min/max from:",wave_min,wave_max
print "===================="
for i in range(num):
	print "waves #",i
	if data[i][0][0] > wave_min:
		print data[i][0][0],">",wave_min,"\n-->","new min = ",data[i][0][0]
		wave_min = int(data[i][0][0])
	if data[i][-1][0] < wave_max:
		print data[i][-1][0],"<",wave_max,"\n-->","new max = ",data[i][-1][0]
		wave_max = int(data[i][-1][0])
		


print "\nusing min:",wave_min
print "using max:",wave_max
"""
Experienced Errors with truncation due to not exact values when finding indices;
fixed by interpolation

for i in range(num):
	print "\ntruncating data #",i
	print wave[i]
	low = np.where(np.logical_and(wave[i] > wave_min-1,wave[i] < wave_min+1))[0]
	print "minimum wave",wave_min,"at index",low
	high = np.where(np.logical_and(wave[i] > wave_max-1,wave[i] < wave_max+1))[0]
	print "maximum wave",wave_max,"at index",high
	wave[i] = wave[i][low:high]
	flux[i] = flux[i][low:high]
"""
print "\nCreating linspace from",wave_min,"to",wave_max,"with",(wave_max-wave_min),"points"
#creates space for wavelengths given min/max parameters
wavelengths = np.linspace(wave_min,wave_max,(wave_max-wave_min))


#Interpolating
fit = []
fit_flux = []
for i in range(num):
	spline = intp.splrep(wave[i],flux[i])
	curr_flux = intp.splev(wavelengths, spline)
	curr_flux /= np.median(curr_flux)
	fit.append([wavelengths,curr_flux])
	fit_flux.append(curr_flux) 


avg_flux = sum(fit_flux)/num

res_flux = []
#residual and RMS 
for i in range(num):
	res_flux.append(fit_flux[i]-avg_flux)

rms = np.sqrt(np.mean(np.power(res_flux,2)))
pos = avg_flux + rms
neg = avg_flux - rms
scatter = np.divide(rms,avg_flux)

############################
print "===================="
print "finished"
print "plotting..."
print "===================="
############################
fig = plt.figure()

#top plot containing the Composite spectrum +/- the RMS
top = fig.add_subplot(211)
plot1 = top.plot(wavelengths,avg_flux,'k')
plot2 = top.plot(wavelengths,pos,'b')
plot3 = top.plot(wavelengths,neg,'r')
plt.title('Composite Spectrum +/- RMS spectrum')
plt.ylabel('Relative Flux')
plt.legend([plot1,plot2,plot3],('Composite Flux','+RMS','-RMS'),'upper right',numpoints=1)

#bottom plot of the Residual
bottom = fig.add_subplot(212)
bottom.plot(wavelengths,scatter, 'k')
plt.title('Residual')
plt.xlabel('Wavelength')
plt.ylabel('Residual RMS')
plt.savefig('Composite_Spectra.png')
plt.show()

