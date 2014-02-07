import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intp
import math
"""
Things to work on:
-Clean by making more functions that are modular in design
-Make code more efficient by minimizing loops and improving code design 
(combine tasks under same loop, binary searches, etc)
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
		data.append(np.loadtxt(f))	#adds spectras info into data array
		path.append(f[14:-4])		#adds SNe name into path array
	except ValueError:
		junk.append(f)	#accounts for some faults in files


"""
this function reads in, line by line, the file containing
SNe names and z values (among other stats) and appends
relevent information (right now name and z value) into
lists 'name' and 'z' respectively
"""
def get_stats(s):			#reads in particular row from stat file
	name.append(s[0])			#gets name of sn from stats
	z.append(s[1])				#gets z-val of sn from stats


"""
this function takes a single spectra 'path' and 'data', 
the lists 'name' and 'z', and de-redshifts the wavelength
of this particular data 
"""
def dez(path,data,name,z):
	for i in range(len(name)):
		if name[i] in path:
			### Found SN corresponding to file
			### and adjusts with proper z
			print "matched",name[i],"to",path
			print "Redshift is",z[i]
			print "...calculating..."
			for j in range(len(data)):
				data[j,0] /= 1+z[i]
				data[j,0] = int(data[j,0])
			break	

	
#Variable initialization
#and
#Read in files to an array and gathers stats
files = glob.glob ('../../data/cfa/*/*.flm')	#array of files
stats = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
name = []	#sn name from stats
z = []		#z from stats

data = []	#holds spectra data from files
junk = []	#holds junk files 

path = []	#holds pathname of files
wave = []	#holds wavelengths			
flux = []	#holds flux

num = 2		#measure num spectra
wave_min = 0	
wave_max = 10000
###########################
### Pre-loop Debug area ###
###########################
print "===================="
print "There are",len(files),"files to read"
print "There are",len(stats),"Supernovae with stats"
print "initiating loop"
print "===================="
###########################

"""
Three Loops for C++ coders on their desktop computers
Seven for the C# coders on their mobile devices
Nine for the Java coders doomed to mis-manage memory
One for the Python Lord on his dark throne
In the Land of Type 1a Supernovae where Dark Energies lie
One Loop to rule them all, One Loop to find them,
One Loop to bring them all and in the darkness composite them
In the Land of Type 1a Supernovae where Dark Energies lie
"""
#loop for getting data and de-redshifting
for i in range(num):
	###########################
	### Per-loop Debug area ###	
	###########################
	#print "loop",i
	#print "===================="
	###########################

	#read in data
	get_data(files[i])	#gets pathname, wavelength, and flux from files
	get_stats(stats[i])	#gets name and z-value from stats file
	#de-redshift data

	#print "working on path ",path[i]
	#print "orig wavelengths", data[i][:,0]
	#print "...processing..."

	dez(path[i],data[i],name,z)
	#print "new wavelengths",data[i][:,0]
	wave.append(data[i][:,0])
	flux.append(data[i][:,1])

print "calculating range from start:",wave_min,wave_max
for i in range(num):
	if data[i][0][0] > wave_min:
		wave_min = data[i][0][0]
		print "new min",wave_min
	if data[i][-1][0] < wave_max:
		wave_max = data[i][-1][0]
		print "new max",wave_max

print "truncating data..."
for i in range(num):
	low = np.where(np.logical_and(wave[i] > wave_min-1,wave[i] < wave_min+1))[0]
	print wave[i]
	print low
	high = np.where(np.logical_and(wave[i] > wave_max-1,wave[i] < wave_max+1))[0]
	print high
	wave[i] = wave[i][low:high]
	flux[i] = flux[i][low:high]


#creates space for wavelengths given min/max parameters
wavelengths = np.linspace(wave_min,wave_max,(wave_max-wave_min)/5)

#Interpolating loop
fit= []
fit_flux = []
for i in range(num):
	spline = intp.splrep(wave[i],flux[i])
	curr_flux = intp.splev(wavelengths, spline)
	fit.append([wavelengths,curr_flux])
	fit_flux.append(curr_flux)


#averaging flux
sum_flux = []

for i in range(len(fit_flux[0])):
	sum_flux.append(0)
for i in range(num):
	sum_flux += fit_flux[i]
avg_flux = sum_flux/num

res_flux = []
#residual and RMS 
for i in range(num):
	res_flux.append(fit_flux[i]-avg_flux)

rms = np.sqrt(np.mean(np.power(res_flux,2)))
pos = fit_flux + rms
neg = fit_flux - rms

############################
### Post-loop Debug area ###
############################
print "===================="
print "finished"
############################
print "plotting..."

fig, (x0, x1) = plt.subplot(nrows=2, sharex=True)
plot1 = x0.plot(wavelength,fit_flux,label = 'comp')
plot2 = x0.plot(wavelength,pos,label = 'rms+')
plot3 = x0.plot(wavelength,neg,label = 'rms-')
legend = x0.legend(loc='lower right')
x0.set_xlim(wave_min,wave_max)
x0.set_title('RMS spectrum')
x0.set_ylabel('Flux')
x1.plot(wavelength,scatter,'o-')
x1.set_ylabel('Residual')
x1.set_ylim(0,1)
legend = ax1.legend(loc='lower right')

plt.subplots_adjust(hspace=0.1)
plt.show()
plt.savefig('all.png')
