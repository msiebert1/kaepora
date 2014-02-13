import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intp
import math



#authors @nkazmi2, @nishimo1 
#notes for Carbon +/-
"""
Jeffery Silverman Paper:
SNe Ia with C --> "bluer colours and lower luminosities at maximum light"
"We concentrate on the region 5600-7800 A in order to cover the spectral range
around Si II 6355 and C II 6580 as well as C II 7234 and O I triplet[near 7770]"
-velocities tend to be around 12000-13000 km/s for C II 6580
-velocities tend to be around 9500 -11500 km/s for C II 7234

Folatelli et al.(2012) says a relationship between color and light-curve width was 
shown to be steeper for objects with C


"""
"""
Things to work on/Known issues(wk3 code):
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
		wave.append(data[i][:,0])
		flux.append(data[i][:,1])
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
def dez(data,path,name,z):
	for i in range(len(name)):				#Goes through list of names
		if name[i] in path:				#Current name is in path list
			#print"redshift wavelengths:",(data[:,0])
			#print "matched SN",name[i],"to",path
			#print "z =",z[i]
			#print "...de-redshifting..."			
			for j in range(len(data)):		#Goes through each wavelength in single spectra
				data[j,0] /= 1+z[i]		#De-redshifts the wavelength 
			#print"de-redshifted values:",(data[:,0]),"\n"
			break	

	
#Reads in files and stats 
#Initializes variables

carbon = np.loadtxt('wk4_carbon.txt',dtype = str)

files = glob.glob ('../../data/cfa/*/*.flm')	#array of files
stats = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
name = []	#sn name from stats
z = []		#z from stats

data = []	#holds spectra data from files
junk = []	#holds junk files 

path = []	#holds pathname of files
wave = []	#holds wavelengths			
flux = []	#holds flux

num = len(carbon)		#measure num spectra
wave_min = 0	
wave_max = 10000


###########################
print "===================="
print "There are",len(files),"files to read"
print "There are",len(stats),"Supernovae with stats"
print "There are",len(carbon),"Carbon positive Supernovae we want"
print "===================="
print "reading in data..."
print "checking for corrupt files..."
###########################
count = [0]*len(carbon)

#for each file, find the matching SN, record the number of spectra for that SN, and extract
#Wavelength/Flux for each one
for file in files:
	for i in range(len(carbon)):
		if carbon[i] in file:
			try:
				#print "adding",carbon[i],"data:\n",file[14:-4]
				data.append(np.loadtxt(file))	#adds spectra data into array 
				path.append(file[14:-4])		#adds SN path name into array
				count[i] += 1
			except ValueError:
				#print "poop file",file
				junk.append(file)			#adds faulty files into array

print "\nsearch complete"		
print "===================="
print "found",len(data),"files to read"
print "found",len(junk),"corrupt files"
print "===================="
print "file breakdown"
print "===================="
for i in range(len(carbon)):
	print carbon[i],"has",count[i],"files"
print "\ntrimming the fat...(removing SNe with no files)"

carbon_new = []
count_new = []
for i in range(len(count)):
	if (count[i] != 0):
		carbon_new.append(carbon[i])
		count_new.append(count[i])

print "searching for relevant SNe statistics...\n"
#for each row of stats, find the corresponding redshift and add it to an array
for stat in stats:
	for i in range(len(carbon_new)):
		if carbon_new[i] in stat[0]:
			print carbon_new[i],"=",stat[0],"with z =",stat[1]
			z.append(stat[1])



"""
~Sanity Check~  The files we currently have:

carbon_new	length 18	contains names of SN we want
count_new	length 18	contains the integer number of files per SN
z		length 18	contains redshift data

data		length 211	contains wavelength and flux
path		length 211	contains path to data

"""
#de-redshifting data
for i in range(len(data)):
	
	dez(data[i],path[i],carbon_new,z)
	


num = len(data)
print "\ncalculating min/max from:",wave_min,wave_max
print "===================="
for i in range(num):
	#print "waves #",i
	if data[i][0][0] > wave_min:
		#print data[i][0][0],">",wave_min,"\n-->","new min = ",data[i][0][0]
		wave_min = int(data[i][0][0])
	if data[i][-1][0] < wave_max:
		#print data[i][-1][0],"<",wave_max,"\n-->","new max = ",data[i][-1][0]
		wave_max = int(data[i][-1][0])



print "\nhighest minimum:",wave_min
print "lowest maximum:",wave_max

wave_min -= (wave_min % 10)
wave_max -= (wave_max % 10)
print "\nusing pretty min:",wave_min 
print "using pretty max:",wave_max 


print "\nCreating linspace from",wave_min,"to",wave_max,"with",(wave_max-wave_min),"points"
#creates space for wavelengths given min/max parameters
wavelengths = np.linspace(wave_min,wave_max,(wave_max-wave_min+1))

print "plots will have data at these wavelength values:",wavelengths

#Interpolating
fit = []
fit_flux = []


for i in range(num):
	spline = intp.splrep(data[i][:,0],data[i][:,1])
	curr_flux = intp.splev(wavelengths, spline)
	curr_flux /= np.median(curr_flux)
	fit.append([wavelengths,curr_flux])
	fit_flux.append(curr_flux) 


avg_flux = sum(fit_flux)/num


res_flux = []
#residual and RMS 
for i in range(num):
	res_flux.append(fit_flux[i]-avg_flux)

rms = np.sqrt(np.mean(np.square(res_flux),axis = 0))
pos = avg_flux + rms
neg = avg_flux - rms
scatter = np.divide(rms,avg_flux)

############################
print "===================="
print "calculations finished"
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
plt.savefig('Carbon.png')
plt.show()

