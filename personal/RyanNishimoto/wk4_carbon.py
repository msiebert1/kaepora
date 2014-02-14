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
as well as keeping a counter of the number of spectra for a 
particular SN
"""
def get_data(f,p,d,c,i):
	try:
		d.append(np.loadtxt(f))		#adds spectra data into array 
		p.append(f[14:-4])		#adds SN path name into array
		c[i] += 1			#number of spectra for particular SN
	except ValueError:
		junk.append(f)			#adds faulty files into array


"""
this function reads in, line by line, the file containing
SNe names and z values (among other stats) and appends
relevent information (right now name and z value) into
lists 'name' and 'z' respectively
"""
def get_stats(s,n,z,t):				#reads in particular row from stat file
	n.append(s[0])				#adds SN name into array
	z.append(s[1])				#adds z-value into array
	t.append(s[2])				#adds time of max light

"""
this function takes a single spectra 'path' and 'data', 
the lists 'name' and 'z', and de-redshifts the wavelength
of this particular data 
"""
def dez(d,p,n,z):
	for i in range(len(n)):				#Goes through list of names
		if n[i] in p:				#Current name is in path list
			#print"redshift wavelengths:",(data[:,0])
			#print "matched SN",name[i],"to",path
			#print "z =",z[i]
			#print "...de-redshifting..."			
			for j in range(len(data)):		#Goes through each wavelength in single spectra
				d[j,0] /= 1+z[i]		#De-redshifts the wavelength 
			#print"de-redshifted values:",(data[:,0]),"\n"
			break	

	
#Reads in files and stats 
#Initializes variables

carbon_pos = np.loadtxt('files/wk4_carbon_neg.txt',dtype = str)
carbon_neg = np.loadtxt('files/wk4_carbon_neg.txt',dtype = str)

files = glob.glob ('../../data/cfa/*/*.flm')	#array of files
stats = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
ages = np.genfromtxt('../../data/cfa/cfasnIa_mjdspec.dat',dtype = None)	#age of each spectra file


name_pos = []	#sn name from stats
name_neg = []

z_pos = []	#z from stats
z_neg = []

tmax_pos = []	#time of maxiumum light from stats
tmax_neg = []

data_pos = []	#holds spectra data from files
data_neg = []

count_pos = [0]*len(carbon_pos)		#array to hold number of spectra per SN
count_neg = [0]*len(carbon_neg)

junk = []	#holds junk files 

path_pos = []	#holds pathname of files
path_neg = []

num_pos = len(carbon_neg)		#measure num spectra
num_neg = len(carbon_neg)

wave_min = 3800
wave_max = 7400


###########################
print "===================="
print "There are",len(files),"files to read"
print "There are",len(stats),"Supernovae with stats"
print "There are",len(carbon_pos),"Carbon positive (CP) Supernovae"
print "There are",len(carbon_neg),"Carbon negative (CN) Supernovae"
print "===================="
print "reading in data..."
print "checking for corrupt files..."
###########################


#for each file, find the matching SN, record the number of spectra for that SN, and extract
#Wavelength/Flux for each one
for file in files:
	for i in range(len(carbon_pos)):
		if carbon_pos[i] in file:
			get_data(file,path_pos,data_pos,count_pos,i)

print "\nsearch complete"		
print "===================="
print "found",len(data_pos),"CP files to read"
print "found",len(data_neg),"CN files to read"
print "found",len(junk),"corrupt files"
print "===================="
print "file breakdown"
print "===================="
for i in range(len(carbon_pos)):
	print "CP sn",carbon_pos[i],"has",count_pos[i],"files"

count_pos = [0]*len(carbon_pos)		#array to hold number of spectra per SN
count_neg = [0]*len(carbon_neg)

#for i in range(len(carbon_neg)):
	#print "CN sn",carbon_neg[i],"has",count_neg[i],"files"
print "\ntrimming the fat...(removing SNe with no files)"




carbon_new = []
count_new = []
for i in range(len(count_pos)):
	if (count_pos[i] != 0):
		carbon_new.append(carbon[i])
		count_new.append(count[i])

print "searching for relevant SNe statistics...\n"
#for each row of stats, find the corresponding redshift and add it to an array
for stat in stats:
	for i in range(len(carbon_new)):
		if carbon_new[i] in stat[0]:
			get_stats(stat,name_pos,z_pos,tmax_pos)



"""
~Sanity Check~  The files we currently have:

carbon_new	length 18	contains names of SN we want
count_new	length 18	contains the integer number of files per SN
z		length 18	contains redshift data

data		length 211	contains wavelength and flux
path		length 211	contains path to data

"""
#de-redshifting data
for i in range(len(data_pos)):
	
	dez(data_pos[i],path_pos[i],carbon_new,z_pos)
	


num = len(data_pos)

#this code will determine the data to take based off of the the given ranges
print "\nfinding spectra in range:",wave_min,wave_max

delete = []
for i in range(num):
	#print "\nwaves #",i,":",data_pos[i][:,0]
	#print "from path:",path_pos[i]
	if data_pos[i][0][0] > wave_min or data_pos[i][-1][0] < wave_max:
		delete.append(i)
		#print "\nwaves #",i,":",data_pos[i][:,0]
		#print "from path:",path_pos[i],"does not make the cut"

for i in range(len(delete)):
	print "remove",data_pos[delete[len(delete)-1-i]]
	del data_pos[delete[len(delete)-1-i]]
	print data_pos[delete[len(delete)-1-i]]

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


for i in range(len(data_pos)):
	spline = intp.splrep(data_pos[i][:,0],data_pos[i][:,1])
	curr_flux = intp.splev(wavelengths, spline)
	curr_flux /= np.median(curr_flux)
	fit.append([wavelengths,curr_flux])
	fit_flux.append(curr_flux) 


avg_flux = sum(fit_flux)/num


res_flux = []
#residual and RMS 
for i in range(len(data_pos)):
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
plt.savefig('plots/Carbon_pos.png')
plt.show()

"""#labeling plot1 carbon +
top = fig.add_subplot(4,1,1)
plot1 = top.plot(wavelengths,avg_flux,'k')
plot2 = top.plot(wavelengths,pos,'b')
plot3 = top.plot(wavelengths,neg,'r')
plt.title('Composite Spectrum: Carbon +')
plt.ylabel('Relative Flux')
plt.legend([plot1,plot2,plot3],('Carbon +','+ RMS','- RMS'),'upper right',numpoints=1)

#labeling plot2 carbon -
mid = fig.add_subplot(4,1,2)
#for carbon -, need a different avg_flux, pos, neg
plot4 = mid.plot(wavelengths,avg_flux,'k')
plot5 = mid.plot(wavelengths,pos,'b')
plot6 = mid.plot(wavelengths,neg,'r')
plt.title('Composite Spectrum: Carbon -')
plt.ylabel('Relative Flux')
plt.legend([plot4,plot5,plot6],('Carbon - ','+ RMS','- RMS'),'upper right',numpoints=1)

#labeling plot3 residual
bottom = fig.add_subplot(4,1,4)
bottom.plot(wavelengths,scatter, 'k')
plt.title('Residual')
plt.xlabel('Wavelength')
plt.ylabel('Residual RMS')

plt.savefig('Composite_Spectra.png')
plt.show()
"""
