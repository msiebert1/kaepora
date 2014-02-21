import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intp
import math

#'#'#'#'#
import Plot_2
#'#'#'#'#


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

TASKS
[x] read/plot both C +/-
[] age range
[] weight error
[] median thingy

[]fix problem with ~7000 wavelengths
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
			for j in range(len(d)):		#Goes through each wavelength in single spectra
				d[j,0] /= 1+z[i]		#De-redshifts the wavelength 
			#print"de-redshifted values:",(data[:,0]),"\n"
			break	

"""
this function cleans the data by removing the files that
are outside of the range of wavelengths we want to look at
"""
def chop_data(d,p):
	delete = []	#temporary variable to track which indices to remove

	for i in range(len(d)):
		if d[i][0][0] > wave_min or d[i][-1][0] < wave_max:
			delete.append(i)
			#print "\nwaves #",i,":",data_pos[i][:,0]
			#print "from path:",path_pos[i],"does not make the cut"
	for i in range(len(delete)):
		#print "remove",d[delete[len(delete)-1-i]]
		del d[delete[len(delete)-1-i]]	#have to remove from end of array to account for moving indexes
		del p[delete[len(delete)-1-i]]
		#print d[delete[len(delete)-1-i]]
"""
this function cleans the lists of SN we are looking at
by removing the ones that have no corresponding file
"""
def chop_sn(sn,count):
	delete = []
	
	for i in range(len(count)):
		#print i
		if (count[i] == 0):
			#print sn[i],"has",count[i],"files"
			delete.append(i)
	for j in range(len(delete)):
		del sn[delete[len(delete)-1-j]]
		del count[delete[len(delete)-1-j]]
"""
This function interpolates the data. It reads in the wavelength and flux
file for both C + and C -. Using built in functions we complete the 
interpolation process and scale the data using the median. 

In later code we will not use the median to scale the data. 
We will scale the the data to agree with a "good spectrum". 
"""

def interpolate(d,ff):
	for i in range(len(d)):
		spline = intp.splrep(d[i][:,0],d[i][:,1])
		curr_flux = intp.splev(wavelengths, spline)
		curr_flux /= np.median(curr_flux)
		#fit.append([wavelengths,curr_flux])
		ff.append(curr_flux) 


#Reads in files, stats(z,time of max, etc), and ages(time of each spectra)
carbon_pos = np.loadtxt('files/wk4_carbon_pos.txt',dtype = str).tolist()
carbon_neg = np.loadtxt('files/wk4_carbon_neg.txt',dtype = str).tolist()


files = glob.glob ('../../data/cfa/*/*.flm')	#array of files
stats = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe
ages = np.genfromtxt('../../data/cfa/cfasnIa_mjdspec.dat',dtype = None)	#age of each spectra file

#Initializes variables
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

num_pos = len(carbon_pos)		#measure num spectra
num_neg = len(carbon_neg)

fit_flux_pos = []		#new fitted/scaled flux (following interpolation)
fit_flux_neg = []

wave_min = 3800
wave_max = 6500


###########################
print "\n===================="
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
			continue
	for j in range(len(carbon_neg)):
		if carbon_neg[j] in file:
			get_data(file,path_neg,data_neg,count_neg,j)
			continue
"""
print "\nsearch complete"		
print "===================="
print "found",len(data_pos),"CP files to read"
print "found",len(data_neg),"CN files to read"
print "found",len(junk),"corrupt files"
print "===================="
print "file breakdown"
print "===================="
"""
"""
for i in range(len(carbon_pos)):
	print "CP sn",carbon_pos[i],"has",count_pos[i],"files"
print "total of",sum(count_pos),"CP files found"
for i in range(len(carbon_neg)):
	print "CN sn",carbon_neg[i],"has",count_neg[i],"files"
print "total of",sum(count_neg),"CN files found"
"""

#remove files outside of range of wavelengths wanted
#print "\nfinding spectra in range:",wave_min,wave_max



#print "\noriginal data/path of length",len(data_pos),"/",len(path_pos)
chop_data(data_pos,path_pos)
#print "copped data/path of length",len(data_pos),"/",len(path_pos)

print "\nChopping C + complete."

#print "\noriginal data/path of length",len(data_neg),"/",len(path_neg)
chop_data(data_neg,path_neg)
print "\nChopping C - complete."


#removes sn names we don't want to look at
#print "\ntrimming the fat...(removing SNe with no files)"
#print "originally have",len(carbon_pos),"CP SNe matching",len(count_pos),"count"
chop_sn(carbon_pos,count_pos)
chop_sn(carbon_neg,count_neg)



print "searching for relevant SNe statistics...\n"
#for each row of stats, find the corresponding redshift and add it to an array
for stat in stats:
	for i in range(len(carbon_pos)):
		if carbon_pos[i] in stat[0]:
			get_stats(stat,name_pos,z_pos,tmax_pos)
	for i in range(len(carbon_neg)):
		if carbon_neg[i] in stat[0]:
			get_stats(stat,name_neg,z_neg,tmax_neg)

###################################

# Weight the spectra
"""
for i in range(len(data_pos)):
	s = np.multiply(data_pos[i][:,1],data_pos[i][:,2])
	e = data_pos[i][:,1]
	w = np.divide(s,e)
	data_pos[i][:,1] = w
	
for i in range(len(data_neg)):
	s2 = np.multiply(data_neg[i][:,1],data_neg[i][:,2])
	e2 = data_neg[i][:,1]
	w2 = np.divide(s2,e2)
	data_neg[i][:,1] = w2
"""
"""
def weight(d):
	for i in range(len(d)):
		s = np.multiply(d[i][:,1],d[i][:,2])
		e = d[i][:,1]
		w = np.divide(s,e)
		d[i][:,1] = w
		print d[i][:,1]
"""
#weight(data_pos)
#weight(data_neg)	

	
###################################
#de-redshifting data

for i in range(len(data_pos)):
	dez(data_pos[i],path_pos[i],carbon_pos,z_pos)
print "De-reddening C + complete."

print len(data_neg),len(z_neg),len(carbon_neg)
for i in range(len(data_neg)):
	dez(data_neg[i],path_neg[i],carbon_neg,z_neg)
print "De-reddening C - complete."


#print "\nhighest minimum:",wave_min
#print "lowest maximum:",wave_max

print "\nCreating linspace from",wave_min,"to",wave_max,"with",(wave_max-wave_min),"points"
#creates space for wavelengths given min/max parameters
wavelengths = np.linspace(wave_min,wave_max,(wave_max-wave_min+1))

print "plots will have data at these wavelength values:",wavelengths

#Interpolating
interpolate(data_pos, fit_flux_pos)
interpolate(data_neg, fit_flux_neg)

avg_flux_p = sum(fit_flux_pos)/len(data_pos)
avg_flux_n = sum(fit_flux_neg)/len(data_neg)


#residual and RMS 
res_flux_p = []
res_flux_n = []

for i in range(len(data_pos)):
	res_flux_p.append(fit_flux_pos[i]-avg_flux_p)

rms_p     = np.sqrt(np.mean(np.square(res_flux_p),axis = 0))
plus_p    = avg_flux_p + rms_p
minus_p   = avg_flux_p - rms_p
scatter_p = np.divide(rms_p,avg_flux_p)


for i in range(len(data_neg)):
	res_flux_n.append(fit_flux_neg[i]-avg_flux_n)

rms_n     = np.sqrt(np.mean(np.square(res_flux_n),axis = 0))
plus_n    = avg_flux_n + rms_n
minus_n   = avg_flux_n - rms_n
scatter_n = np.divide(rms_n,avg_flux_n)

print scatter_p

#Write data to file
datafile = open('Carbon_Results.txt','w')
for i in range(len(data_pos)):
	row = str(wavelengths[i]) + '\t' + str(avg_flux_p[i]) + '\t' + str(plus_p[i]) + '\t' + str(minus_n[i]) + '\n'
	datafile.write(row)
datafile.close()


"""
#RMS function needs lots of debugging.
#Will work on this on a later date. 
def rms(d,avg,r_f,f_f,p,n,s):
	for i in range(len(d)):
		r_f.append(f_f[i]-avg)
	rms = np.sqrt(np.mean(np.square(r_f),axis = 0))
	p = avg + rms
	print type(p), len(p)
	n = avg - rms
	s = np.divide(rms,avg)

rms(data_pos,avg_flux_p,res_flux_p,fit_flux_pos,plus_p,minus_p,scatter_p)
rms(data_neg,avg_flux_n,res_flux_n,fit_flux_neg,plus_n,minus_n,scatter_n)
print type(plus_p), len(plus_p)
plus_p = np.asarray(plus_p)
print type(plus_p), len(plus_p)
minus_p = np.array(minus_p)
"""

#print len(plus_p), len(wavelengths)

############################
print "===================="
print "calculations finished"
print "plotting..."
print "===================="
############################


#'#'#'#'#
title    = "Gram Test"
x_axis   = "Wavelength ($\AA$)"
y_axis   = "Scaled Flux"
legend1  = "Red and Unscaled"
legend2  = "Error +"
legend3  = "Error -"
x_range  = wavelengths
y_range  = avg_flux_p
Plot_2.main(x_range, y_range ,title, x_axis, y_axis, legend1, legend2, legend3)
#'#'#'#'#
"""
### file location  ### name[0], Data File: X_range , Y_range , Error_Plus , Error_minus
Carbon_Results.txt
#../../data/cfa/sn2006E/sn2006E-20060126.53-fast.flm
### Name of plot   ### name[1]
Carbon_Test.png
### title for plot ### name[2]
Test Plot
### x axis label   ### name[3]
$\mathtt{Wavelength}$ ($\AA$)
### y axis label   ### name[4]
Scaled Flux
### legend label   ### name[5], name[6], name[7]
Red and Unscaled
Error +
Error -

"""

#fig = plt.figure()
"""
#top plot containing the Composite spectrum +/- the RMS
top = fig.add_subplot(211)
plot1 = top.plot(wavelengths,avg_flux_p,'k')
plot2 = top.plot(wavelengths,plus_p,'b')
plot3 = top.plot(wavelengths,minus_p,'r')
plt.title('Composite Spectrum +/- RMS spectrum')
plt.ylabel('Relative Flux')
plt.legend([plot1,plot2,plot3],('Composite Flux','+RMS','-RMS'),'upper right',numpoints=1)

#bottom plot of the Residual
bottom = fig.add_subplot(212)
bottom.plot(wavelengths,scatter_p, 'k')
plt.title('Residual')
plt.xlabel('Wavelength')
plt.ylabel('Residual RMS')
plt.savefig('plots/Carbon_pos.png')
plt.show()
"""
"""

#labeling plot1 carbon +
top = fig.add_subplot(3,1,1)
plot1 = top.plot(wavelengths,avg_flux_p,'k')
plot2 = top.plot(wavelengths,plus_p,'b')
plot3 = top.plot(wavelengths,minus_p,'r')
plt.title('Composite Spectrum: Carbon +')
plt.ylabel('Relative Flux')
plt.legend([plot1,plot2,plot3],('Carbon +','+ RMS','- RMS'),'upper right',numpoints=1)

#labeling plot2 carbon -
mid = fig.add_subplot(3,1,2)
plot4 = mid.plot(wavelengths,avg_flux_n,'k')
plot5 = mid.plot(wavelengths,plus_n,'b')
plot6 = mid.plot(wavelengths,minus_n,'r')
plt.title('Composite Spectrum: Carbon -')
plt.ylabel('Relative Flux')
plt.legend([plot4,plot5,plot6],('Carbon - ','+ RMS','- RMS'),'upper right',numpoints=1)

#labeling plot3 residual
bottom = fig.add_subplot(3,1,3)
plot7 = bottom.plot(wavelengths,scatter_p, 'k')
plot8 = bottom.plot(wavelengths,scatter_n, 'b')
plt.title('Residual')
plt.xlabel('Wavelength')
plt.ylabel('Residual RMS')
plt.legend([plot7,plot8],('Carbon +','Carbon -'),'upper right',numpoints=1)

plt.savefig('Composite_Spectra.png')
plt.show()
"""
"""

#labeling plot1 carbon +
top = fig.add_subplot(2,1,1)
plot1 = top.plot(wavelengths,avg_flux_p,'k')
plot4 = top.plot(wavelengths,avg_flux_n,'b')
plt.title('Composite Spectrum: Carbon +')
plt.ylabel('Relative Flux')
plt.legend([plot1,plot4],('Carbon +','Carbon -'),'upper right',numpoints=1)

#labeling plot3 residual
bottom = fig.add_subplot(2,1,2)
plot7 = bottom.plot(wavelengths,scatter_p, 'k')
plot8 = bottom.plot(wavelengths,scatter_n, 'b')
plt.title('Residual')
plt.xlabel('Wavelength')
plt.ylabel('Residual RMS')
plt.legend([plot7,plot8],('Carbon +','Carbon -'),'upper right',numpoints=1)

plt.show()
"""
