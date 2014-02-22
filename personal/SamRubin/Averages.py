#!/usr/bin/env python -i
'''
Created on Feb 14, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table

class supernova(object):
    """Attributes can be added"""

SN_Array = []
compare_spectrum = []

#names = np.loadtxt("", usecols = (0)) #"need to grab all the other data too"
#for row in names:
#    SN = supernova()
#    SN.input = row[0]
#    SN_Array.append(SN)
#print len(names), #"supernovae found"
file_list = []
file_list = glob.glob("../../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../SamRubin/MaxSpectra.dat", dtype = 'str')

delta_t = Table.read('../AaronBeaudoin/week5/SpectraTime.dat',format='ascii')
file_name = delta_t["File"]
ages = delta_t["delta_t"]
sn_ages=[]
for i in range(len(max_light)):
	age=0
	for j in range(len(file_name)):
		if max_light[i] in file_name[j]:
			age = ages[j]
	sn_ages.append(age)	
#selected correct ones for averaging
# for i in range(len(files)):
# 	file=files[i]
# 	file_list.append(file["col2"])

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading supernovae from database..."
names = []
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = supernova()
    if row[2] != None:
		SN.filename = row[0]
		SN.name = row[1]
		SN.redshift = row[2]
		SN.minwave = row[3]
		SN.maxwave = row[4]
		SN.age = sn_ages
		SN_Array.append(SN)
		names.append(SN.name)

print len(SN_Array), "items found"
print len(sn_ages)

names,temp = np.unique(names,return_index=True)
SN_Array =  [SN_Array[i] for i in temp]
print len(SN_Array), "unique spectra found"
# rand = np.random.randint(0,len(SN_Array),10)
# all_array = SN_Array
# SN_Array=[]
# for i in rand:
#     tempSN = all_array[i]
#     SN_Array.append(tempSN)
#"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1

for SN in SN_Array[0:10]:
    k = 1
    for filename in max_light:   
        if SN.filename in filename:
            data = np.loadtxt(filename)
            SN.wavelength = data[:,0]
            SN.flux = data[:,1]
            SN.redshifts = np.zeros(len(data[:,1]))
            SN.redshifts.fill(SN.redshift)
            error = data[:,2]
            if not all(x == 0.0 for x in error):
                SN.error = error
            break
    print j, 'supernovae fluxed'
    j += 1
print "Done searching. Checking data..."   
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]

print "Done checking. ", len(SN_Array), "items remain"       

#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

bad_range_Array = []

#averages with weights based on the given errors in .flm files
def average(compare_spectrum,SN):
	residual = compare_spectrum.error
	avg_flux = compare_spectrum.flux
	avg_red = compare_spectrum.redshifts
	mean_flux = compare_spectrum.flux
	redshifts = compare_spectrum.redshifts
	lowindex = find_nearest(compare_spectrum.wavelength,np.min(SN.wavelength))
	highindex = find_nearest(compare_spectrum.wavelength,np.max(SN.wavelength))
#should be np.where(compare_spectrum.wavelength == np.min(SN.wavelength) if data aligned)
	for i in lowindex+np.arange(highindex-lowindex):
#        avg_flux[i] = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
		try:
			if ((SN.error[i] != 0) & (compare_spectrum.error[i] != 0)):
				fluxes = np.array([compare_spectrum.flux[i],SN.flux[i]])
				weights = np.array([1./compare_spectrum.error[i], 1./SN.error[i]])
				redshifts = np.array([compare_spectrum.redshifts[i],SN.redshifts[i]])
				mean_flux[i] = np.average(fluxes, weights=weights)
				avg_flux[i] = np.average(fluxes,weights=weights)
				avg_red[i] = np.average(redshifts, weights = weights) 
				#compare_spectrum.error[i] = math.sqrt((compare_spectrum.flux[i]-avg_flux[i])**2 + (SN.flux[i]-avg_flux[i])**2)
				compare_spectrum.redshifts[i] = np.average(redshifts, weights=weights)
			else:
				break
		except IndexError:
			print "No flux data?"
			i += 1
			break
            
	compare_spectrum.flux = avg_flux
	compare_spectrum.mean = mean_flux
	compare_spectrum.redshifts = avg_red
# Add residual formula?
	return compare_spectrum

"""scale"""
print "Scaling.."
q = 1
compare_spectrum = SN_Array[0]
#scales, averages, weighted average
for SN in SN_Array:
    low_wave = 3000
    high_wave = 7000
#     plt.plot(SN.wavelength, SN.flux,label=SN.name)
    if np.min(SN.wavelength)>= low_wave:
        low_wave = np.min(SN.wavelength)
    if np.max(SN.wavelength) <= high_wave:
        high_wave = np.max(SN.wavelength)
    temp = np.abs(SN.wavelength-low_wave)
    lowindex = np.where(temp == np.min(temp))
    temp = np.abs(SN.wavelength-high_wave)
    highindex = np.where(temp == np.min(temp))
    lowindex = lowindex[0]
    highindex = highindex[0]
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
#    plt.plot(SN.wavelength, SN.flux,label=SN.name)
    print lowindex, "low index", highindex, "high index"
    factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
    scale_factor = np.mean(factors)
    SN.flux[lowindex:highindex] *= scale_factor
    SN.error[lowindex:highindex] *= scale_factor
    print q, "items scaled"
    
    compare_spectrum = average(compare_spectrum,SN)
        
    q += 1
	
print len(compare_spectrum.wavelength), len(compare_spectrum.redshifts)
"""composite"""
plt.figure(1)
plt.subplot(311)
plt.plot(compare_spectrum.wavelength, compare_spectrum.flux,label='Weighted Composite')
plt.plot(compare_spectrum.wavelength, compare_spectrum.mean,label='Composite')
plt.xlabel('Wavelength (A)')
plt.ylabel('Scaled Flux')
plt.subplot(312)
plt.plot(compare_spectrum.wavelength, compare_spectrum.error, label='Error')
plt.ylabel('Error')
plt.subplot(313)
plt.plot(compare_spectrum.wavelength, compare_spectrum.redshifts, label='Redshift')
plt.subplot(311)
plt.legend(prop={'size':10})
plt.savefig('Average.png')
plt.show()
#!/usr/bin/env python -i
'''
Created on Feb 14, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table

class supernova(object):
    """Attributes can be added"""

SN_Array = []
compare_spectrum = []

#names = np.loadtxt("", usecols = (0)) #"need to grab all the other data too"
#for row in names:
# SN = supernova()
# SN.input = row[0]
# SN_Array.append(SN)
#print len(names), #"supernovae found"
file_list = []
file_list = glob.glob("../../data/cfa/*/*.flm")
max_light = []
max_light = np.loadtxt("../SamRubin/MaxSpectra.dat", dtype = 'str')

delta_t = Table.read('../AaronBeaudoin/week5/SpectraTime.dat',format='ascii')
file_name = delta_t["File"]
ages = delta_t["delta_t"]
sn_ages=[]
#for i in range(len(max_light)):
# age=0
# for j in range(len(file_name)):
# if max_light[i] in file_name[j]:
# age = ages[j]
# sn_ages.append(age)
#selected correct ones for averaging
# for i in range(len(files)):
# file=files[i]
# file_list.append(file["col2"])

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
print "Reading supernovae from database..."
names = []
for row in cur.execute('SELECT Filename, SN, Redshift, MinWave, MaxWave FROM Supernovae'):
    SN = supernova()
    if row[2] != None:
        SN.filename = row[0]
        SN.name = row[1]
        SN.redshift = row[2]
        SN.minwave = row[3]
        SN.maxwave = row[4]
        #SN.age = sn_ages #this appends the entire list of ages to each supernova which is totally wrong
        SN_Array.append(SN)
        names.append(SN.name)
SN_Array = SN_Array[20:]
print len(SN_Array), "items found"
print len(sn_ages)

#names,temp = np.unique(names,return_index=True)
#SN_Array = [SN_Array[i] for i in temp]
#print len(SN_Array), "unique spectra found"
# rand = np.random.randint(0,len(SN_Array),10)
# all_array = SN_Array
# SN_Array=[]
# for i in rand:
# tempSN = all_array[i]
# SN_Array.append(tempSN)
#"""Adds flux to the stuff"""
print "Matching flux data..."
j = 1

for SN in SN_Array[0:100]:
    k = 1
    for filename in max_light:
        if SN.filename in filename:
            data = np.loadtxt(filename)
            SN.wavelength = data[:,0]
            SN.flux = data[:,1]
            SN.redshifts = np.zeros(len(data[:,1]))
            SN.redshifts.fill(SN.redshift)
            SN.residual = np.zeros(len(data[:,1]))
            #SN.age = sn_ages[i]
            #somewhere here is where the ages should be assigned
            error = data[:,2]
            if not all(x == 0.0 for x in error):
                SN.error = error
            break
    print j, 'supernovae fluxed'
    j += 1
print "Done searching. Checking data..."
SN_Array = [SN for SN in SN_Array if hasattr(SN, 'error')]

print "Done checking. ", len(SN_Array), "items remain"

#gets as close as possible to matching the compare spectrum wavelength values
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

bad_range_Array = []

#averages with weights based on the given errors in .flm files
def average(compare_spectrum,SN):
        if SN == SN_Array[0]:
            compare_spectrum == SN
            return compare_spectrum
avg_flux = compare_spectrum.flux
avg_red = compare_spectrum.redshifts
mean_flux = compare_spectrum.flux
redshifts = compare_spectrum.redshifts
# ages = compare_spectrum.age
# avg_ages = compare_spectrum.age
lowindex = find_nearest(compare_spectrum.wavelength,np.min(SN.wavelength))
highindex = find_nearest(compare_spectrum.wavelength,np.max(SN.wavelength))
#should be np.where(compare_spectrum.wavelength == np.min(SN.wavelength) if data aligned)
for i in lowindex+np.arange(highindex-lowindex):
# avg_flux[i] = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
try:
if SN.error[i] != 0 and compare_spectrum.error[i] != 0:
fluxes = np.array([compare_spectrum.flux[i],SN.flux[i]])
weights = np.array([1./compare_spectrum.error[i], 1./SN.error[i]])
redshifts = np.array([compare_spectrum.redshifts[i],SN.redshifts[i]])
#ages = np.array([compare_spectrum.age[i],SN.age[i]])
                                mean_flux[i] = np.average(fluxes, weights=weights)
avg_flux[i] = np.average(fluxes,weights=weights)
avg_red[i] = np.average(redshifts, weights = weights)
compare_spectrum.redshifts[i] = np.average(redshifts, weights=weights)
# avg_ages[i] = np.average(ages, weights=weights)
                                compare_spectrum.error[i] = math.sqrt((weights[0]**2)*(compare_spectrum.error[i])**2 + (weights[1]**2)*(SN.error[i])**2)
else:
break
except IndexError:
print "No flux data?"
                        bad_range_Array.append(SN)
break
            
compare_spectrum.flux = avg_flux
compare_spectrum.mean = mean_flux
compare_spectrum.redshifts = avg_red
        #compare_spectrum.ages = avg_ages
# Add residual formula?
return compare_spectrum

"""scale"""
print "Scaling.."
q = 1
compare_spectrum = SN_Array[0]
#scales, averages, weighted average
plt.figure(1)
out_of_range_low = []
out_of_range_high = []
for SN in SN_Array:
    low_wave = 3000
    high_wave = 7000
# plt.plot(SN.wavelength, SN.flux,label=SN.name)
    if np.min(SN.wavelength)>= low_wave:
        low_wave = np.min(SN.wavelength)
    if np.max(SN.wavelength) <= high_wave:
        high_wave = np.max(SN.wavelength)
    #check the number of scaled versus the number before scaling, if != then expand range and run again
    temp = np.abs(SN.wavelength-low_wave)
    lowindex = np.where(temp == np.min(temp))
    temp = np.abs(SN.wavelength-high_wave)
    highindex = np.where(temp == np.min(temp))
    lowindex = lowindex[0]
    highindex = highindex[0]
    low = SN.wavelength[lowindex]
    high = SN.wavelength[highindex]
    SN.flux /= np.median(SN.flux)
# plt.plot(SN.wavelength, SN.flux,label=SN.name)
    print lowindex, "low index", highindex, "high index"
    factors = compare_spectrum.flux[lowindex:highindex] / SN.flux[lowindex:highindex]
    scale_factor = np.mean(factors)
    SN.flux[lowindex:highindex] *= scale_factor
    SN.error[lowindex:highindex] *= scale_factor
    plt.subplot(311)
    plt.plot(SN.wavelength, SN.flux)
    print q, "items scaled at factor", scale_factor
    compare_spectrum = average(compare_spectrum,SN)
    q += 1
print compare_spectrum.error
    
        
print len(compare_spectrum.wavelength), len(compare_spectrum.redshifts)
"""composite"""
plt.figure(1)
plt.subplot(311)
plt.plot(compare_spectrum.wavelength, compare_spectrum.flux,label='Weighted Composite')
plt.plot(compare_spectrum.wavelength, compare_spectrum.mean,label='Composite')
plt.xlabel('Wavelength (A)')
plt.ylabel('Scaled Flux')
plt.xlim(4000,7000)
plt.ylim(0,2.8)
plt.subplot(312)
plt.plot(compare_spectrum.wavelength, compare_spectrum.error, label='Error')
plt.ylabel('Error')
plt.plot(compare_spectrum.wavelength, compare_spectrum.residual, label='Residual')
plt.xlim(4000,7000)
plt.subplot(313)
plt.plot(compare_spectrum.wavelength, compare_spectrum.redshifts, label='Redshift')
plt.ylabel('Avg. Redshift')
#plt.subplot(414)
#plt.plot(compare_spectrum.wavelength, compare_spectrum.ages, label='Age')
#plt.ylabel('Age')
plt.subplot(311)
plt.legend(prop={'size':10})
plt.savefig('Average.png')
plt.show()