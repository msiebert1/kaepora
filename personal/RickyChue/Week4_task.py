# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:17:39 2014

@author: Brian, Ricky, Lunan
"""

import numpy as np
import os, math
import matplotlib.pyplot as plt

rootdir = '/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/data/cfa/' # directory where data is


filename = []   # filename of data set
z = []          # redshift for particular SN

redshift = np.genfromtxt(rootdir + "cfasnIa_param.dat", dtype = None, skiprows = 44, unpack = True)   # loads redshift and other SN data
julian = np.genfromtxt(rootdir + "cfasnIa_mjdspec.dat", dtype = None, skiprows = 1, unpack = True)    # loads SN and julian date of max brightness
# 2603

maxbright = []

# Extracts name of SN and Julian Date of max brightness. (271)
for m in range(len(redshift)):
    if redshift[m][2] != 99999.9:
        maxbright.append([redshift[m][0], redshift[m][2], redshift[m][1]]) # SN name, julian date, z

# Finds the spectra at max brightness.
date_dif = []
for m in range(len(julian)):
    for n in range(len(maxbright)):
        if maxbright[n][0] in julian[m][0]:
            date_dif.append([maxbright[n][0], julian[m][0], math.fabs(julian[m][1] - maxbright[n][1]), maxbright[n][2]]) # SN name, SN data file name, time from max bright, z
            continue

"""
for m in range(len(files) - 1):
    if files[m][0] == files[m + 1][0]:
        min_time = min(files[m][2], files[m + 1][2])
        print files[m][0], files[m][1], min_time
"""

# finds the unique names of SNe 
SN = []
for m in range(len(date_dif)):
    SN.append(date_dif[m][0])
SN = list(set(SN))

# finds the time closest to maximum brightness
time = []
min_time = []
for m in range(len(SN)):
    for n in range(len(date_dif)):
        if SN[m] in date_dif[n][0]:
            time.append(date_dif[n][2])
            continue
    min_time.append([SN[m], min(time)])
    time = []
    
# print min_time

# associates the max brightness with the file name
max_files = []
for m in range(len(min_time)):
    for n in range(len(date_dif)):
        if min_time[m][0] == date_dif[n][0] and min_time[m][1] == date_dif[n][2] and min_time[m][1] <= 5:
            max_files.append([min_time[m][0], date_dif[n][1], date_dif[n][3]]) # SN name, data file name, z
            continue
        
SNdata = []     # top-level list of all data components
data = []       # wavelength, flux data sets

j = 0
for m in range(len(max_files)):
    print rootdir + "sn" + max_files[m][0] + "/" + max_files[m][1]
    data = np.genfromtxt(rootdir + "sn" + max_files[m][0] + "/" + max_files[m][1], skiprows = 0, unpack = True)
    if data[2][2] != 0:
        SNdata.append([max_files[m][0], max_files[m][2], data]) # fills SN data list with individual SN data set: SN name, z, data
        continue
    j = j + 1

print j

# SNdata[SN set][0=filename, 1=redshift, 2=spectral data][0=wavelength,1=flux]
# print SNdata[0][0] # SN data file
# print SNdata[0][1] # redshift, z
# print SNdata[0][2][0][1] # wavelength

deredshift_data = []    # wavelength data adjusted for redshift

#print len(SNdata)       # prints number of data files opened

for m in range(len(SNdata)):    # performs redshift adjustment on wavelength data
    deredshift_data.append((SNdata[m][2][0])/(1 + SNdata[m][1]))
#    print SNdata[m][0], deredshift_data[m]

#print len(deredshift_data)

"""
wave_min = []   # collection of the minimum wavelengths in each SN data set
wave_max = []   # collection of the maximum wavelengths in each SN data set

# finds the overlap of the data sets
for m in range(len(deredshift_data)):       # finds max and min wavelengths in each data set
    wave_min.append(min(deredshift_data[m]))
    wave_max.append(max(deredshift_data[m]))
"""
v_range = 300 # km/s
xmin = 5890*(1 - v_range/3e5) #5884    # 5890A - 300 km/s
xmax = 5896*(1 + v_range/3e5) #5902    # 5896A + 300 km/s
xstep = xmax - xmin
tstep = 10000

xmin2 = 3500
xmax2 = 7500
tstep2 = 100000

x_val = np.linspace(xmin, xmax, tstep) # creates linearly evenly spaced wavelength values within overlap
x_val2= np.linspace(xmin2, xmax2, tstep2)

#print x_val

interp_data = []    # interpolated data points for operlap
interp_data2 = []

for m in range(len(deredshift_data)):   # creates interpolated data points for flux at the overlapped wavelength values
    interp_data.append(np.interp(x_val, SNdata[m][2][0], SNdata[m][2][1]))
    interp_data2.append(np.interp(x_val2, SNdata[m][2][0], SNdata[m][2][1]))
    
median = []         # median values for each interpolated flux set
scaled_data = []    # scaled flux data using median value

median2 = []
scaled_data2 = []

for m in range(len(deredshift_data)):   # finds median value and scales interpolated data to it
    median.append(np.median(interp_data[m]))
    scaled_data.append(interp_data[m]/median[m])
    median2.append(np.median(interp_data2[m]))
    scaled_data2.append(interp_data2[m]/median2[m])

"""
print "median =", median
print len(x_val)
print len(scaled_data)
print len(interp_data)
print len(SNdata)
"""

ave_var = []
for m in range(len(SNdata)):
    ave_var.append(1/(np.median(SNdata[m][2][2]))**2)

#print ave_var

comp_data = np.average(np.array(scaled_data), weights = np.array(ave_var), axis = 0)  # composite data set, takes mean values of interpolated flux data at each wavelength
comp_data2 = np.average(np.array(scaled_data2), weights = np.array(ave_var), axis = 0) 

delta_data = []     # difference from the mean value
delta_data2 = []
for m in range(len(scaled_data)):   # finds difference between composite data and interpolated data sets
    delta_data.append(comp_data-scaled_data[m]) 
for m in range(len(scaled_data2)):   # finds difference between composite data and interpolated data sets
    delta_data2.append(comp_data2-scaled_data2[m]) 


rms_data = np.sqrt(np.mean(np.square(delta_data), axis = 0))    # finds root-mean-square of differences at each wavelength within overlap
rms_data2 = np.sqrt(np.mean(np.square(delta_data2), axis = 0))   

unity = []
for m in range(len(scaled_data)):
    unity.append(np.interp(x_val, [x_val[0], x_val[tstep - 1]], [scaled_data[m][0], scaled_data[m][tstep - 1]]))

EW_term = []
for m in range(len(scaled_data)):
    EW_term.append((1-scaled_data[m]/unity[m])*(xstep/tstep))

EW = []
for m in range(len(EW_term)):
    EW.append(sum(EW_term[m])) 

#print EW_SN_dat_wt[0][3]

EW_lo = []
EW_hi = []
scaled_lo = []
scaled_hi = []
SN_hi = []
SN_lo = []
var_lo = []
var_hi = []

EW_lo2 = []
EW_hi2 = []
scaled_lo2 = []
scaled_hi2 = []
SN_hi2 = []
SN_lo2 = []
var_lo2 = []
var_hi2 = []

for m in range(len(EW)):
    if EW[m] <= np.median(EW):
        EW_lo.append(EW[m])
        scaled_lo.append(scaled_data[m])
        SN_lo.append(SNdata[m][0])
        var_lo.append(ave_var[m])
        EW_lo2.append(EW[m])
        scaled_lo2.append(scaled_data2[m])
        SN_lo2.append(SNdata[m][0])
        var_lo2.append(ave_var[m])
    else:
        EW_hi.append(EW[m])
        scaled_hi.append(scaled_data[m])
        SN_hi.append(SNdata[m][0])
        var_hi.append(ave_var[m])
        EW_hi2.append(EW[m])
        scaled_hi2.append(scaled_data2[m])
        SN_hi2.append(SNdata[m][0])
        var_hi2.append(ave_var[m])

lo_data = np.average(np.array(scaled_lo), weights = np.array(var_lo), axis = 0)  # composite data set, takes mean values of interpolated flux data at each wavelength
hi_data = np.average(np.array(scaled_hi), weights = np.array(var_hi), axis = 0)  # composite data set, takes mean values of interpolated flux data at each wavelength

lo_data2 = np.average(np.array(scaled_lo2), weights = np.array(var_lo2), axis = 0)
hi_data2 = np.average(np.array(scaled_hi2), weights = np.array(var_hi2), axis = 0) 

comp_unity = np.interp(x_val, [x_val[0], x_val[tstep - 1]], [comp_data[0], comp_data[tstep - 1]])
comp_EW_term = ((1-comp_data/comp_unity)*(xstep/tstep))

# print EW_term

# plots composite and rms data set
plt.subplot(2, 1, 1)
plt.title('Composite Spectra')
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
plt.plot(x_val2, comp_data2, label = "Composite")
plt.plot(x_val2, hi_data2, label = "Composite Hi")
plt.plot(x_val2, lo_data2, label = "Composite Lo")
plt.plot(x_val2, comp_data2+rms_data2, label = "+ RMS")
plt.plot(x_val2, comp_data2-rms_data2, label = "- RMS")
plt.legend()

plt.subplot(2, 1, 2)
plt.title('Composite Spectra of Na doublet')
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
plt.plot(x_val, comp_data, label = "Composite")
plt.plot(x_val, hi_data, label = "Composite Hi")
plt.plot(x_val, lo_data, label = "Composite Lo")
plt.plot(x_val, comp_data+rms_data, label = "+ RMS")
plt.plot(x_val, comp_data-rms_data, label = "- RMS")
plt.legend()

plt.savefig('/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/RickyChue/W4_plot.png')

"""
plt.subplot(2, 1, 2)
plt.title('EW terms')
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
plt.plot(x_val, comp_EW_term, label = "EW")
#plt.legend()
"""


"""
# plots residual rms of data
plt.subplot(2, 1, 2)
plt.title('Residuals')
plt.xlabel('wavelength (A)')
plt.ylabel('residual RMS')
plt.plot(x_val, rms_data/comp_data*100, label = "residual rms")
#plt.savefig('/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/RickyChue/Week3_spectrum.png')
#plt.legend()
"""

# other plots:
"""
# plots de-redshifted, interpolated, unscaled data within overlap
plt.subplot(1, 4, 1)
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
for m in range(len(deredshift_data)):
    plt.plot(x_val, interp_data[m], label = SNdata[m][0])
#plt.legend()

# plots de-redshifted, interpolated, scaled data within overlap
plt.subplot(1, 4, 2)
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
for m in range(len(deredshift_data)):
    plt.plot(x_val, scaled_data[m], label = SNdata[m][0])
#plt.legend()

# plots composite data set
plt.subplot(1, 4, 3)
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
plt.plot(x_val, comp_data, label = "Composite")
plt.legend()

# plots difference from composite
plt.subplot(1, 4, 4)
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
for m in range(len(deredshift_data)):
    plt.plot(x_val, delta_data[m], label = SNdata[m][0])
#plt.legend()

# plots rms of data sets from composite
plt.subplot(1, 4, 4)
plt.xlabel('wavelength (A)')
plt.ylabel('rms/comp %')
plt.plot(x_val, rms_data/comp_data*100, label = "rms")
#plt.legend()

"""









# print redshift[0]
# print redshift
# print len(redshift)


"""

for m in range(len(min_time)):
    for n in range(len(date_dif)):
        for root, dirs, files in os.walk(rootdir):  # looks in all folders within rootdir
            for file in files:
                if min_time[m][0] == date_dif[n][0] and min_time[m][1] == date_dif[n][2] and file == date_dif[n][1]:       # loads files that matches SN and max brightness time
                    data = np.genfromtxt(os.path.join(root, file), skiprows = 0, unpack = True)    # opens data file
                    filename = os.path.splitext(os.path.basename(os.path.join(root, file)))[0]     # takes full directory name of file, removes directory folder names, then removes .flm extension
                    SNname = os.path.basename(root)    # name of SN, used to match redshift data with SN
                    print filename                     # displays name of SN, useful for seeing last good data read in case of a corrupted data file
                    for m in range(len(redshift)):     # matches SN name with associated redshift
                    #                 print SNname[2:], redshift[m][0]
                        if SNname[2:] == redshift[m][0]:
                             z = redshift[m][1]
                             #                     print z
                             SNdata.append([filename, z, data]) # fills SN data list with individual SN data set
                             print SNdata

"""
















# older coding blocks:

"""
freq1 = d1[0]
flux1 = d1[1]
""
SNdata = [filename,freq1,flux1]
print SNdata[0]
print SNdata[1]
print SNdata[2]

""
for files in os.walk(rootdir):
#    print files
    for file in files[:]:
        print file
        if file.endswith(".flm"):
             print file
"""


"""

for file in os.listdir(rootdir):
    if file.endswith(".flm"):
        print file

"""

"""
deredshift_wave = []

#    for n in range(len(SNdata[m][2][0])):
#        deredshift_wave.append((SNdata[m][2][0][n])/(1 + SNdata[m][1]))
"""
