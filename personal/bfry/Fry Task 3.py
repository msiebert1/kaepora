# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:17:39 2014

@author: Brian and Ricky
"""

import numpy as np
import os
import matplotlib.pyplot as plt

rootdir = './data/cfa' # directory where data is

SNdata = []     # top-level list of all data components
data = []       # wavelength, flux data sets
filename = []   # filename of data set
z = []          # redshift for particular SN

redshift = np.genfromtxt("cfasnIa_param.dat", dtype = None, skiprows = 44, unpack = True)   # loads redshift and other SN data

# print redshift[0]
# print redshift
# print len(redshift)

for root, dirs, files in os.walk(rootdir):  # looks in all folders within rootdir
    for file in files:
        if file.endswith(".flm"):       # loads files with .flm extensions
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
#             print SNdata

# SNdata[SN set][0=filename, 1=redshift, 2=spectral data][0=wavelength,1=flux]
# print SNdata[0][0] # SN data file
# print SNdata[0][1] # redshift, z
# print SNdata[0][2][0][1] # wavelength

deredshift_data = []    # wavelength data adjusted for redshift

print len(SNdata)       # prints number of data files opened

for m in range(len(SNdata)):    # performs redshift adjustment on wavelength data
    deredshift_data.append((SNdata[m][2][0])/(1 + SNdata[m][1]))
#    print SNdata[m][0], deredshift_data[m]

# print len(deredshift_data)

wave_min = []   # collection of the minimum wavelengths in each SN data set
wave_max = []   # collection of the maximum wavelengths in each SN data set

# finds the overlap of the data sets
for m in range(len(deredshift_data)):       # finds max and min wavelengths in each data set
    wave_min.append(min(deredshift_data[m]))
    wave_max.append(max(deredshift_data[m]))

xmin = max(wave_min)    # finds the maximum minimum value (lower wavelength bound)
xmax = min(wave_max)    # finds the minimum maximum value (upper wavelength bound)

x_val = np.linspace(xmin, xmax, 100000) # creates linearly evenly spaced wavelength values within overlap

interp_data = []    # interpolated data points for operlap

for m in range(len(deredshift_data)):   # creates interpolated data points for flux at the overlapped wavelength values
    interp_data.append(np.interp(x_val, SNdata[m][2][0], SNdata[m][2][1]))

median = []         # median values for each interpolated flux set
scaled_data = []    # scaled flux data using median value

for m in range(len(deredshift_data)):   # finds median value and scales interpolated data to it
    median.append(np.median(interp_data[m]))
    scaled_data.append(interp_data[m]/median[m])

"""
print "median =", median
print len(x_val)
print len(scaled_data)
print len(interp_data)
print len(SNdata)
"""

comp_data = np.mean(scaled_data, axis = 0)  # composite data set, takes mean values of interpolated flux data at each wavelength

delta_data = []     # difference from the mean value

for m in range(len(scaled_data)):   # finds difference between composite data and interpolated data sets
    delta_data.append(comp_data-scaled_data[m]) 

rms_data = np.sqrt(np.mean(np.square(delta_data), axis = 0))    # finds root-mean-square of differences at each wavelength within overlap

# plots composite and rms data set
plt.subplot(2, 1, 1)
plt.title('Composite Spectra')
plt.xlabel('wavelength (A)')
plt.ylabel('relative flux')
plt.plot(x_val, comp_data, label = "Composite")
plt.plot(x_val, comp_data+rms_data, label = "+ RMS")
plt.plot(x_val, comp_data-rms_data, label = "- RMS")
plt.legend()

# plots residual rms of data
plt.subplot(2, 1, 2)
plt.title('Residuals')
plt.xlabel('wavelength (A)')
plt.ylabel('residual RMS')
plt.plot(x_val, rms_data/comp_data*100, label = "residual rms")
#plt.legend()



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
