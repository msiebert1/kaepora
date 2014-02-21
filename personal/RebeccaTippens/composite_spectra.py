from __future__ import division #Future statement to use New-style division
from math import floor, ceil
from random import sample
import matplotlib.pyplot as plt
import os #Uses operating system-dependent functionality
import glob #Enables filename pattern matching
import numpy as np 
import scipy as sy 
import sqlite3 as sq3
import time

"""
#Declares path to root directory
root = "/Users/Rebecca/astr596/data"
"""

#Connect to the database
path = "/Users/Rebecca/astr596/personal/MichaelSchubert/SNe.db"
con = sq3.connect(path)

#Creates a cursor object to execute commands
cur = con.cursor()

#Returns the 10 SNe with the highest redshifts
cur.execute("SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 30")

root = "/Users/Rebecca/astr596/data/cfa/"

SN_data = {} #Creates empty dictionary for SN data

#Keeps track of how many spectra load, how many don't
good = 0 #Initializes total spectra count to zero
bad = 0 #Initializes bad spectra count to zero
bad_files = []

min_waves = [] #Creates empty array for minimum deredshifted wavelength from each file
max_waves = [] #Creates empty array for maximum deredshifted wavelength from each file

#Puts filenames in list
for row in cur:
    file_name = row[0]
    SN_name = "sn" + row[1]
    file_path = os.path.join(root, SN_name, file_name)
    z = row[2]
    try: #Makes sure data loads, deredshifting and minmax functions work properly
        wave, flux = np.loadtxt(file_path, usecols=(0,1), unpack=True)
        deredshifted_wave = wave/(1+z)
        SN_data[file_name] = [z, deredshifted_wave, flux]
        min_waves.append(min(deredshifted_wave))
        max_waves.append(max(deredshifted_wave))
        good += 1 #Counts number of files for which everything executed properly
    except: #Stores number and names of files that didn't load correctly
        bad +=1
        bad_files.append(file_name)

if bad == 0:
    print "All files loaded successfully." #Returns this message if all files loaded correctly
else:
    print str(good) + " files read successfully." #Says how many files loaded correctly (if not all)
    print "The following " + str(bad) + " file path(s) produced load errors:" #Says how many files produced load errors
    for item in bad_files:
        print item #Prints of names of files which produced load errors

min_wave = min(min_waves) #Finds overall min deredshifted wavelength
max_wave = max(max_waves) #Finds overall max deredshifted wavelength

integer_wave = np.arange(ceil(min_wave), floor(max_wave), dtype=int, step=1)

table_for_composite = []

for key in SN_data.keys():
    deredshifted_wave = SN_data[key][1]
    flux = SN_data[key][2]
    interpolated_flux = np.interp(integer_wave, deredshifted_wave, flux)
    SN_data[key].append(integer_wave)
    SN_data[key].append(interpolated_flux)
    table_for_composite.append(interpolated_flux)

composite_flux = np.mean(table_for_composite, axis = 0) #Averages flux values to make composite

dist_from_mean_table = []

for key in SN_data.keys():
    interpolated_flux = SN_data[key][4]
    dist_from_mean = composite_flux - interpolated_flux
    SN_data[key].append(dist_from_mean)
    dist_from_mean_table.append(dist_from_mean)

RMS = np.sqrt(np.mean(np.square(dist_from_mean_table), axis = 0))

#Plots composite spectra
plt.subplot(211)
plt.plot(integer_wave, composite_flux, color="c", label = "Composite")
plt.plot(integer_wave, (composite_flux + RMS), color="m", label = "+ RMS")
plt.plot(integer_wave, (composite_flux - RMS), color="#ff6600", label = "- RMS")
plt.yscale("log")
plt.xlabel("Wavelength " + "(" + u"\u212B" + ")") #need to replace unicode?
plt.ylabel("log (Flux)")
plt.legend(loc=4)

#Plots residuals
plt.subplot(212)
plt.plot(integer_wave, (RMS/composite_flux)*100, label = "Residual RMS")
plt.xlabel("Wavelength " + "(" + u"\u212B" + ")") #need to replace unicode?
plt.ylabel("Residual RMS")
plt.savefig("composite_spectra.pdf", format="PDF")
plt.show()