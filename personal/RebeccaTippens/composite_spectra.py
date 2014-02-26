from __future__ import division #Enables new-style division
from math import floor, ceil
from random import sample
import matplotlib.pyplot as plt
import numpy as np 
import sqlite3 as sq3
import os

#Connect to the database
path = "/Users/Rebecca/astr596/personal/MichaelSchubert/SNe.db"
con = sq3.connect(path)

#Creates a cursor object to execute commands
cur = con.cursor()

#Returns the 10 SNe with the highest redshifts
cur.execute("SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 40")

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

integer_wave = np.arange(ceil(min_wave), floor(max_wave), dtype=int, step=2)

table_for_composite = []

for key in SN_data.keys():
    deredshifted_wave = SN_data[key][1]
    flux = SN_data[key][2]
    interpolated_flux = np.interp(integer_wave, deredshifted_wave, flux)
    median_flux = np.median(interpolated_flux)
    scaled_flux = interpolated_flux/median_flux
    SN_data[key].append(integer_wave)
    SN_data[key].append(scaled_flux)
    table_for_composite.append(scaled_flux)

composite_flux = np.mean(table_for_composite, axis = 0) #Averages flux values to make composite

dist_from_mean_table = []

for key in SN_data.keys():
    interpolated_flux = SN_data[key][4]
    dist_from_mean = composite_flux - interpolated_flux
    SN_data[key].append(dist_from_mean)
    dist_from_mean_table.append(dist_from_mean)

RMS = np.sqrt(np.mean(np.square(dist_from_mean_table), axis = 0))

np.savetxt('composite_spectra.flm', np.transpose([integer_wave, composite_flux]), fmt="%d %26.18e")

#Plots composite spectra
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(integer_wave, composite_flux, color="#9acd32", label = "Composite")
ax1.plot(integer_wave, (composite_flux + RMS), color="#ff6600", label = "+ RMS")
ax1.plot(integer_wave, (composite_flux - RMS), color="c", label = "- RMS")
ax1.set_xlim(min(integer_wave), max(integer_wave))
ax1.set_xlabel("Rest Wavelength " + "(" + u"\u212B" + ")") #need to replace unicode?
ax1.set_ylabel("Relative Flux")
ax1.legend(loc = "best")

#Legend Options ---
#upper right: loc=1
#upper left: loc=2
#lower left: loc=3
#lower right: loc=4
#right: loc=5
#center left: loc=6
#center right: loc=7
#lower center: loc=8
#upper center: loc=9
#center: loc=10

#Plots 
ax2 = fig.add_subplot(212)
ax2.plot(integer_wave, (RMS/composite_flux)*100, color="m", label = "Residual RMS")
ax2.set_xlim(min(integer_wave), max(integer_wave))
ax2.set_ylabel("Residual RMS")
ax2.legend(loc = "best")
fig.set_tight_layout(True)
plt.savefig("composite_spectra.pdf", format="PDF")
plt.show()