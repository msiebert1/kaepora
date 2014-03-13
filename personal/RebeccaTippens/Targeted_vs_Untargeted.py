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

#Sets root directory
root = "/Users/Rebecca/astr596/data/cfa/"

#--------------------------- Targeted objects code

cur.execute("SELECT * FROM Supernovae WHERE targeted = 1")

targ_SN_data = {} #Creates empty dictionary for SN data

#Keeps track of how many spectra load, how many don't
targ_good = 0 #Initializes total spectra count to zero
targ_bad = 0 #Initializes bad spectra count to zero
targ_bad_files = []

targ_min_waves = [] #Creates empty array for minimum deredshifted wavelength from each file
targ_max_waves = [] #Creates empty array for maximum deredshifted wavelength from each file

#Puts filenames in list
for row in cur:
    targ_file_name = row[0]
    targ_SN_name = "sn" + row[1]
    targ_file_path = os.path.join(root, targ_SN_name, targ_file_name)
    targ_z = row[2]
    try: #Makes sure data loads, deredshifting and minmax functions work properly
        targ_wave, targ_flux = np.loadtxt(targ_file_path, usecols=(0,1), unpack=True)
        targ_deredshifted_wave = targ_wave/(1+targ_z)
        targ_SN_data[targ_file_name] = [targ_z, targ_deredshifted_wave, targ_flux]
        targ_min_waves.append(min(targ_deredshifted_wave))
        targ_max_waves.append(max(targ_deredshifted_wave))
        targ_good += 1 #Counts number of files for which everything executed properly
    except: #Stores number and names of files that didn't load correctly
        targ_bad +=1
        targ_bad_files.append(targ_file_name)

if targ_bad == 0:
    print "All " + str(targ_good) + " targeted SN files loaded successfully." #Returns this message if all files loaded correctly
else:
    print str(targ_good) + " targeted SN files read successfully." #Says how many files loaded correctly (if not all)
    print "The following " + str(targ_bad) + " targeted SN file path(s) produced load errors:" #Says how many files produced load errors
    for item in targ_bad_files:
        print item #Prints of names of files which produced load errors

targ_min_wave = min(targ_min_waves) #Finds overall min deredshifted wavelength
targ_max_wave = max(targ_max_waves) #Finds overall max deredshifted wavelength

targ_integer_wave = np.arange(ceil(targ_min_wave), floor(targ_max_wave), dtype=int, step=2)

targ_table_for_composite = []

for key in targ_SN_data.keys():
    targ_deredshifted_wave = targ_SN_data[key][1]
    targ_flux = targ_SN_data[key][2]
    targ_interpolated_flux = np.interp(targ_integer_wave, targ_deredshifted_wave, targ_flux)
    targ_median_flux = np.median(targ_interpolated_flux)
    targ_scaled_flux = targ_interpolated_flux/targ_median_flux
    targ_SN_data[key].append(targ_integer_wave)
    targ_SN_data[key].append(targ_scaled_flux)
    targ_table_for_composite.append(targ_scaled_flux)

targ_composite_flux = np.mean(targ_table_for_composite, axis = 0) #Averages flux values to make composite

targ_dist_from_mean_table = []

for key in targ_SN_data.keys():
    targ_interpolated_flux = targ_SN_data[key][4]
    targ_dist_from_mean = targ_composite_flux - targ_interpolated_flux
    targ_SN_data[key].append(targ_dist_from_mean)
    targ_dist_from_mean_table.append(targ_dist_from_mean)

targ_RMS = np.sqrt(np.mean(np.square(targ_dist_from_mean_table), axis = 0))

np.savetxt('Targeted.flm', np.transpose([targ_integer_wave, targ_composite_flux]), fmt="%d %26.18e")

#--------------------------- Untargeted objects code

cur.execute("SELECT * FROM Supernovae WHERE targeted = 0")

untarg_SN_data = {} #Creates empty dictionary for SN data

#Keeps track of how many spectra load, how many don't
untarg_good = 0 #Initializes total spectra count to zero
untarg_bad = 0 #Initializes bad spectra count to zero
untarg_bad_files = []

untarg_min_waves = [] #Creates empty array for minimum deredshifted wavelength from each file
untarg_max_waves = [] #Creates empty array for maximum deredshifted wavelength from each file

#Puts filenames in list
for row in cur:
    untarg_file_name = row[0]
    untarg_SN_name = "sn" + row[1]
    untarg_file_path = os.path.join(root, untarg_SN_name, untarg_file_name)
    untarg_z = row[2]
    try: #Makes sure data loads, deredshifting and minmax functions work properly
        untarg_wave, untarg_flux = np.loadtxt(untarg_file_path, usecols=(0,1), unpack=True)
        untarg_deredshifted_wave = untarg_wave/(1+untarg_z)
        untarg_SN_data[untarg_file_name] = [untarg_z, untarg_deredshifted_wave, untarg_flux]
        untarg_min_waves.append(min(untarg_deredshifted_wave))
        untarg_max_waves.append(max(untarg_deredshifted_wave))
        untarg_good += 1 #Counts number of files for which everything executed properly
    except: #Stores number and names of files that didn't load correctly
        untarg_bad +=1
        untarg_bad_files.append(untarg_file_name)

if untarg_bad == 0:
    print "All " + str(untarg_good) + " untargeted SN files loaded successfully." #Returns this message if all files loaded correctly
else:
    print str(untarg_good) + " untargeted SN files read successfully." #Says how many files loaded correctly (if not all)
    print "The following " + str(untarg_bad) + " untargeted SN file path(s) produced load errors:" #Says how many files produced load errors
    for item in untarg_bad_files:
        print item #Prints of names of files which produced load errors

untarg_min_wave = min(untarg_min_waves) #Finds overall min deredshifted wavelength
untarg_max_wave = max(untarg_max_waves) #Finds overall max deredshifted wavelength

untarg_integer_wave = np.arange(ceil(untarg_min_wave), floor(untarg_max_wave), dtype=int, step=2)

untarg_table_for_composite = []

for key in untarg_SN_data.keys():
    untarg_deredshifted_wave = untarg_SN_data[key][1]
    untarg_flux = untarg_SN_data[key][2]
    untarg_interpolated_flux = np.interp(untarg_integer_wave, untarg_deredshifted_wave, untarg_flux)
    untarg_median_flux = np.median(untarg_interpolated_flux)
    untarg_scaled_flux = untarg_interpolated_flux/untarg_median_flux
    untarg_SN_data[key].append(untarg_integer_wave)
    untarg_SN_data[key].append(untarg_scaled_flux)
    untarg_table_for_composite.append(untarg_scaled_flux)

untarg_composite_flux = np.mean(untarg_table_for_composite, axis = 0) #Averages flux values to make composite

untarg_dist_from_mean_table = []

for key in untarg_SN_data.keys():
    untarg_interpolated_flux = untarg_SN_data[key][4]
    untarg_dist_from_mean = untarg_composite_flux - untarg_interpolated_flux
    untarg_SN_data[key].append(untarg_dist_from_mean)
    untarg_dist_from_mean_table.append(untarg_dist_from_mean)

untarg_RMS = np.sqrt(np.mean(np.square(untarg_dist_from_mean_table), axis = 0))

np.savetxt('Untargeted.flm', np.transpose([untarg_integer_wave, untarg_composite_flux]), fmt="%d %26.18e")

#--------------------------- Plot both types

"""
#This code plots the targeted/untargeted composite spectra side by side with their residuals
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.set_title("Targeted SNe")
ax1.plot(targ_integer_wave, targ_composite_flux, color="#9acd32", label = "Composite")
ax1.plot(targ_integer_wave, (targ_composite_flux + targ_RMS), color="#ff6600", label = "+ RMS")
ax1.plot(targ_integer_wave, (targ_composite_flux - targ_RMS), color="c", label = "- RMS")
ax1.set_xlim(min(targ_integer_wave), max(targ_integer_wave))
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

ax2 = fig.add_subplot(223)
ax2.plot(targ_integer_wave, (targ_RMS/targ_composite_flux)*100, color="m", label = "Residual RMS")
ax2.set_xlim(min(targ_integer_wave), max(targ_integer_wave))
ax2.set_ylabel("Residual RMS")
ax2.legend(loc = "best")

ax3 = fig.add_subplot(222)
ax3.set_title("Untargeted SNe")
ax3.plot(untarg_integer_wave, untarg_composite_flux, color="#9acd32", label = "Composite")
ax3.plot(untarg_integer_wave, (untarg_composite_flux + untarg_RMS), color="#ff6600", label = "+ RMS")
ax3.plot(untarg_integer_wave, (untarg_composite_flux - untarg_RMS), color="c", label = "- RMS")
ax3.set_xlim(min(untarg_integer_wave), max(untarg_integer_wave))
ax3.set_xlabel("Rest Wavelength " + "(" + u"\u212B" + ")") #need to replace unicode?
ax3.set_ylabel("Relative Flux")
ax3.legend(loc = "best")

ax4 = fig.add_subplot(224)
ax4.plot(untarg_integer_wave, (untarg_RMS/untarg_composite_flux)*100, color="m", label = "Residual RMS")
ax4.set_xlim(min(untarg_integer_wave), max(untarg_integer_wave))
ax4.set_ylabel("Residual RMS")
ax4.legend(loc = "best")

fig.set_tight_layout(True)
plt.savefig("Targeted_vs_Untargeted_SidebySide.pdf", format="PDF")
plt.show()
"""

#This code plots the targeted/untargeted composite spectra together
#It needs to account for phase (try binning by a few days when grabbing from the database)
fig = plt.figure()
ax1 = fig.add_subplot(111)
xmin = min(min(targ_integer_wave), min(untarg_integer_wave))
xmax = max(max(targ_integer_wave), max(untarg_integer_wave))
ax1.plot(targ_integer_wave, targ_composite_flux, color="#ff6600", label = "Targeted")
ax1.plot(untarg_integer_wave, untarg_composite_flux, color="c", label = "Untargeted")
ax1.set_xlim(xmin, xmax)
ax1.set_xlabel("Rest Wavelength " + "(" + u"\u212B" + ")") #need to replace unicode?
ax1.set_ylabel("Relative Flux")
ax1.legend(loc = "best")
plt.savefig("Targeted_vs_Untargeted.pdf", format="PDF")
plt.show()