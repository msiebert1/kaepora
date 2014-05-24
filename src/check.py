import glob
import composite
import Plotting
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import prep

import sqlite3 as sq3
from scipy import interpolate as intp
import math
import msgpack as msg
import msgpack_numpy as mn

"""
*Mean't to be used from src folder

This code is designed to check each spectra file, one by one, given a particular source,
by plotting the original and the interpolated spectra together, along with the 
inverse varience.  

comments can be given to note telluric absorption, errors in reddening, spikes in fluxes,
clipping, or any other problems with either the original or interpolated spectra
"""

"""
ROUGH OUTLINE OF CODE
-select SOURCE (bsnip,cfa,csp,other,uv)
-pull all FILENAMES from selected source
-select START index [0 to (len(source)-1))]

MAIN LOOP:	
	-read data at index (from START)
	-create separate interpolated data INTERP
	-get inverse varience values INVAR
	-plot ORIG and INTERP data 
	 with INVAR
	-Observe and Comment:
		-[enter](blank input)
		 	plot is fine(don't add to list)
		-type comment(telluric, spikes, etc...)
		 	added to 2 column list-> |FILENAME|COMMENT|
		-'q' or 'quit' to end, saves ending index END
	-Continues onto next file

	When Quitting:
		save list to file in format: 
		('SOURCE_'+'START'+'_'+'END')
		ex. ('cfa_0_200.dat')

===========
guidelines:
===========
[ ] - need to implement
[x] - implemented

#comment on the functionality of the blocks of code
##commented out code that does checks on the functionality of code
==========
variables: (work in progress)
==========
SOURCE		which source data is being looked at(bsnip,cfa,etc.)
FILENAMES	path to the spectra files
FILES		truncated path list to match the data being looked at
DATA		data from the spectra files
LIST		holds n column list that has the filenames, comment on errors

START		index started at
END		index ended at

"""

#Select which data you want to check, will loop until a correct data source is selected
while (True):
	source = str(raw_input("Please select which data to check(type one of the following):\nbsnip\ncfa\ncsp\nother\nuv\n:"))
	if(source != 'bsnip' and source !='cfa' and source !='csp' and source !='other' and source !='uv'):
		print "\ninvalid data selection\n"
	else:
		print "checking",source
		break

#accounts for different file organizations for sources
if source == 'bsnip':
	path = "../data/spectra/"+source+"/*.flm"
if source == 'cfa':
	path = "../data/spectra/"+source+"/*/*.flm"
if source == 'csp':
	path = "../data/spectra/"+source+"/*.dat"
if source == 'other':
	path = "../data/spectra/"+source+"/*.dat"
if source == 'uv':
	path = "../data/spectra/"+source+"/*.flm"

#read in all filenames to extract data from 
filenames = glob.glob(path)

###checks that all the files were correctly read in
##for files in filenames:
##	print files

#select an index value to start at so you don't have to check a set of data
#from the beginning every time you run the program
while (True):
	print "valid index range: 0 to",(len(filenames)-1)
	start = int(raw_input("Please select what index you would like to start at:\n"))
	if (start < 0 and start > (len(filenames)-1)):
		print "\ninvalid index"
	else:
		print "starting at index:",start		
		break

"""
#trims filenames so that cuts out the path and leaves only the file
for i in range(len(filenames)):
	if source == 'cfa':
		filenames[i] = filenames[i].split('/')[5]
	else:
		filenames[i] = filenames[i].split('/')[4]

###checks that all filenames were correctly split
##for files in filenames:
##	print files
"""

data = []
wave = []
flux = []
badlist = []
badcomment = []
interp = []
invar = []
end = -1

"""
[ ] - need to implement
[x] - implemented

MAIN LOOP CHECKLIST:

	For each spectra file:

	[x]read data at index (from START)
	[x]except value error for bad file
	[x]make sure FILE lines up with DATA
	[]create separate interpolated data INTERP
	[]get inverse varience values INVAR
	[]plot ORIG and INTERP data 
	 with INVAR
	[]Observe and Comment(for BADLIST and BADCOMMENT):
		[][enter](blank input)
		 	plot is fine(don't add to list)
		[]type comment(telluric, spikes, clipping, etc...)
		 	add FILES[i] to BADLIST and 
			keyboard input to BADCOMMENT
		[]'q' or 'quit' to end, saves ending index to END
	[x]Continues onto next file if not quitting

	if quitting:
		[]combine BADLIST and BADCOMMENT into 2-column list
		[]save list to txt file in format: 
		('SOURCE_'+'START'+'_'+'END')
		ex. ('cfa_0_200.txt')

"""
files = []

#main loop, see MAIN LOOP CHECKLIST above for breakdown
for i in range(len(filenames)):
	if (i >= start):
		try:
			###checks that correct data files are appended 
			###when taking into account the index offset
			##print i,filenames[i]
			data.append((np.loadtxt(filenames[i])))
			files.append(filenames[i])

		except ValueError:
			print "found bad file! at index:",i
			badlist.append(filenames[i])
			badcomment.append("bad file")


###checks to make sure the right data is being looked at
##print "looking at",files[0],"with wave vals:\n",data[0][:,0]	#all wavelengths of data at index 0
##print "looking at",files[0],"with flux vals:\n",data[0][:,1]	#all flux values of data at index 0

###double-check to make sure right number of data are appended
##print len(filenames),"- starting index",start,"=",len(data)
##print badlist,badcomment





