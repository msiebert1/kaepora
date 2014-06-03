import glob
import Plotting
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import time
import sqlite3 as sq3
from scipy import interpolate as intp
import math
import msgpack as msg
import msgpack_numpy as mn
from prep import *
from datafidelity import *

"""
*Meant to be used from src folder

This code is designed to check each spectra file, one by one, given a particular source,
by plotting the original and the interpolated spectra together, along with the error.  

comments can be given to note telluric absorption, *Clipping issues, and other features that
should be weighted differently, as well as any other issues with the spectra/data
"""

"""
ROUGH OUTLINE OF CODE
-select SOURCE (bsnip,cfa,csp,other,uv)
-pull all FILENAMES from selected source
-select START index [0 to (len(source)-1))]
	-can check the index and corresponding file with input '-1'

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
(important) variables: (work in progress)
==========
SOURCE		which source data is being looked at(bsnip,cfa,etc.)
FILENAMES	path to the spectra files
FILES		truncated path list to match the START (Mostly useless(harmless))***

DATA		data from the spectra files
ORIG_WAVE	separate original wavelength from DATA
ORIG_FLUX	separate original flux from DATA
INTERP		interpolated data
INTERP_WAVE	separate interpolated wavelength from INTERP
INTERP_FLUX	separate interpolated flux from INTERP
INVAR		inverse variance from Datafidelity code

BADINDEX	holds the index of the bad file
BADFILE		holds filename that has problems
BADCOMMENT	holds the comment associated with the bad file
BADLIST		combines BADFILE and BADCOMMENT into one table

START		index started at
END		index ended at

"""
######################
## Helper Functions ##
######################
def tolists(filename,index,comment):
	if source == 'cfa':
		badfile.append(filenames[i].split('/')[5])
	else:
		badfile.append(filenames[i].split('/')[4])

	badcomment.append(comment)
	badindex.append(i)


####################
## File Selecting ##
####################

#Select which data you want to check, will loop until a correct data source is selected
while (True):
	source = str(raw_input("Please select which data to check(type one of the following):\nbsnip\ncfa\ncsp\nother\nuv\n:"))
	if(source != 'bsnip' and source !='cfa' and source !='csp' and source !='other' and source !='uv'):
		print "\ninvalid data selection\n"
	else:
		print "checking",source
		break

#accounts for different file organizations for sources
if source == 'bsnip' or source == 'uv':
	path = "../data/spectra/"+source+"/*.flm"
if source == 'cfa':
	path = "../data/spectra/"+source+"/*/*.flm"
if source == 'csp' or source == 'other':
	path = "../data/spectra/"+source+"/*.dat"


#read in all filenames to extract data from 
filenames = glob.glob(path)

###checks that all the files were correctly read in
##for files in filenames:
##	print files

#select an index value to start at so you don't have to check a set of data
#from the beginning every time you run the program
while (True):
	print "valid index range: 0 to",(len(filenames)-1)
	try:
		start = int(raw_input("Please select what index you would like to start at?(-1 for list)\n:"))
		if start == -1:
			step = int(raw_input("Please enter step size(how many spectra per [enter])\n:"))
			print "PRESS [ENTER] TO SHOW NEXT",step,"(q to quit)"
			for j in range(len(filenames)):
				if source == 'cfa':
					print j,filenames[j].split('/')[5]
				else:
					print j,filenames[j].split('/')[4]
				if (j % step == 0):
					flag = raw_input();
					if flag == 'q':
						break
		if (start < 0 or start > (len(filenames)-1)):
			print "INDEX OUT OF RANGE\n"
		else:
			print "starting at index:",start		
			break
	except ValueError:
		print "NOT AN INTEGER\n"


data = []
wave = []
flux = []
badfile = []
badcomment = []
badindex = []
interp = []
invar = []
#ending index, remains -1 if gone through whole list
end = -1

"""
[ ] - need to implement
[x] - implemented
[!] - implemented but needs work
MAIN LOOP CHECKLIST:

	For each spectra file:

	[x]read data at index (from START)
	[x]except value error for bad file
	[x]make sure FILE lines up with DATA
	[x]create separate interpolated data INTERP
	[x]get inverse varience values INVAR
#	[!]plot ORIG and INTERP data (problems using Plotting.py)
#	[!]plot error (issue with error column of interpolated data,
	   right now reversing the inverse variance to error)
	[x]Observe and Comment(for BADFILE and BADCOMMENT):
		[x][enter](blank input)
		 	plot is fine(don't add to list)
		[x]COMMENT(telluric,clipping,etc...)
		 	add FILES[i] to BADFILE and 
			keyboard input to BADCOMMENT
			and index to BADINDEX
		[x]'q' to end, saves ending index to END
			[x]combine BADFILE/BADCOMMENT/BADINDEX into 3-column table
#			[!]save table (no safeguard for overwriting file of same
			   name, but issue is local; another issue is commiting file
			   of same name when both files have different comments)
			   FORMAT:'checked_SOURCE_START-END'
			   ex.:'checked_cfa_0-200'
	[x]Continues onto next file if not quitting

	


"""
files = []
###############
## MAIN LOOP ## see MAIN LOOP CHECKLIST above for breakdown
###############

#explanation for offset for my own sanity and
#fear of mislabeling data
"""
try appending next data, can't if i is not lined up.

offset = start so that it'll take into account starting
index and adjust the data being looked at correctly

(start = 1 instead of 0, won't run first time through,
want to look at 2nd filename, so it'll append filename[1]
to data, data length now 1.  access data[0] (i-start or 
1-1=0) <--should work fine

CASE:bad data read (ExceptValueError):
Given: start at 5, filenames[6] is bad, move on to 7
i >=5, read in filename 5, data length 1, access data
(5-5 = 0th index [good]). i=6 filename[6] is bad, move
on.  i = 7 now, offset = start+1 = 6.  read filenames[7]
append to data (size 2).  Want to access data[1], so
7-6 = 1[GOOD]

"""
offset = start
for i in range(len(filenames)):
    if (i >= start):
        try:
			###checks that correct data files are appended 
			###when taking into account the index offset
			##print i,filenames[i]

			#gets wave and flux from current file
            data.append((np.loadtxt(filenames[i])))
			#keeps track of files looking at (due to index offset)
            if source == 'cfa':
                files.append(filenames[i].split('/')[5])
            else:
			files.append(filenames[i].split('/')[4])
            print files
		
			#separate wave/flux/error for easier manipulation
            orig_wave = data[i-offset][:,0]
            orig_flux = data[i-offset][:,1]
            try:
				# check if data has error array
                    orig_error = data[i-offset][:,2] 
            except IndexError:
				# if not, set default
                    orig_error = np.array([0])

			##get invar, to use in interp, and separate wave/flux/errors
            invar = genivar(orig_wave,orig_flux,orig_error)
#                print invar                           
            interp = Interpo(orig_wave,orig_flux,invar)
            interp_wave = interp[0]
            interp_flux = interp[1]
            interp_ivar = interp[2]
#                print interp_wave,interp_flux,interp_ivar


		##plotting orig, interp, error
            """
Plotting.py code that isn't compatible with this code right now
Relative_Flux = [orig_wave, orig_flux, interp_wave, interp_flux]  # Want to plot a composite of multiple spectrum
Variance = invar
Residuals     = []
Spectra_Bin   = []
Age           = []
Delta         = []
Redshift      = [] 
Show_Data     = [Relative_Flux,Variance,Residuals,Spectra_Bin,Age,Delta,Redshift]
image_title   = source,i            # Name the image (with location)
title         = "checking",filenames[i]
		#
		## Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
		##                   0              1          2            3    4      5         6,                 7
		####################################^should be varience?...
                name_array = ["original", " ", "interpolated", " "]
                Names = name_array
                Plots = [0,1] # the plots you want to create
		#
		## The following line will plot the data
		#
                xmin         = orig_wave[0]-50 
                xmax         = orig_wave[len(orig_wave)-1]+50
                # Use the plotting function here
                #Plotting.main(Show_Data , Plots , image_title , title, Names,xmin,xmax)
"""
            xmin = orig_wave[0]-50 
            xmax = orig_wave[len(orig_wave)-1]+50
            print "FILE",i,":",filenames[i]
		
                #test plotting (when Plotting code is not working properly)
                #plt.figure(1)
            plt.subplot(2,1,1)                
            plt.plot(orig_wave,orig_flux,'b',label = 'Original')
            plt.plot(interp_wave, interp_flux,'r',label = 'Interpolated')
            plt.xlim(xmin,xmax)
            plt.xlabel('Rest Wavelength')
            plt.ylabel('Flux')
		
            plt.legend(loc="upper right")
            plt.title(files[i-offset])
		
            plt.subplot(2,1,2)
            plt.plot(orig_wave,invar**-0.5,label = 'Original')
            plt.plot(interp_wave,interp_ivar**-0.5,label = 'Clipped')
            plt.xlim(xmin,xmax)
            plt.xlabel('Rest Wavelength')
            plt.ylabel('Error')
            plt.legend()
            '''
		print "Checking Original data"
		print orig_wave
		print orig_flux
		print invar
		print "checking Interpolated data"
		print interp_wave
		print interp_flux
		print interp_ivar
		'''
            plt.show()
                #print "spectra is plotted"

            comment = str(raw_input("Please comment on this spectra\n([enter](blank) = no error, 'q' or 'quit' to stop)\n:"))
			#no error, don't record anything
            if comment == '':
			print "nothing to see here"
			print "move along, move along"
			#done checking, record ending index and stop loop
            elif comment == 'q' or comment == 'quit':
			end = i
			break
			#comment made, record it and the file to respective lists
            else:
			tolists(filenames[i],i,comment)
			print "COMMENT:",comment
				
				
        except ValueError:
		print "found bad file! at index:",i
		tolists(filenames[i],i,comment)
		#can't read file ->messes up indexing and this corrects for this
		offset += 1

#did not quit early, update end to last index
if end == -1:
	end = len(filenames)-1
###checks to make sure the right data is being looked at
##print "looking at",files[0],"with wave vals:\n",data[0][:,0]	#all wavelengths of data at index 0
##print "looking at",files[0],"with flux vals:\n",data[0][:,1]	#all flux values of data at index 0

###double-check to make sure right number of data are appended
##print len(filenames),"- starting index",start,"=",len(data)
##print badfile,badcomment

##############
## QUITTING ##
##############
badlist =  Table([badindex,badfile,badcomment])
badlist_filename = "checked"+"_"+source+"_"+str(start)+"-"+str(end)+time.strftime("(%H,%M,%S)")
print "REVIEW:"
print "BADLIST FILENAME:",badlist_filename
print "LIST:\n",badlist
while(True):
	save = str(raw_input("save?(y/n)\n:"))
	if save == 'y':
		print "saving..."
		ascii.write(badlist,badlist_filename)
		break
	elif save == 'n':
		safety = str(raw_input("are you sure you want to quit?(y/n)\n:"))
		if safety == 'y':
			print "quitting without saving..."
			break
		elif safety =='n':
			continue

	else:
		print "not valid input"

