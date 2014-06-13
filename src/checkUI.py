import glob
import Plotting
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import numpy as np
import os

from scipy import interpolate as intp
from prep import *
from datafidelity import *
import time


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
    
    badfile.append(os.path.split(filename)[1])
    badcomment.append(comment)
    badindex.append(index)
    

"""
this part is useless

def trimfile(filename):
	if source == 'cfa':
		return os.path.split(filename)[1]
	else:
		return os.path.split(filename)[1]  """   

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
if source == 'bsnip' :
    path = "../data/spectra/"+source+"/*.flm"
    sne = ReadExtin('extinctionbsnip.dat')         
if source == 'uv':
    path = "../data/spectra/"+source+"/*.flm"
    sne = ReadExtin('extinctionuv.dat')
if source == 'cfa':
    path = "../data/spectra/"+source+"/*/*.flm"
    sne = ReadExtin('extinction.dat')
    z = ReadParam()
if source == 'csp' :
    path = "../data/spectra/"+source+"/*.dat"
    sne = ReadExtin('extinctioncsp.dat')
if source == 'other' :
    path = "../data/spectra/"+source+"/*.dat"
    sne = ReadExtin('extinctionother.dat')

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
                print j , os.path.split(filenames[j])[1]
#				if source == 'cfa':
#					print j,filenames[j].split('/')[5]
#				else:
#					print j,filenames[j].split('/')[4]
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
orig_wave = []
orig_flux = []
orig_error = []
badfile = []
badcomment = []
badindex = []
interp = []
interp_wave = []
interp_flux = []
interp_error = []
invar = []
#ending index, remains -1 if gone through whole list
end = []

files = []
offset = start
i = start
#######################
## CLASS FOR BUTTONS ##
#######################
class Index:
	ind = i	
	def next(self, event):
		if self.ind == len(filenames)-1:
			self.ind = 0
		else:
			self.ind += 1
		process(self.ind,filenames)
		end.append(self.ind)

	def prev(self, event):
		if self.ind == 0:
			self.ind = len(filenames)-1
		else:
			self.ind -= 1
		process(self.ind,filenames)
		return self.ind
	def comm(self, event):
		print self.ind, filenames[self.ind]
		comment = str(raw_input("Comments on this spectra\n:"))
		tolists(filenames[self.ind],self.ind,comment)
		print "FILE",self.ind,':',filenames[self.ind]
		print "COMMENT:",comment


#plot first file we want to look at

##################
# Main Functions #
##################
def process(index, f):
	i = index
	try:

		#gets wave and flux from current file
            data = (np.loadtxt(f[i]))
#            print f[i]
            filename = os.path.split(f[i])[1]
		#separate wave/flux/error for easier manipulation
            print filename
            orig_wave = data[:,0]
            orig_flux = data[:,1]
            try:
			# check if data has error array
			orig_error = data[:,2] 
            except IndexError:
			# if not, set default
			orig_error = np.array([0])

            # deredenning
#            new_flux = dered(sne, filename, orig_wave, orig_flux)
            # deredshifting
#            if source != 'csp' :
#                new_wave = orig_wave/(1.+z)
		#get invar, to use in interp, and separate wave/flux/errors
            invar = genivar(orig_wave,orig_flux,orig_error)
            err = invar**-0.5
            place = np.where((orig_wave > 5800.0 ) & (orig_wave < 6000.0 ))
            print err[place]  
                     
            interp = Interpo(orig_wave,orig_flux,invar)
            interp_wave = interp[0]
            interp_flux = interp[1]
            interp_ivar = interp[2]
            interp_err = interp_ivar**-0.5
            
            place = np.where((interp_wave > 5800.0 ) & (interp_wave < 6000.0 ))
            print interp_err[place]
            
            # Just added : averaging weighted spectra
#            index = np.where(np.logical_not(np.isnan(interp_flux)))[0]
#            avgflux  = np.average(interp_flux[index], weights = interp_ivar[index], axis = 0)
#            print avgflux
#            weight_flux1 = orig_flux * invar 
#            weight_flux2 = interp_flux*interp_ivar 
            
            # Just Added: Clipping effect on flux
            index = np.where( interp_ivar == 0)
#            print index
#            interp_flux[index] = float('NaN')
            # the following replace the bad points to smoothed spectra
            var = np.ones(len(index), float)
            interp_flux[index] = gsmooth(interp_wave[index], interp_flux[index], var, 0.005)
            
            
		##plotting orig, interp, error
            xmin = orig_wave[0]-50
            xmax = orig_wave[len(orig_wave)-1]+50
		
            print "FILE",i,":",filename
		
            main.cla()
            main.set_xlim(xmin,xmax)
            main.plot(orig_wave,orig_flux,'b',label = 'Original')
            main.plot(interp_wave, interp_flux ,'r',label = 'Clipped')            
            main.set_xlabel('Rest Wavelength')
            main.set_ylabel('Flux')
            main.legend()
            main.set_title(filename)
            sub.cla()
            sub.set_xlim(xmin,xmax)
            sub.plot(orig_wave,err,label = 'Original')
            sub.plot(interp_wave,interp_err,label = 'Clipped')
            sub.legend()		
            sub.set_xlabel('Rest Wavelength')
            sub.set_ylabel('Error')
		
            plt.draw()
		#print "spectra is plotted"
		
		#NEW COMMENTING CODE
	except ValueError:
		print "found bad file! at index:",i
		tolists(filenames[i],i,'badfile')
		#can't read file ->messes up indexing and this corrects for this




###############
## MAIN CODE ##
###############
main = plt.subplot(2,1,1)  
sub = 	plt.subplot(2,1,2)
process(start,filenames)
callback = Index()
axprev = plt.axes([0.7, 0.025, 0.1, 0.05])
axnext = plt.axes([0.81, 0.025, 0.1, 0.05])
axcomm = plt.axes([.55,.025,.1,.05])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)
bprev = Button(axprev, 'Prev')
bprev.on_clicked(callback.prev)
bcomm = Button(axcomm, 'Comment')
bcomm.on_clicked(callback.comm)

plt.show()


end = end[len(end)-1]
print "ending at index:",end
##############
## QUITTING ##
##############
badlist =  Table([badindex,badfile,badcomment])
badlist_filename = "checked"+"_"+source+"_"+str(start)+"_"+str(end)+time.strftime("(%H:%M:%S)")
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

