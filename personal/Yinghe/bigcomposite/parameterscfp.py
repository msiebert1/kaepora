"""
-Reads in file pathnames for spectra based on data origin
	-bsnip 	ex. -> sn1989a-19890427-o1i1.flm	already has adjusted flux?
	-cfa	ex. -> sn1993ac-19931016.49-mmt.flm	has weight in 3rd column
	-csp	ex. -> SN04dt_040914_b01_DUP_WF.dat	Redshift already accounted for (as well as earth's atmosphere?)

-Gets SN name from pathnames in form SNyearxx
	-different name truncation code depending on data origin

-Uses Astroquery to get parameter information
	-put each parameter into own lists
	-Writes lists into ascii table via astropy so paramenter 
	 data can be used in compprep.py


Table should look something like this:
SN_name	Host_Galaxy	R(v)	Redshift	B_val	V_val	Carbon_pos/neg

##################################
CURRENTLY ONLY DOES CSP DATA FILES
CURRENTLY RETRIEVES ONLY B/V VALS
##################################
"""
import glob
#import astroquery
from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np


#list of files
#bsnip = glob.glob ('../../data/bsnip/*.flm')
#cfa = glob.glob ('../../data/cfa/*')
csp = glob.glob ('../../../data/csp/*.dat')

#print csp
#Truncates cfa pathname into needed format and removes cfa.dat files
#del csp[0]
#del csp[0]
cspmini= [] # a shorter list

for i in range(len(csp)):    
    csp[i] = 'sn20'+csp[i][20:24]
#    print csp[i]
        
for sn in csp: # delete duplicate items
    if sn not in cspmini:
        cspmini.append(sn)
    
for i in range(len(cspmini)): # For single-letter SNs
    if cspmini[i].endswith('_') :        
        cspmini[i] = cspmini[i][:-1]
        
print cspmini    
    
b = []
v = []
for i in range(len(cspmini)):
	print "looking at SN",cspmini[i]
	ext = IrsaDust.get_extinction_table(cspmini[i])
	print ext
	print ext[1][0],ext[1][3]
	print ext[2][0],ext[2][3]
	b.append(ext[1][3])
	v.append(ext[2][3])

param = Table([cspmini,b,v],names=('sn','B','V'))
ascii.write(param,'extinctioncsp.dat')
	
