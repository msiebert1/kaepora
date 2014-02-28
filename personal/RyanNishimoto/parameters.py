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
CURRENTLY ONLY DOES CFA DATA FILES
CURRENTLY RETRIEVES ONLY B/V VALS
##################################
"""
import glob
import astroquery
from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np


#list of files
bsnip = glob.glob ('../../data/bsnip/*.flm')
cfa = glob.glob ('../../data/cfa/*')
csp = glob.glob ('../../data/csp/*.flm')

#Truncates cfa pathname into needed format and removes cfa.dat files
del cfa[0]
del cfa[0]
for i in range(len(cfa)):
	cfa[i]= cfa[i][15:]
	print cfa[i]
b = []
v = []
for i in range(len(cfa)):
	print "looking at SN",cfa[i]
	ext = IrsaDust.get_extinction_table(cfa[i])
	print ext
	print ext[1][0],ext[1][3]
	print ext[2][0],ext[2][3]
	b.append(ext[1][3])
	v.append(ext[2][3])

param = Table([cfa,b,v],names=('sn','B','V'))
ascii.write(param,'extinction.dat')
	
