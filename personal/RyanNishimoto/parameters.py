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


#already dereddened, corrected for telluric absorbion
#in rest frame wavelength and flux density
#csp_files = glob.glob ('../../data/csp/*.flm')
#csp_sn = []

def query_bsnip():
	#Creates extinction.dat from bsnip data
	#list of files
	bsnip = glob.glob ('../../data/bsnip/*.flm')
	bsnip_sn = []
	for i in range (len(bsnip)):
		bsnip[i] = bsnip[i].split('-')[0]
		bsnip[i] = bsnip[i].split('/')
		bsnip[i] = bsnip[i][len(bsnip[i])-1]
		if i > 0:#removes redundant sn names
			if bsnip[i-1] != bsnip[i]:
				bsnip_sn.append(bsnip[i])
		b = []
		v = []
	for i in range(len(bsnip_sn)):
		print "looking at",bsnip_sn[i]
		#certain cases for specfic sn names
		if bsnip_sn[i] == 'sn2007s1':
			print "WARNING"
			del bsnip_sn[i:]
			break
		ext = IrsaDust.get_extinction_table(bsnip_sn[i])
		print ext
		print ext[1][0],ext[1][3]
		print ext[2][0],ext[2][3]
		b.append(ext[1][3])
		v.append(ext[2][3])
	#makes table in format 'sn','B','V'
	param = Table([bsnip_sn,b,v])
	ascii.write(param,'bsnip_extinc.dat')

def query_cfa():
	#gets cfa folder pathnames
	cfa = glob.glob ('../../data/cfa/*')
	#Truncates cfa pathname into needed format and removes cfa.dat files
	del cfa[0]
	del cfa[0]
	for i in range(len(cfa)):
		cfa[i]= cfa[i][15:]
		print cfa[i]
	b = []
	v = []
	for i in range(len(cfa)):
		print "looking at",cfa[i]
		ext = IrsaDust.get_extinction_table(cfa[i])
		print ext
		print ext[1][0],ext[1][3]
		print ext[2][0],ext[2][3]
		b.append(ext[1][3])
		v.append(ext[2][3])
	#makes table in format 'sn','B','V'
	param = Table([cfa,b,v])
	ascii.write(param,'extinctioncfa.dat')


def bsnip_edits():
	bsnip = glob.glob ('../../data/bsnip/*.flm')
	bsnip_sn = []
	b = []
	v = []
	for i in range (len(bsnip)):
			bsnip[i] = bsnip[i].split('-')[0]
			bsnip[i] = bsnip[i].split('/')
			bsnip[i] = bsnip[i][len(bsnip[i])-1]
			if i > 0:#removes redundant sn names
				if bsnip[i-1] != bsnip[i]:
					bsnip_sn.append(bsnip[i])
	cut_index = 0
	flag = False
	for i in range(len(bsnip_sn)):
			if (not flag):
				print "passing",bsnip_sn[i]
			#certain cases for specfic sn names
			if bsnip_sn[i] == 'sn2007s1':#want to start at first changed sn
				print "start with",len(bsnip_sn),"sn"
				print "START at index:",i
				cut_index = i
				flag = True #start here 
			if(flag):
				print "looking at",bsnip_sn[i]
				if bsnip_sn[i] == "sn2007s1":
					bsnip_sn[i] = "SNF20071021-000"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s1":
					bsnip_sn[i] = "SNF20080514-002"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008r3":
					bsnip_sn[i] = "sn1989a"###doesn't work
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s3":
					bsnip_sn[i] = "SNF20080825-006"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s4":
					bsnip_sn[i] = "sn1989a"###can't find correct name
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s5":
					bsnip_sn[i] = "SNF20080909-030"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s8":
					bsnip_sn[i] = "sn1989a"###look for correct names
					print "changed to",bsnip_sn[i]
				
				ext = IrsaDust.get_extinction_table(bsnip_sn[i])
				print ext
				print ext[1][0],ext[1][3]
				print ext[2][0],ext[2][3]
				b.append(ext[1][3])
				v.append(ext[2][3])


	print "stopped successfullly?...cutting table"
	del bsnip_sn[:cut_index]
	param = Table([bsnip_sn,b,v])
	print param
	ascii.write(param,'bsnip_extinc_added.dat')
def test():
	bsnip = glob.glob ('../../data/bsnip/*.flm')
	bsnip_sn = []
	b = []
	v = []
	for i in range (len(bsnip)):
			bsnip[i] = bsnip[i].split('-')[0]
			bsnip[i] = bsnip[i].split('/')
			bsnip[i] = bsnip[i][len(bsnip[i])-1]
			if i > 0:#removes redundant sn names
				if bsnip[i-1] != bsnip[i]:
					bsnip_sn.append(bsnip[i])
	cut_index = 0
	flag = False
	for i in range(len(bsnip_sn)):
			if (not flag):
				print "passing",bsnip_sn[i]
			#certain cases for specfic sn names
			if bsnip_sn[i] == 'sn2008s4':
				print "start with",len(bsnip_sn),"sn"
				print "START at index:",i
				cut_index = i
				flag = True #start here 
			if(flag):
				print "looking at",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008r3":
					bsnip_sn[i] = "ROTSE3 J125642.7+273041"###doesn't work
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s3":
					bsnip_sn[i] = "SNF20080825-006"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s4":
					bsnip_sn[i] = "SNF20080825-010"###can't find correct name
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s5":
					bsnip_sn[i] = "SNF20080909-030"
					print "changed to",bsnip_sn[i]
				if bsnip_sn[i] == "sn2008s8":
					bsnip_sn[i] = "SNF20080920-000"###look for correct names
					print "changed to",bsnip_sn[i]
				
				ext = IrsaDust.get_extinction_table(bsnip_sn[i])
				print ext
				print ext[1][0],ext[1][3]
				print ext[2][0],ext[2][3]
				b.append(ext[1][3])
				v.append(ext[2][3])


	print "stopped successfullly?...cutting table"
	del bsnip_sn[:cut_index]
	param = Table([bsnip_sn,b,v])
	print param
	#ascii.write(param,'bsnip_extinc_added.dat')
test()
