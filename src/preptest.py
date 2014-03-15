# -*- coding: utf-8 -*-
"""
Created on Sun Mar 09 22:39:51 2014

@author: QuantumMonkey
"""
import os
import glob
from specutils import extinction as ex
import astroquery
from astroquery.ned import Ned
from astroquery.irsa_dust import IrsaDust
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from math import floor,ceil

#################################################################################################
########################### Read in Files #######################################################

#con = sq3.connect('SNe.db')
#cur = con.cursor()

#list of files
spectra_files = glob.glob ('../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []
#holds file name
file_name = []

junk_data = []

#number of spectra to modify
#num = len(spectra_files)
num = 20

#get data, pathnames
for i in range(num):
	try:
         spectra_data.append(np.loadtxt(spectra_files[i]))         
         file_path.append(spectra_files[i][14:-4])
         file_name.append(spectra_files[i][20:-4])    
#         print file_name
             
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)

#print num

############################# this is for testing typical files####################################
def testonesn(spectra_data,file_name) :
    for i in range(num) :
        if 'sn2008Z' in file_name[i] :
            return i

            
def ReadIn(file):
    f = open(file,'r')
    d = []
    for line in f.readlines():
        line=line.decode('utf-8-sig')
        line=line.encode('UTF-8')
    #    print line.split()
        d.append([float(value) for value in line.split()])    
    f.close()
    return d
        
def ImpArr(d) :
    wave = []
    flux = []  
    var = []         
    for i in range(len(d)) :
        wave.append(d[i][0])
        flux.append(d[i][1])
        var.append(d[i][2])
#    print 'wave',wave    
    new = np.array([wave,flux,var])    
    return new  

#############################################################################################################
######################### Now processing ! ##################################################################

from prep import *

navglist = [] 
      
for i in range(num) :  #go through selected spectra data
#    i = testonesn(spectra_data,file_name)
#    print i
    spectrum = spectra_data[i]	#declares new spectrum from list
    
#sn =  '../data/cfa/sn2001V/sn2001V-20010325.40-mmt.flm'   
#data = ReadIn(sn)
#spectrum = ImpArr(data)
#print spectrum 
#file_name=sn[27:-4]
#print file_name
    try :
        data = compprep(spectrum,file_name[i])
        navglist.append(data[1]) # the S/N ratio
#        print data
    except ValueError:
        print "ignore",file_name[i]
        navglist.append('NaN')
#    wave = data[0][0]
#    print wave
#    flux = data[0][1]
#    var = data[0][2]    
    


#######################################################################################################################
################### The rest is just for output testing################################################################
################### No need to implement this part to the database code ###############################################


# output data into a file 

#    output = 'testdata/modified-%s.dat'%(file_name[i])
#    ascii.write(newdata, output)
   
    # plot spectra 

#    plt.subplot(1,2,1)
#    plt.plot(wave,flux)
#    plt.xlim(3000,7000)
#    plt.subplot(1,2,2)
#    plt.plot(wave,var)
#    plt.xlim(3000,7000)

#plt.show()
#plt.savefig('test_host.png')


################## Finishline of output testing #######################################################################
#######################################################################################################################


# Output of noise

#print navglist
ntable = Table([file_name,navglist],names=('spectra','S/N ratio'))
#print ntable   
ascii.write(ntable,'snr.dat')

# Plotting of noise
listi = range(num)
#print listi
plt.plot(listi,navglist)
plt.show()


