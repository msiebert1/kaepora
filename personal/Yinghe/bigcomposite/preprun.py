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
spectra_files = glob.glob ('../../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []
#holds file name
file_name = []

junk_data = []

#number of spectra to modify
num = 20

#get data, pathnames
for i in range(num):
	try:
         spectra_data.append(np.loadtxt(spectra_files[i]))         
         file_path.append(spectra_files[i][14:-4])
         file_name.append(spectra_files[i][27:-4])    
#         print file_path
             
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)


#############################################################################################################
######################### Now processing ! ##################################################################

from prep import *

data = compprep(spectra_data,file_name,num)

print data

#######################################################################################################################
################### The rest is just for output testing################################################################
################### No need to implement this part to the database code ###############################################


# output data into a file 

#    output = 'testdata/modified-%s.dat'%(file_name[i])
#    ascii.write(newdata, output)
   
    # plot spectra 
#    x = newdata['col1']
#    y = newdata['col2']
#   z = newdata['col3']
#    print z
#    plt.subplot(1,2,1)
#    plt.plot(x,y)
#    plt.xlim(3000,7000)
#    plt.subplot(1,2,2)
#    plt.plot(x,z)
#    plt.xlim(3000,7000)

#plt.show()
#plt.savefig('test_host.png')


################## Finishline of output testing #######################################################################
#######################################################################################################################


# Output of noise

    #print navglist
#    ntable = Table([file_name,navglist],names=('spectra','noise'))   
#    ascii.write(ntable,'noise.dat')