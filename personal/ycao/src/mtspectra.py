#!/usr/bin/env python -i
# Week 3 task
# Generate a composite spectrum

import matplotlib.pyplot as plt
import numpy as np
import scipy
import glob
import fnmatch
import os

datadir = '../../../data/cfa/'

#Load in all SN data directories in teh snpaths
#snpaths = []
snnames = []

for dirs,subdirs,files in os.walk(datadir):    
    for subdirs in fnmatch.filter(subdirs, 'sn*'):
#         snpaths.append(os.path.join(dirs,subdirs))
        snnames.append(subdirs)

#Find the maximum spectra for each SNe
#--Reading datadir/cfasnIa_param.dat, Colum1, 2, 3: Name, redshift, MJD at B-band maximum light
#--Then converst UT time to MJD to match *.flm filename with tmax 
#--Find the spectrum at maximum light.

#Sub-routine to load parameter file
def loadPara(parafile):
    f = open(parafile)
    lines = f.readlines()
    f.close()

    names = [] # SN name
    zhels = [] # heliocentric redshift
    tmaxs = [] # MJD at B-band maximum light (99999.9=no estimate)


    for line in lines:
        p = line.split()
        if p[0][0] != '#':
            names.append('sn'+p[0])
            zhels.append(float(p[1]))
            tmaxs.append(float(p[2]))
    return names,zhels,tmaxs

names,zhels,tmaxs = loadPara(datadir+"cfasnIa_param.dat")

# Convert UT to MJD
#def ut2jd()
from astropy.time import Time

#Select the first 25 as the samples 
#sns = snnames[0:24]
sns = snnames[0:30]
print 'Our sample: ',sns

# from string import Template
import string

for sn in sns:
    index = names.index(sn)
    tmax = tmaxs[index]
    if tmax == 99999.9: 
        continue

    t = Time(tmax, format='mjd', scale='utc')
    date = t.iso[0:10]
    daymax = int(''.join([i for i in date if i.isdigit()]))

    list = glob.glob(datadir+sn+"/*.flm")
    
    dates=[]
    for file in list:
        hdind = str.find(file,'-') + 1
        yyyy = file[hdind:hdind+4]
        mm = file[hdind+4:hdind+6]
        dd = file[hdind+6:hdind+8]
        date = int(yyyy+mm+dd)
        dates.append(date)
        
    dates = np.array(dates)

    days = dates - daymax
    mdays = np.ma.array(days,mask = days < 0.)# Calculate the delt_day and mask negtive values

    fn = np.nanargmin(mdays)

    file = list[fn]
    
    if days[fn] > 8:
        continue

    print daymax, dates[fn], days[fn]

        
