#!/usr/bin/env python -i
# Week 3 task
# Generate a composite spectrum

import matplotlib.pyplot as plt
import numpy as np
import scipy
import glob
import fnmatch
import os
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec

datadir = '../../../data/cfa/'

#Load in all SN data directories in teh snpaths
#snpaths = []
snnames = []

for dirs,subdirs,files in os.walk(datadir):    
    for subdirs in fnmatch.filter(subdirs, 'sn*'):
#         snpaths.append(os.path.join(dirs,subdirs))
        snnames.append(subdirs)

snnames = np.array(snnames)
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

nbnames,vnebs,verrs = loadPara("../nebular.dat")
# Convert UT to MJD
#def ut2jd()
from astropy.time import Time

#Select the first 25 as the samples 
#sns = snnames[0:24]

# nsample = 50

# rand = np.random.randint(0,len(snnames),len(snnames))
# sns = snnames[rand]
sns = nbnames
nsample = 24

# from string import Template
import string

snsample = []
files = []
zsample = []
vnebsample = []
count = 0
for sn,vneb in zip(nbnames,vnebs):
    if count == nsample:
        break

    try:
        index = names.index(sn)
    except ValueError:
        continue

    index = names.index(sn)
    name = names[index]
    zhel = zhels[index]
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

#     fn = np.nanargmin(mdays)
    fn = np.nanargmax(mdays)
    
#     if days[fn] > 14:
    if days[fn] < 150:
        continue

    file = list[fn]

    snsample.append(name)
    zsample.append(zhel)
    files.append(file)
    vnebsample.append(vneb)
    
    count += 1


# snsample,vnebsample,zsample,files

def vbin(vmin,vmax):
    nfile=[]
    nz = []
    nv = []
    nname = []
    for file,z,vneb,name in zip(files, zsample,vnebsample,snsample):
        if vneb >= vmin and vneb < vmax:
           nfile.append(file)
           nz.append(z)
           nv.append(vneb)
           nname.append(name)
    return nfile,nz,nv,nname

import spectrum 
from spectrum import *

##Set plotting...
pltdir = '../plots/'


############################################

def mtspectra(vmin,vmax):
    files,zhels,nv,nname= vbin(vmin,vmax)
    spectrum(files,zhels)

# files,zhels,nv = vbin(0,5000)
# spectrum(files,zhels)

    plt.figtext(0.15,0.35,'Samples:',fontsize=10)

    plt.figtext(0.15,0.32,nname,fontsize=8)
    plt.figtext(0.15,0.29,nv,fontsize=8)

    plt.savefig(pltdir+'averspectra'+str(vmin)+'-'+str(vmax)+'.png')
    plt.show()
    return


