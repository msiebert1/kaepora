#!/usr/bin/env python -i
# Week 3 task
# Generate a composite spectrum

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import interpolate
from matplotlib.ticker import AutoMinorLocator
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

#Select the first 25 as the samples 
#snpaths = snpaths[0:24] 
sns = snnames[0:20]

#Find the maximum spectra for each SNe
#--Reading datadir/cfasnIa_param.dat, Column 3: MJD at B-band maximum light
#--Then converst to siderial time to match with *.flm filename 
#--Find the spectrum at maximum light.

#Sub-routine to load parameter file
def loadPara(parafile):
    f = open(parafile)
    lines = f.readlines()
    f.close()

    name = []
    tmax = []

    for line in lines:
        p = line.split()
        if p[0][0] != '#':
            name.append('sn'+p[0])
            tmax.append(float(p[2]))
        
    list = [name,tmax]
    return list


parafile = datadir+"cfasnIa_param.dat"
paras = loadPara(parafile)
names = paras[0]
tmaxs = paras[1]

# for sn in sns:
#     index = names.index(sn)
#     tmax = tmaxs[index]

#     list = glob.glob(sn+"/*.flm")

