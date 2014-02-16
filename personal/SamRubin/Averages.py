'''
Created on Feb 14, 2014

@author: Alien
'''
import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math

class supernova(object):
    """Attributes can be added"""

SN_Array = []
names = np.loadtxt("input data", usecols = (0)) """need to grab all the other data too"""
for row in names:
    SN = supernova()
    SN.input = row[0]
    SN_Array.append(SN)
print len(names), "supernovae found"

file_list = glob.glob("../../data/cfa/*/*.flm")

con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()
for SN in SN_Array:
    for row in cur.execute('SELECT Filename, SN, Redshift, Phase, MinWave, MaxWave, Dm15, M B, B mMinusV M FROM Supernovae'):
        SN.filename = row[0]
        SN.name = row[1]
        SN.redshift = row[2]
        SN.phase = row[3]
        SN.minwave = row[4]
        SN.maxwave = row[5]
        SN.dm15 = row[6]
        SN.mb = row[7]
        SN.bminusv = row[8]
        SN_Array.append(SN)
        


high_wave = 1000
low_wave = 3000
waverange = high_wave-low_wave
compare_spectrum = SN_Array[0]
avg_flux = []

for SN in SN_Array:
    for i in range(len(waverange)):
        average = np.sum(SN.flux[i]/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0)/np.sum(1/SN.error[i]**2 for SN in SN_Array if SN.error[i] != 0 and SN.error[i] != None)
        avg_flux.append(average)
        compare_spectrum.flux[i] = average
        
