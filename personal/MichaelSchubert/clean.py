from __future__ import division
import numpy as np
import os
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import time

#make msgpack aware of msgpack_numpy
mn.patch()

tstart = time.clock()

def find_key(name):
    sn = name.split('-')
    #case for the SNF2008******* files, note the upper case requirement
    if sn[0][:3] == 'snf':
        keylist = [sn[0], sn[1]]
        key = '-'.join(keylist).upper()
    else:
        key = sn[0][2:]
    return key

#julian day - mjd =~ phase
#key,     0         1           2       3        4        5       6     7        8        9       10        11       12       13
#SN,  zhel,  tmax(B),  +/-  ref.,  Dm15, +/-  ref.,  M_B   +/-,   B-V,   +/-,  Bm-Vm, +/-,   Phot. ref
f = open('../../data/cfa/cfasnIa_param.dat')
lines = f.readlines()
sndict = {}
print "Populating dictionary"
for line in lines:
    if not line.startswith('#'):
        data = line.split()
        sndict[data[0]] = data[1:]
f.close()

#get time spectra was taken
f = open('../../data/cfa/cfasnIa_mjdspec.dat')
lines = f.readlines()
datedict = {}
for line in lines:
    if not line.startswith('#'):
        data = line.split()
        datedict[data[0]] = data[1]
f.close()

#initialize database
print "Creating table"
con = sq3.connect('SNe.db')
#cur = con.cursor()

#make sure no prior table in db to avoid doubling/multiple copies of same data
con.execute("""DROP TABLE IF EXISTS Supernovae""")
con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename 
                    TEXT PRIMARY KEY, SN Text, Redshift REAL, Phase REAL,
                    MinWave REAL, MaxWave REAL, RedMin REAL, RedMax REAL,
                    Dm15 REAL, M_B REAL, B_mMinusV_m REAL, Spectra BLOB)""")

root = '../../data/cfa'
spectra = {}
bad_files = []
filenames = []

#traverse folder containing spectra and load
print "Adding information to table"
for path, subdirs, files in os.walk(root):
    for name in files:
        filenames.append(name)
        f = os.path.join(path, name)
        if f.endswith('.flm'):
            #ignore spectra that produce loaderrors
            try:
                data = np.loadtxt(f)
            except:
                bad_files.append(f)
                continue

            #find the key to access param list dictionary
            key = find_key(name)
            #corner case for the redshiftless/dataless sn2011 files, no data = array of -1's
            if 'sn2011' not in name:
                vals = sndict[key]
                mjd = datedict[name]
                if vals[1] == '99999.9':
                    phase = None
                else:
                    phase = float(mjd) - float(vals[1])
            else:
                vals = [None]*14
                phase = None

            redshift = vals[0]
            if vals[4] == '9.99':
                dm15 = None
            else:
                dm15 = vals[4]
            #some dumb floating point rounding error happens if just pass in m_b = -99.99
            if vals[7] == '-99.99':
                m_b = None
            else:
                m_b = vals[7]
            if vals[11] == '-9.99':
                bm_vm = None
            else:
                bm_vm = vals[11]

            waves = data[:, 0]
            min_wave = waves[0]
            max_wave = waves[len(waves) - 1]
            if redshift is not None:
                red_min = min_wave/(1+float(redshift))
                red_max = max_wave/(1+float(redshift))
            else:
                red_min = None
                red_max = None

            spec = msg.packb(data)
            con.execute("""INSERT INTO Supernovae(Filename, SN, Redshift,
                                Phase, MinWave, MaxWave, RedMin, RedMax, Dm15, 
                                M_B, B_mMinusV_m, Spectra) 
                                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (name, key, 
                                redshift, phase, min_wave, max_wave, red_min, 
                                red_max, dm15, m_b, bm_vm, buffer(spec)))
        else:
            continue
con.commit()
tend = time.clock()

print tend - tstart


