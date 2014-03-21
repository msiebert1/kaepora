from __future__ import division
import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import time


mn.patch()

def read_cfa_or_bsnip(fname):
    """
    Returns a numpy array with spectra from a cfa or bsnip source
    """
    spectra = np.loadtxt(fname)
    return spectra

def read_csp(fname):
    """
    Returns a spectra from a csp source as well as the associated information.
    Info is a list with fields [SN, File, Redshift, Date Max, Date Obs, Epoch]
    """
    spectra = np.loadtxt(fname)
    with open(fname) as f:
        info = [f.next().rstrip().split()[1] for x in range(6)]

    return spectra, info

def read_cfa_info(data_file, dates_file):
    """
    Reads Cfa SN information from separate files.
    Output dict format:
    #key,     0         1           2       3        4        5       6     7        8        9       10        11       12       13              14
    #SN,  zhel,  tmax(B),  +/-  ref.,  Dm15, +/-  ref.,  M_B   +/-,   B-V,   +/-,  Bm-Vm, +/-,   Phot. ref   Obs Date
    """
    with open(dates_file) as dates:
        lines = dates.readlines()
        date_dict = {}
        for line in lines:
            if not line.startswith('#'):
                cleandate = line.split()

                if cleandate:
                    date_dict[cleandate[0]] = cleandate[1]

    with open(data_file) as data:
        lines = data.readlines()
        sndict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                sndict[sndata[0]] = sndata[1:]

    return sndict, date_dict

def read_bsnip_data(data_file):
    #  Supernova Name (1)
    #  SNID (Sub)Type (2)
    #  Host Galaxy
    #  Host Morphology (3)
    #  Heliocentric Redshift, cz (km/s) (4)
    #  Milky Way Reddening, E(B-V) (mag) (5)
    #  UT date of discovery
    with open(data_file) as data:
        data = pd.read_fwf('obj_info_table.txt', names=('SN name', 'Type',
                'Host Galaxy', 'Host Morphology', 'Redshift', 'Reddening',
                'Discovery Date'), colspecs=((0,10),(10,17),(17,51),(51,58),
                (58,63),(63,69),(69,97)))
        dmat = data.as_matrix()
        bsnip_dict = {}
        for line in dmat:
            key = line[0].split()[1].lower()
            rs = line[4]
            vals = [rs]
            bsnip_dict[key] = vals
        return bsnip_dict

def find_SN(fname, source=None, csplist=None):
    """
    Returns SN name, either from cspdict if source is a csp spectra
    or from slicing the file name if source is Cfa or bsnip
    """
    if source == 'csp':
        snname = csplist[0]
        return snname[2:]
    else:
        snname = fname.replace('_', '-').split('-')
        if snname[0][:3] == 'snf':
            namelist = [snname[0], snname[1]]
            snname = '-'.join(namelist).upper()
        else:
            snname = snname[0][2:]

        return snname


#change this depending on where script is
sndict, date_dict= read_cfa_info('../data/cfa/cfasnIa_param.dat',
                                                       '../data/cfa/cfasnIa_mjdspec.dat')

ts = time.clock()
con = sq3.connect('SNe.db')

#make sure no prior table in db to avoid doubling/multiple copies of same data
con.execute("""DROP TABLE IF EXISTS Supernovae""")
con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename
                    TEXT PRIMARY KEY, SN Text, Source Text, Redshift REAL, Phase REAL,
                    MinWave REAL, MaxWave REAL, Dm15 REAL, M_B REAL,
                    B_mMinusV_m REAL, Targeted INTEGER, Snr REAL,
                    Interpolated_Spectra BLOB)""")

#change this depending on where script is
root = '../data'
bad_files = []
bad_interp = []
bsnip_vals = read_bsnip_data('obj_info_table.txt')
print "Adding information to table"
for path, subdirs, files in os.walk(root):
    for name in files:
        f = os.path.join(path, name)
        if f.endswith('.flm') or f.endswith('.dat'):
            if 'cfasnIa' in f:
                continue
            try:
                if 'csp' in path:
                    spectra, info = read_csp(f)
                    sn_name = find_SN(name, 'csp', info)
                else:
                    spectra = read_cfa_or_bsnip(f)
                    sn_name = find_SN(name)
            except:
                bad_files.append(f)
                continue

            #finds cfa data for particular sn if applicable
            if 'cfa' in f:
                if 'sn2011' not in name:
                    sn_cfa = sndict[sn_name]
                else:
                    snd = None
                    sn_cfa = [None] * 14

            #csp source
            if 'csp' in f:
                source = 'csp'
                redshift = info[2]
                phase = float(info[4]) - float(info[3])
                Dm15 = None
                m_b = None
                bm_vm = None

            #cfa spectra
            elif 'cfa' in f:
                source = 'cfa'
                redshift = sn_cfa[0]
                if  sn_cfa[1] == '99999.9':
                    phase = None
                else:
                    #try/except catches and fixes sn2011 errors
                    try:
                        phase = float(date_dict[name]) - float(sn_cfa[1])
                    except:
                        phase = None
                if sn_cfa[4] == '9.99':
                    Dm15 = None
                else:
                    Dm15 = sn_cfa[4]

                if sn_cfa[7] == '-99.99':
                    m_b = None
                else:
                    m_b = sn_cfa[7]

                if sn_cfa[11] == '-9.99':
                    bm_vm = None
                else:
                    bm_vm = sn_cfa[11]

            #bsnip spectra
            else:
                c = 299792.458
                source = 'bsnip'
                data = bsnip_vals[sn_name.lower()]
                redshift = data[0]/c
                phase = None
                Dm15 = None
                m_b = None
                bm_vm = None

            waves = spectra[:, 0]
            min_wave = waves[0]
            max_wave = waves[len(spectra) - 1]
            spec = msg.packb(spectra)

            try:
                interp_spec, sig_noise = prep.compprep(spectra, sn_name, redshift, source)
            except:
                print "Interp failed"
                bad_interp.append(name)
                bad_files.append(name)
                interp_spec, sig_noise = None, None

            interped  = msg.packb(interp_spec)
            con.execute("""INSERT INTO Supernovae(Filename, SN, Source, Redshift, Phase,
                                MinWave, MaxWave, Dm15, M_B, B_mMinusV_m, Snr, Interpolated_Spectra)
                                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (name, sn_name, source, redshift,
                                phase, min_wave, max_wave, Dm15, m_b, bm_vm, sig_noise, buffer(interped)))
con.commit()
te = time.clock()
print 'bad files', bad_files, len(bad_files)
print 'bad interps', bad_interp, len(bad_interp)
print te - ts