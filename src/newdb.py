from __future__ import division
import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import re
import math
import time
import photometry as phot

mn.patch()
global c
c = 299792.458


def read_cfa_or_bsnip_or_uv(fname):
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
    # key,     0         1           2       3        4        5       6     7        8        9
    # SN,  zhel,  tmax(B),  +/-  ref.,  Dm15, +/-  ref.,  M_B   +/-,   B-V,
     10        11       12       13              14
    +/-,  Bm-Vm, +/-,   Phot. ref   Obs Date
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
    """
    Returns a dict keyed on 'snname-longdate' ie. 1989a-19890427
    """
    data1 = pd.read_fwf('../data/info_files/obj_info_table.txt', names=('SN name', 'Type',
                                        'Host Galaxy', 'Host Morphology', 'Redshift', 'Reddening',
                                        'Discovery Date'), colspecs=((0, 10), (10, 17), (17, 51),
                                        (51, 58), (58, 63), (63, 69), (69, 97)))

    data2 = pd.read_fwf('../data/info_files/spec_info_table.txt',
                        names=('SN name', 'Year', 'Month', 'Day',
                                       'MJD of Spectrum', 'Phase of Spectrum',),
                        colspecs=((0, 9), (10, 15), (15, 18), (18, 25),
                                         (25, 36), (36, 43)))
    dmat1 = data1.as_matrix()
    dmat2 = data2.as_matrix()
    dmat2 = dmat2.tolist()
    #format month as '04' instead of '4'
    #do the same for day (5.000 -> 05.000)
    for line in dmat2:
        if line[2] < 10:
            line[2] = ''.join(['0', str(line[2])])
        if line[3] < 10:
            line[3] = ''.join(['0', str(line[3])])
    for line in dmat2:
        fulldate = ''.join([str(line[1]), str(line[2]), str(line[3])[:2]])
        line.append(fulldate)

    bsnip_dict_rs = {}
    bsnip_dict = {}
    #split matrix lines to populate dicts
    for line in dmat1:
        key = line[0].split()[1].lower()
        rs = line[4]
        bsnip_dict_rs[key] = rs
    for line in dmat2:
        key1 = line[0].split()[1].lower()
        key2 = line[len(line)-1]
        full_key = key1+'-'+key2
        bsnip_dict[full_key] = [key1, bsnip_dict_rs[key1], line[len(line)-2]]
    return bsnip_dict


def find_SN(fname, source=None, csplist=None):
    """
    Returns SN name, either from cspdict if source is a csp spectra
    or from slicing the file name if source is Cfa or bsnip
    """
    if source == 'csp':
        snname = csplist[0]
        return snname[2:]
    elif source == 'other':
        snname = fname.replace('_', '-').split('-')
        if snname[0][:2] == 'sn':
            return snname[0][2:]
        else:
            return snname[0]
    else:
        snname = fname.replace('_', '-').split('-')
        if snname[0][:3] == 'snf':
            namelist = [snname[0], snname[1]]
            snname = '-'.join(namelist).upper()
        else:
            snname = snname[0][2:]

        return snname


def build_morph_dict():
    """
    Builds a dictionary of the form {sn_name: galaxy morphology}
    """
    with open('../data/info_files/host_info.dat') as f:
        txt = f.readlines()

    clean = [x.strip() for x in txt]
    outlist = [item for item in clean if not item.startswith('#')]
    morphlist = filter(None, outlist)

    morph_dict = {}
    for entry in morphlist:
        ents = entry.split()
        morph_dict[ents[0]] = ents[1]

    return morph_dict


def build_vel_dict():
    """
    Builds a dictionary of the form {sn_name: velocity}
    """
    with open('../data/info_files/foley_master_data') as f:
        txt = f.readlines()

    clean = [x.strip() for x in txt]
    vel_dict = {}
    for entry in clean:
        ents = entry.split()
        if ents:
            vel_dict[ents[0]] = ents[2]
    return vel_dict


def build_gas_dict():
    """
    Builds a dictionary of the form {sn_name: gas rich (0/1)}
    """
    with open('../data/info_files/Gas-Rich-Poor.txt') as f:
        txt = f.readlines()
        clean = [x.strip() for x in txt]
        gas_dict = {}
        for entry in clean:
            ents = entry.split()
            if ents:
                gas_dict[ents[0]] = ents[1]
    return gas_dict


def build_carbon_dict():
    """
    Builds a dictionary of the form {sn_name: carbon }
    """
    with open('../data/info_files/carbon_presence.txt') as f:
        txt = f.readlines()
        clean = [x.strip() for x in txt]
        carbon_dict = {}
        for entry in clean:
            if not entry.startswith('#'):
                ents = entry.split()
                if ents:
                    carbon_dict[ents[0].lower()] = ents[1]
    return carbon_dict

def build_residual_dict():
    """
    Builds a dictionary of the form {sn_name: hubble_residual }
    """
    with open('../data/info_files/ryan_av.txt') as f:
        txt = f.readlines()
        clean = [x.strip() for x in txt]
        residual_dict = {}
        for entry in clean:
            if not entry.startswith('#'):
                ents = entry.split()
                if ents:
                    residual_dict[ents[0].lower()] = ents[2]
    return residual_dict


def build_redshift_dict(bsnipdict, cfadict):
    """
    Creates a dictionary of the form snname: redshift from the cfa and bsnip data
    """
    rsd = {}
    for item in cfadict:
        rsd[item] = float(cfadict[item][0])
    for item in bsnipdict:
        if item not in rsd:
            rsd[item] = float(bsnipdict[item])
    return rsd


def create_short_bsnip_dict(bsnipdict):
    newdict = {}
    for k, v in bsnipdict.iteritems():
        newdict[k.split('-')[0]] = v[1]/c
    return newdict


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#build necessary dictionaries
sndict, date_dict = read_cfa_info('../data/spectra/cfa/cfasnIa_param.dat',
                                  '../data/spectra/cfa/cfasnIa_mjdspec.dat')
bsnip_vals = read_bsnip_data('obj_info_table.txt')
short_bsnip_dict = create_short_bsnip_dict(bsnip_vals)
rsd = build_redshift_dict(short_bsnip_dict, sndict)
morph_dict = build_morph_dict()
vel_dict = build_vel_dict()
gas_dict = build_gas_dict()
carbon_dict = build_carbon_dict()
residual_dict = build_residual_dict()
ts = time.clock()

#con = sq3.connect('SNe.db')
con = sq3.connect('SNe_3.db')

#make sure no prior table in db to avoid doubling/multiple copies of same data

##version 1
# con.execute("""DROP TABLE IF EXISTS Supernovae""")
# con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename
#                     TEXT PRIMARY KEY, SN Text, Source Text, Redshift REAL,
#                     Phase REAL, MinWave REAL, MaxWave REAL, Dm15 REAL,
#                     M_B REAL, B_mMinusV_m REAL, Velocity REAL,
#                     Morphology INTEGER, Carbon TEXT, GasRich INTEGER, snr REAL,
#                     Interpolated_Spectra BLOB)""")

##version 2
# con.execute("""DROP TABLE IF EXISTS Supernovae""")
# con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename
#                     TEXT PRIMARY KEY, SN Text, Source Text, Redshift REAL,
#                     Phase REAL, MinWave REAL, MaxWave REAL, Dm15 REAL,
#                     M_B REAL, B_mMinusV_m REAL, Velocity REAL,
#                     Morphology INTEGER, Carbon TEXT, GasRich INTEGER, snr REAL, 
#                     Hubble_Res Real, Interpolated_Spectra BLOB)""")

con.execute("""DROP TABLE IF EXISTS Supernovae""")
con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename
                    TEXT PRIMARY KEY, SN Text, Source Text, Redshift REAL,
                    Phase REAL, MinWave REAL, MaxWave REAL, Dm15 REAL,
                    M_B REAL, B_mMinusV_m REAL, Velocity REAL,
                    Morphology INTEGER, Carbon TEXT, GasRich INTEGER, snr REAL, 
                    Hubble_Res Real, Interpolated_Spectra BLOB, Photometry BLOB)""")

#read all bsnip to find corrected
corr_list = []
allspec = []
for paths, subdirs, files in os.walk('../data/spectra/bsnip'):
    for name in files:
        if name.endswith('.flm'):
            allspec.append(name)
            if 'corrected' in name:
                corr_list.append(name)
clean1 = [re.sub('\-corrected.flm', '', x) for x in corr_list]
clean2 = [re.sub('\.flm', '', x) for x in allspec]
#find files that have both corrected and raw
res = set(clean1).intersection(set(clean2))
#have_both contains files to ignore
have_both = []
for s in res:
    s += '.flm'
    have_both.append(s)

#change this depending on where script is
root = '../data/spectra'
bad_files = []
bad_interp = []
shiftless = []

#ignore bad files
spectra_ignore = [
    'SN06mr_061113_g01_BAA_IM.dat',
    'SN07af_070310_g01_BAA_IM.dat',
    'SN07ai_070310_g01_BAA_IM.dat',
    'SN07af_070310_g01_BAA_IM.dat',
    'SN07ai_070312_g01_BAA_IM.dat',
    'sn1996ai-19960620.17-fast.flm',
    'sn1994T-19940612.25-mmt.flm',
    'sn1994T-19940613.19-fast.flm',
    'sn2000dn-20001006.25-fast.flm',
    'sn2006bt-20060502.33-fast.flm',
    'sn2006bt-20060503.33-fast.flm',
    'sn1994Q-19940612.32-mmt.flm',
    'sn2006H-20060206.21-fast.flm',
    'sn2002dj-20020614.17-fast.flm',
    'sn1999ek-19991030.47-fast.flm',
    'sn2003cq-20030408.26-fast.flm',
    'sn1996ai-19960620.17-fast.flm',
    ]
print "Adding information to table"
count = 1
for path, subdirs, files in os.walk(root):
    for name in files:
        #ignore bad spectra and uncorrected where corrected exist
        if name in spectra_ignore or name in have_both:
            continue

        #reset all fields to prevent carrying over from last loop iteration
        Dm15 = None
        m_b = None
        bm_vm = None
        phase = None
        redshift = None
        source = None

        f = os.path.join(path, name)
        if f.endswith('.flm') or f.endswith('.dat'):
            if 'cfasnIa' in f:
                continue
            try:
                if 'csp' in path:
                    spectra, info = read_csp(f)
                    sn_name = find_SN(name, 'csp', info)
                else:
                    spectra = read_cfa_or_bsnip_or_uv(f)
                    if 'other' in path:
                        sn_name = find_SN(name, 'other')
                    else:
                        sn_name = find_SN(name)
            except:
                bad_files.append(f)
                continue

            print count, sn_name
            count += 1

            #other source
            if 'other' in f:
                print f
                source = 'other'
                if sn_name in rsd:
                    redshift = rsd[sn_name]

            #cfa source
            elif 'cfa' in f:
                if 'sn2011' not in name:
                    sn_cfa = sndict[sn_name]
                else:
                    snd = None
                    sn_cfa = [None] * 14

            #csp source
            if 'csp' in f:
                print f
                source = 'csp'
                redshift = float(info[2])
                phase = float(info[4]) - float(info[3])

            #uv source
            elif 'uv' in path:
                print f
                print name
                source = 'uv'
                if sn_name in rsd:
                    redshift = rsd[sn_name]

            #cfa spectra
            elif 'cfa' in f:
                print f
                source = 'cfa'
                redshift = float(sn_cfa[0])
                #if statements assign values to fields if they exist, else none is assigned
                if sn_cfa[1] == '99999.9':
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
            elif 'bsnip' in f:
                source = 'bsnip'
                """
                certain files have words in their name that cause the
                split to be misalligned from the typical splits.
                """
                if is_number(name.split('-')[1][:8]):
                    longdate = name.split('-')[1][:8]
                else:
                    longdate = name.split('-')[2][:8]

                data = bsnip_vals[sn_name.lower()+'-'+longdate]
                print f
                if math.isnan(data[1]):
                    #skippy shitty redshiftless spectra
                    shiftless.append(name)
                    continue
                redshift = data[1]/c
                #if phase exists, add to db
                if not math.isnan(data[2]):
                    phase = data[2]

            waves = spectra[:, 0]
            min_wave = waves[0]
            max_wave = waves[len(spectra) - 1]
            spec = msg.packb(spectra)

            try:
                interp_spec, sig_noise = prep.compprep(spectra, sn_name, redshift, source)
            except Exception, e:
                #raise
                print e
                print "Interp failed"
                bad_interp.append(name)
                bad_files.append(name)
                interp_spec, sig_noise = None, None

            if sn_name in carbon_dict:
                carbon = carbon_dict[sn_name]
            elif sn_name == 'SNF20080909-030':
                carbon = carbon_dict['2008s5']
            elif sn_name == 'SNF20080514-002':
                carbon = carbon_dict['2008s1']
            elif sn_name == 'SNF20071021-000':
                carbon = carbon_dict['2007s1']
            else:
                carbon = None

            if sn_name in morph_dict:
                morph = morph_dict[sn_name]
            else:
                morph = None

            if sn_name in vel_dict:
                vel = vel_dict[sn_name]
            else:
                vel = None

            if sn_name in gas_dict:
                gasrich = gas_dict[sn_name]
            else:
                gasrich = None

            if sn_name in residual_dict:
                hubble_residual = residual_dict[sn_name]
            else:
                hubble_residual = None

            #add redshift to redshift dictionary if not already present
            if name not in rsd and redshift is not None:
                rsd[name] = redshift

            interped = msg.packb(interp_spec)

            sn_phot = phot.get_photometry(sn_name) 
            phot_blob = msg.packb(sn_phot)

            ##version 1
            # con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
            #                     Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
            #                     B_mMinusV_m, Velocity, Morphology, Carbon,
            #                     GasRich, snr, Interpolated_Spectra)
            #                     VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            #             (name, sn_name, source, redshift, phase,
            #              min_wave, max_wave, Dm15, m_b, bm_vm, vel,
            #              morph, carbon, gasrich, sig_noise, buffer(interped)))
            ## version 2
            # con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
            #                     Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
            #                     B_mMinusV_m, Velocity, Morphology, Carbon,
            #                     GasRich, snr, Hubble_Res, Interpolated_Spectra)
            #                     VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            #             (name, sn_name, source, redshift, phase,
            #              min_wave, max_wave, Dm15, m_b, bm_vm, vel,
            #              morph, carbon, gasrich, sig_noise, hubble_residual, buffer(interped)))
            ## version 3
            con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
                                Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
                                B_mMinusV_m, Velocity, Morphology, Carbon,
                                GasRich, snr, Hubble_Res, Interpolated_Spectra, Photometry)
                                VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                        (name, sn_name, source, redshift, phase,
                         min_wave, max_wave, Dm15, m_b, bm_vm, vel,
                         morph, carbon, gasrich, sig_noise, hubble_residual, buffer(interped), buffer(phot_blob)))

con.commit()
te = time.clock()
print 'bad files', bad_files, len(bad_files)
print 'bad interps', bad_interp, len(bad_interp)
print 'no available redshift', shiftless
print te - ts
