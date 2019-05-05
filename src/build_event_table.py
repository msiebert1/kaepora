from astropy.io import ascii
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import photometry as phot
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import scipy.interpolate as inter
import matplotlib.gridspec as gridspec
import pandas as pd

c = 299792.458

#there are a lot of redundant functions from build_spectral_table.py here

def find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,swift_dict,hst_dict,cfa_dict,short_bsnip_dict):
    events = []

    for SN in salt:
        events.append(SN['SimbadName'].lower())

    for SN in salt2:
        events.append(SN['SimbadName'].lower())

    for SN in mlcs31:
        events.append(SN['SimbadName'].lower())

    for SN in mlcs17:
        events.append(SN['SimbadName'].lower())

    for SN in lcparams:
        events.append(SN['SimbadName'].lower())

    for SN in host_data:
        events.append('sn' + SN['SN'].lower())

    for SN in av_dict.keys():
        events.append(SN.lower())

    for SN in cfa_dict.keys():
        events.append(SN.lower())

    for SN in short_bsnip_dict.keys():
        events.append(SN.lower())

    for SN in hst_dict.keys():
        events.append(SN.lower())

    for SN in swift_dict.keys():
        events.append(SN.lower())

    other_events = ['sn1989a', 'sn1989m', 'sn1993y', 'sn1993z', 'sn1999dg', 'sn1999do', 'sn2000dr', 'sn2000dx', 'sn2001dl', 'sn2001dt', 'sn2001dw', 'sn2002cv', 'sn2002eb', 
                    'sn2002eh', 'sn2002el', 'sn2003gs', 'sn2003gt', 'sn2003he', 'sn2003hs', 'sn2004bl', 'sn2004br', 'sn2004bv', 'sn2004bw', 'sn2004e', 'sn2004go', 'sn2005de', 
                    'sn2005di', 'sn2005dm', 'sn2005er', 'sn2005gj', 'sn2006dm', 'sn2007aj', 'sn2007cs', 'sn2007fr', 'sn2007ge', 'sn2007gi', 'sn2007gk', 'sn2007s1', 'sn2008bt', 
                    'sn2008cl', 'sn2008dt', 'sn2008dx', 'sn2008ec', 'sn2008ei', 'sn2008ha', 'sn2008hs', 'sn2008s1', 'sn2008s5', 'sn2005A', 'sn2005M', 'sn2005W', 'sn2006dd', 
                    'sn2006D', 'sn2006X', 'sn2007A', 'sn2007N', 'sn2007S', 'sn2008C', 'sn2008gl', 'sn2008hu', 'sn2008R', 'sn2009aa', 'sn2009ab', 'sn2009ad', 'sn2009ag', 
                    'sn2009D', 'sn2009F', 'sn2009Y','sn1994e', 'sn1994b', 'sn1998dj', 'sn1994j', 'sn2002ec', 'sn2005lt', 'snsnf20080522-011', 'sn1994u', 'sn2008bz', 'sn1998dw', 
                    'sn2002ep', 'sn1994x', 'sn2008fr', 'sn2004w', 'sn1997fc', 'sn2007e', 'sn2006es', 'sn2008fg', 'sn2008fj', 'sn2006eb', 'sn1995a', 'sn1995c', 'sn1999gf', 'sn2006ch', 
                    'sn2003x', 'sn1998en', 'sn1990g', 'snsnf20080522-000', 'sn1995t', 'sn1993aj', 'sn1993ai', 'sn2004gl', 'sn2005dh', 'sn2002dr', 'sn1993ab', 'sn2005do', 'sn2007v', 
                    'sn2006dy', 'sn2001fg', 'sn2006dw', 'sn2006dh', 'sn2006di', 'sn2006do', 'sn1995l', 'sn2008s4', 'sn2007m', 'sn1996o', 'sn2001a', 'sn2003v', 'sn2005bu', 'sn2002gx', 
                    'sn2002gg', 'sn2002gf', 'sn2002gc', 'sn2002gb', 'snsnf20080720-001', 'sn1996p', 'sn2005as', 'sn2004gw', 'sn2008dr', 'sn2008ds', 'sn2001cg', 'sn2000ej', 'sn1992ap', 
                    'sn2008db', 'sn2002fi', 'snsnf20080623-001', 'sn1998cl', 'sn1998cm', 'sn2006nr', 'sn1998cd', 'sn1997t', 'sn1999aq', 'sn1997fb', 'sn2004fw', 'sn1994ab', 'sn2005x', 
                    'sn2008gw', 'sn2007fq', 'sn2005p', 'sn2000dd', 'sn2004fg', 'sn2000df', 'sn2008ge', 'sn2008gh', 'sn2005f', 'sn2004y', 'sn2004bq', 'sn1991bj', 'sn1990m', 'sn1991bd', 
                    'sn1991bf', 'sn2003bf', 'sn2000j', 'sn1991bb', 'sn1991bc', 'sn2000q', 'sn1990r', 'snsnf20080514-002', 'sn2007ry', 'sn2001eg', 'sn1999c', 'sn2006ay', 'sn2008r3', 
                    'sn2004eq', 'sn2001ew', 'sn2001eu', 'sn2001er', 'sn2008bv', 'sn2008bw', 'sn2002dx', 'sn2003dw', 'sn1991k', 'sn2007sa', 'sn1991ak', 'sn1991b', 'sn1991am', 'sn1991at', 
                    'sn1991ay', 'sn2008ez', 'sn2008ey', 'sn2008er', 'sn2003p', 'sn2003ls', 'sn2008s3', 'sn2008gy', 'sn2002hl', 'sn2004di', 'sn2003lb', 'sn2008ee', 'sn2001dn', 'sn2003ax', 
                    'sn2003ay', 'sn1998fc', 'sn2003au', 'sn2003av', 'sn2003ij', 'sn2003ik', 'sn2005ec', 'sn2003an', 'sn2003ah', 'sn1999fz', 'sn1992m', 'sn2000al', 'sn2006ct', 'sn2004cv', 
                    'sn2003hj', 'sn2006ce', 'sn2004cb', 'sn2008hq', 'sn2008hr', 'sn2007qd', 'sn2008hy', 'sn2008hz', 'sn2002jo', 'sn2002bp', 'sn1993c', 'sn2002bk', 'sn2003bi', 'sn2003bh', 
                    'sn2008hj', 'sn2008hk', 'sn2004bz', 'sn2001dd', 'sn2001dm', 'sn2008cd', 'sn2008cf', 'sn2008ca', 'sn2008cb', 'sn2001ds', 'sn2008ct', 'ASASSN-14lp','sn2007sr','sn2008Q',
                    'sn2009an','sn2009dc','sn2009ig','sn2009Y','sn2010ev','sn2011aa','sn2011ao','sn2011by','sn2011fe','sn2011iv','sn2012cg','sn2012dn','sn2012fr','sn2012ht','sn2013aa',
                    'sn2013cg','sn2014J','sn2016ccz','sn2016coj','sn2016eiy','sn2016ekg','sn2016eoa','sn2016fff','sn2016gsb','sn2016itd','sn2017cbv','sn2017erp']

    for event in other_events:
        events.append(event)

    events = set(events)
    return events

def ReadExtin(file):
    #table containing B and V values for determining extinction -> dereddening due to milky way
    sne = np.genfromtxt(file, dtype=None)
    return sne

def MW_extinction_dict():
    sne_cfa = ReadExtin('extinction.dat')
    sne_bsnip = ReadExtin('extinctionbsnip.dat')
    sne_csp = ReadExtin('extinctioncsp.dat')
    sne_uv = ReadExtin('extinctionuv.dat')
    sne_other = ReadExtin('extinctionother.dat')
    sne_swift = ReadExtin('extinctionswiftuv.dat')
    sne_hst = ReadExtin('extinctionhst.dat')

    ext_dict = {}
    R_V = 3.1
    for j in range(len(sne_cfa)):
        sn = sne_cfa[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = R_V*(float(sne_cfa[j][1]) - float(sne_cfa[j][2]))
    for j in range(len(sne_bsnip)):
        sn = sne_bsnip[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = R_V*(float(sne_bsnip[j][1]) - float(sne_bsnip[j][2]))
    for j in range(len(sne_csp)):
        sn = sne_csp[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = R_V*(float(sne_csp[j][1]) - float(sne_csp[j][2]))
    for j in range(len(sne_uv)):
        sn = sne_uv[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = R_V*(float(sne_uv[j][1]) - float(sne_uv[j][2]))
    for j in range(len(sne_other)):
        sn = sne_other[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = R_V*(float(sne_other[j][1]) - float(sne_other[j][2]))
    for j in range(len(sne_swift)):
        if sne_swift[j][0] != 'ASASSN-14lp':
            sn = sne_swift[j][0].lower()
        else:
            sn = sne_swift[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = float(sne_swift[j][1])
    for j in range(len(sne_hst)):
        if sne_hst[j][0] != 'ASASSN-14lp':
            sn = 'sn' + sne_hst[j][0].lower()
        else:
            sn = sne_hst[j][0].lower()
        if sn not in ext_dict.keys():
            ext_dict[sn] = float(sne_hst[j][1])

    return ext_dict

def build_salt_dict(salt):
    salt_dict = {}
    for SN in salt:
        ra = SN['RA']
        dec = SN['DEC']
        zCMB = float(SN['zCMB'])
        e_zCMB = float(SN['e_zCMB'])
        Bmag = float(SN['Bmag'])
        e_Bmag = float(SN['e_Bmag'])
        s = float(SN['s'])
        e_s = float(SN['e_s'])
        c = float(SN['c'])  
        e_c = float(SN['e_c'])
        mu = float(SN['mu'])
        e_mu = float(SN['e_mu'])

        salt_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,Bmag,e_Bmag,s,e_s,c,e_c,mu,e_mu]

    return salt_dict

def build_salt2_dict(salt2):
    salt2_dict = {}
    for SN in salt2:
        ra = SN['RA']
        dec = SN['DEC']
        zCMB = float(SN['zCMB'])
        e_zCMB = float(SN['e_zCMB'])
        Bmag = float(SN['Bmag'])
        e_Bmag = float(SN['e_Bmag'])
        c = float(SN['c'])  
        e_c = float(SN['e_c'])
        mu = float(SN['mu'])
        e_mu = float(SN['e_mu'])
        x1 = float(SN['x1'])
        e_x1 = float(SN['e_x1'])

        salt2_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,Bmag,e_Bmag,x1,e_x1,c,e_c,mu,e_mu]
        
    return salt2_dict

def build_salt2_hub_res_dict(salt2_hub_res):
    salt2_hub_res_dict = {}
    for SN in salt2_hub_res:
        res = float(SN['MURES'])
        salt2_hub_res_dict['sn'+SN['CID'].lower()] = res
        
    return salt2_hub_res_dict

def build_mlcs31_dict(mlcs31):
    mlcs31_dict = {}
    for SN in mlcs31:
        ra = SN['RA']
        dec = SN['DEC']
        zCMB = float(SN['zCMB'])
        e_zCMB = float(SN['e_zCMB'])
        mu = float(SN['mu'])
        e_mu = float(SN['e_mu'])
        delta = float(SN['Delta'])
        e_delta = float(SN['e_Delt'])
        av = float(SN['AV'])
        e_av = float(SN['e_AV'])

        mlcs31_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,mu,e_mu,delta,e_delta,av,e_av]
        
    return mlcs31_dict

def build_mlcs17_dict(mlcs17):
    mlcs17_dict = {}
    for SN in mlcs17:
        ra = SN['RA']
        dec = SN['DEC']
        zCMB = float(SN['zCMB'])
        e_zCMB = float(SN['e_zCMB'])
        mu = float(SN['mu'])
        e_mu = float(SN['e_mu'])
        delta = float(SN['Delta'])
        e_delta = float(SN['e_Delt'])
        av = float(SN['AV'])
        e_av = float(SN['e_AV'])

        mlcs17_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,mu,e_mu,delta,e_delta,av,e_av]
        
    return mlcs17_dict

def build_host_dict(host_data):
    host_dict = {}
    for SN in host_data:
        ra = SN['RA']
        dec = SN['DEC']
        glon = float(SN['GLON'])
        glat = float(SN['GLAT'])
        cz = float(SN['cz'])
        czLG = float(SN['czLG'])
        czCMB = float(SN['czCMB'])
        mtype = SN['MType']
        xpos = float(SN['Xpos'])
        ypos = float(SN['Ypos'])
        t1 = float(SN['t1'])
        filt = SN['Filt']
        Ebv = float(SN['E(B-V)'])
        
        host_dict['sn' + SN['SN'].lower()] = [ra,dec,glon,glat,cz,czLG,czCMB,mtype,xpos,ypos,t1,filt,Ebv]
        
    return host_dict

def build_lcparams_dict(lcparams):
    lcparams_dict = {}
    for SN in lcparams:
        ra = SN['RA']
        dec = SN['DEC']
        zCMB = float(SN['zcmb'])
        zhel = float(SN['zhel'])
        mb = float(SN['mb'])
        e_mb = float(SN['e_mb'])
        c = float(SN['c'])  
        e_c = float(SN['e_c'])
        x1 = float(SN['x1'])
        e_x1 = float(SN['e_x1'])
        logMst = float(SN['logMst'])
        e_logMst = float(SN['e_log'])
        tmax = float(SN['tmax'])
        e_tmax = float(SN['e_tmax'])
        cov_mb_s = float(SN['cov_mb_s'])
        cov_mb_c = float(SN['cov_mb_c'])
        cov_s_c = float(SN['cov_s_c'])
        bias = float(SN['bias'])

        lcparams_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,zhel,mb,e_mb,c,e_c,x1,e_x1,logMst,e_logMst,tmax,e_tmax,cov_mb_s,cov_mb_c,cov_s_c,bias]
        
    return lcparams_dict

def build_av_dict(file):
     with open(file) as f:
        lines = f.readlines()

        av_dict = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                if 'sn' + l[1].lower() not in av_dict.keys():
                    av_dict['sn' + l[1].lower()] = float(l[18])

     return av_dict

def build_delt_dict(file):
     with open(file) as f:
        lines = f.readlines()
        delt_dict = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                delt_dict['sn' + l[1].lower()] = [float(l[16]), float(l[17])]
     return delt_dict

def build_cfa_dict(data_file):
    with open(data_file) as data:
        lines = data.readlines()
        cfa_dict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                cfa_dict['sn' + sndata[0].lower()] = sndata[1:]
    return cfa_dict

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
                    # date_dict[cleandate[0]] = cleandate[1]
                    date_dict[cleandate[0]] = cleandate[1]

    with open(data_file) as data:
        lines = data.readlines()
        sndict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                sndict['sn'+ sndata[0].lower()] = sndata[1:]

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
        # bsnip_dict[full_key] = [key1, bsnip_dict_rs[key1], line[len(line)-2]]
        bsnip_dict[full_key] = [key1, bsnip_dict_rs[key1], line[len(line)-2], line[len(line)-3]]
    return bsnip_dict

def create_short_bsnip_dict(bsnipdict):
    newdict = {}
    for k, v in bsnipdict.iteritems():
        newdict['sn' + k.split('-')[0]] = v[1]/c
    return newdict

def build_redshift_dict(bsnipdict, cfadict):
    """
    Creates a dictionary of the form snname: redshift from the cfa and bsnip data
    """
    rsd = {}
    for item in cfadict:
        rsd[item.lower()] = float(cfadict[item][0])
    for item in bsnipdict:
        if item not in rsd:
            rsd[item.lower()] = float(bsnipdict[item])
    return rsd

def build_NED_redshift_dict(NED_file):
    with open(NED_file) as data:
        lines = data.readlines()
        NED_red_dict = {}
        for line in lines:
            sn_name = line.split()[0]
            redshift = float(line.split()[1])
            NED_red_dict['sn'+sn_name.lower()] = redshift

    return NED_red_dict

def build_swift_dict(data_file):
    with open(data_file) as data:
        lines = data.readlines()
        swift_dict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                if sndata[0].lower().startswith('sn'):
                    swift_dict[sndata[0].lower()] = [sndata[1],sndata[3]]
                else:
                    swift_dict[sndata[0]] = [sndata[1],sndata[3]]

    return swift_dict

def build_hst_dict(data_file):
    with open(data_file) as data_meta:
        data_dict = {}
        for line in data_meta.readlines():
            if line.split()[0] != 'asassn-14lp':
                data_dict['sn' + line.split()[0].lower()] = line.split()[1:]
            else:
                data_dict[line.split()[0].lower()] = line.split()[1:]
    return data_dict

def build_NED_host_dict(data_file):
    with open(data_file) as data:
        lines = data.readlines()
        NED_host_dict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                if sndata[1] != 'None':
                    NED_host_dict['sn' + sndata[0].lower()] = sndata[1]
                else:
                    NED_host_dict['sn' + sndata[0].lower()] = None
    return NED_host_dict

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
            # vel_dict[ents[0]] = ents[2]
            vel_dict['sn'+ ents[0].lower()] = [ents[2], ents[3]]
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
                # gas_dict[ents[0]] = ents[1]
                gas_dict['sn'+ ents[0].lower()] = ents[1]
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
                    carbon_dict['sn'+ ents[0].lower()] = ents[1]
    return carbon_dict

def build_residual_dict():
    """
    Builds a dictionary of the form {sn_name: hubble_residual }
    """
    with open('../data/info_files/ryan_res.txt') as f:
        txt = f.readlines()
        clean = [x.strip() for x in txt]
        residual_dict = {}
        for entry in clean:
            if not entry.startswith('#'):
                ents = entry.split()
                if ents:
                    residual_dict['sn'+ ents[0].lower()] = ents[2]
    return residual_dict

def build_tmax_dict(cfa_dict):
    #NOTES: 
    #csp from spectra files convert to mjd if no mlcs
    #cfa from cfasnIa_param.dat if no mlcs
    #bsnip no data, phase from spec_info_table.txt if no mlcs
    #uv and other no data
    #prioritize mlcs from lowzrv25_all.fitres

    # here compile csp, cfa, and mlcs into dict
    with open('../data/info_files/csp_mjd_max.txt') as file:
        lines=file.readlines()
        result=[]
        tmax_dict_csp = {}
        for x in lines:
            tmax_dict_csp[x.split()[0]] = float(x.split()[1])
    
    tmax_dict_cfa = {}
    for sn in cfa_dict:
        if cfa_dict[sn][1] != '99999.9':
            tmax_dict_cfa[sn] = float(cfa_dict[sn][1])

    with open('../data/info_files/lowz_rv25_all.fitres') as f:
        lines = f.readlines()

        tmax_dict_mlcs25 = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                if 'sn' + l[1].lower() not in tmax_dict_mlcs25.keys():
                    tmax_dict_mlcs25['sn' + l[1].lower()] = float(l[14])

    tmax_dict = {}
    for t in tmax_dict_mlcs25:
        tmax_dict[t] = tmax_dict_mlcs25[t]

    for t in tmax_dict_cfa:
        if t not in tmax_dict.keys():
            tmax_dict[t] = tmax_dict_cfa[t]

    for t in tmax_dict_csp:
        if t not in tmax_dict.keys():
            tmax_dict[t] = tmax_dict_csp[t]
    return tmax_dict


def dm15_from_fit_params(events, fit_dict, cfa_dict, stretch='N/A'):
    dm15_arr = []
    e_dm15_arr = []
    param_arr = []
    e_param_arr = []
    event_arr = []
    cfa_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None]
    fit_none = [None,None,None,None,None,None,None,None,None,None]
    lowrv_none = [None,None]
    for event in events:
        if event != 'sn2001da' and event != 'sn2001cp': #outliers
            if event in fit_dict and event in cfa_dict:
                dm15_cfa = cfa_dict.get(event, cfa_none)[4]
                if dm15_cfa != '9.99' and dm15_cfa != None:
                    dm15_cfa = float(dm15_cfa)
                    dm15_arr.append(dm15_cfa)
                    e_dm15_arr.append(float(cfa_dict.get(event, cfa_none)[5]))
                    if stretch is 'delta_lowrv':
                        param_arr.append(fit_dict.get(event, lowrv_none)[0])
                        e_param_arr.append(fit_dict.get(event,lowrv_none)[1])
                        event_arr.append(event)
                    else:
                        param_arr.append(fit_dict.get(event, fit_none)[6]) #stretch param is 6th index of dicts
                        e_param_arr.append(fit_dict.get(event,fit_none)[7])
                        event_arr.append(event)

    # weights = 1./np.asarray(e_param_arr)
    # coeffs = np.polyfit(dm15_arr, param_arr, 2, w=weights)

    # x = np.linspace(np.amin(dm15_arr), np.amax(dm15_arr), 10000)
    # y = coeffs[0]*x**2. + coeffs[1]*x + coeffs[2]
    # print coeffs
    # print y[0], y[-1]

    weights = 1./np.asarray(e_dm15_arr)
    coeffs = np.polyfit(param_arr, dm15_arr, 2, w=weights)

    # x = np.linspace(np.amin(param_arr), np.amax(param_arr), 10000)
    if stretch is 's':
        x = np.linspace(.2, 1.3, 10000)
    if stretch is 'x1':
        x = np.linspace(-6., 4., 10000)
    if stretch is 'delta' or stretch is 'delta_lowrv':
        x = np.linspace(-1., 2., 10000)
    y = coeffs[0]*x**2. + coeffs[1]*x + coeffs[2]
    # print coeffs
    # print y[0], y[-1]

    nsig = 3.
    clip_param_arr = []
    clip_e_param_arr = []
    clip_dm15_arr = []
    clip_e_dm15_arr = []
    sig_away = []
    resids = []
    inv_resids = []
    inv_sig_away = []
    clip_event_arr = []

    # interp_func = interp1d(y, x, bounds_error=False, fill_value=None)
    # inter_vals = interp_func(param_arr)

    interp_func = interp1d(x, y, bounds_error=False, fill_value=None)
    inter_vals = interp_func(param_arr)
    # plt.plot(x, y, 'ro')
    # plt.plot(param_arr, dm15_arr, 'bo')
    # plt.plot(param_arr, inter_vals, 'go')
    # plt.show()

    inv_interp_func = interp1d(y, x, bounds_error=False, fill_value="extrapolate")
    inv_inter_vals = inv_interp_func(dm15_arr)

    for i, p in enumerate(dm15_arr):
        resids.append(p - inter_vals[i])
        sig_away.append(np.absolute(p - inter_vals[i])/e_dm15_arr[i])

        inv_resids.append(param_arr[i] - inv_inter_vals[i])
        inv_sig_away.append(np.absolute(param_arr[i] - inv_inter_vals[i])/e_param_arr[i])

        if np.absolute(p - inter_vals[i])/e_dm15_arr[i] > nsig and np.absolute(param_arr[i] - inv_inter_vals[i])/e_param_arr[i] > nsig:
        # if np.absolute(p - inter_vals[i])/e_dm15_arr[i] > nsig:
            clip_dm15_arr.append(p)
            clip_e_dm15_arr.append(e_dm15_arr[i])
            clip_param_arr.append(param_arr[i])
            clip_e_param_arr.append(e_param_arr[i])
            clip_event_arr.append(event_arr[i])

    new_param_arr = []
    new_e_param_arr = []
    new_dm15_arr = []
    new_e_dm15_arr = []
    new_resids = []
    # print inv_sig_away
    for i, p in enumerate(dm15_arr):
        if sig_away[i] < nsig or inv_sig_away[i] < nsig:
        # if sig_away[i] < nsig:
            new_dm15_arr.append(p)
            new_e_dm15_arr.append(e_dm15_arr[i])
            new_param_arr.append(param_arr[i])
            new_e_param_arr.append(e_param_arr[i])
            new_resids.append(resids[i])

    se = ((1./np.asarray(new_e_dm15_arr))*(np.asarray(new_resids)**2.))/(np.average(1./np.asarray(new_e_dm15_arr)))
    # se = (np.asarray(new_resids)**2.)
    mse = np.sum(se)/(len(new_resids) - 3.)
    rmse = np.sqrt(mse)
    # print 'RMSE: ', rmse
    # print 'Outlier Fraction: ', float(len(clip_dm15_arr))/float(len(dm15_arr))
    
    # RMSE:  0.0718342630003
    # Outlier Fraction:  0.0752688172043
    # [ 1.7240621  -4.66995066  3.93727553]
    # sn2006et
    # sn1994t
    # sn1995ac
    # sn2000dk
    # sn2002bf
    # sn2005ls
    # sn1999aa

    # RMSE:  0.066249618878
    # Outlier Fraction:  0.0504201680672
    # [ 0.01828958 -0.13430543  1.02585001]
    # sn1994t
    # sn1995ac
    # sn2002dj
    # sn1999dq
    # sn2007bz
    # sn1999aa

    # RMSE:  0.070183174903
    # Outlier Fraction:  0.0230769230769
    # [-0.16315158  0.76583297  1.12697325]
    # sn2003iv
    # sn1995ac
    # sn1999aa

    # RMSE:  0.07406417245
    # Outlier Fraction:  0.139423076923
    # [-0.04957799  0.56115546  1.11213983]
    # sn2007af
    # sn2006eu
    # sn1997e
    # sn1999by
    # sn2003iv
    # sn2006x
    # sn2006n
    # sn1999ej
    # sn1995ac
    # sn1995al
    # sn2005be
    # sn2002dj
    # sn2004dt
    # sn2000dk
    # sn2006hb
    # sn2005el
    # sn2000cx
    # sn2007if
    # sn2005hk
    # sn2002cx
    # sn2001fh
    # sn2006al
    # sn1997br
    # sn2005lu
    # sn2006cm
    # sn2006cj
    # sn1999aa
    # sn1999ac
    # sn2006bz

    # weights = 1./np.asarray(new_e_param_arr)
    # coeffs = np.polyfit(new_dm15_arr, new_param_arr, 2, w=weights)
    weights = 1./np.asarray(new_e_dm15_arr)
    coeffs = np.polyfit(new_param_arr, new_dm15_arr, 2, w=weights)
    # new_x = np.linspace(np.amin(new_dm15_arr), np.amax(new_dm15_arr), 10000)
    # new_y = coeffs[0]*new_x**2. + coeffs[1]*new_x + coeffs[2]
    new_x = np.linspace(np.amin(new_param_arr), np.amax(new_param_arr), 10000)
    new_y = coeffs[0]*new_x**2. + coeffs[1]*new_x + coeffs[2]
    # print coeffs
    # for e in clip_event_arr:
    #     print e
    # print 
    # print new_y[0], new_y[-1]

    # plotted in lc_fitting_relationships.py
    # plt.rc('font', family='serif')
    # fig, ax = plt.subplots(1,1)
    # fig.set_size_inches(10, 8, forward = True)
    # plt.minorticks_on()
    # plt.xticks(fontsize = 20)
    # plt.yticks(fontsize = 20)
    # plt.tick_params(
    #     which='major', 
    #     bottom='on', 
    #     top='on',
    #     left='on',
    #     right='on',
    #     length=10)
    # plt.tick_params(
    #   which='minor', 
    #   bottom='on', 
    #   top='on',
    #   left='on',
    #   right='on',
    #   length=5)
    # plt.plot(y, x, 'k', linewidth=4)
    # plt.plot(new_y, new_x, 'g', linewidth=4)
    # plt.plot(x, y, 'k', linewidth=4)
    # plt.plot(new_x, new_y, 'k', linewidth=4)
    # # plt.ylim(.6,2.2)
    # plt.ylabel('$\Delta m_{15}$ (B) (mag)', fontsize = 30)
    # if stretch is 's':
    #   plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, yerr= new_e_dm15_arr, fmt='o', color='#7570b3', ms=10)
    #   # plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, fmt='o', color='#7570b3', ms=10)
    #   plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, yerr=clip_e_dm15_arr, fmt='x', color='black', ms=10, mew=1)
    #   # plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, fmt='x', color='black', ms=30, mew=3)
    #   plt.xlim(.4, 1.3)
    #   plt.ylim(.6, 2.5)
    #   plt.xlabel('s', fontsize = 30)
    #   # plt.savefig('../../../Paper_Drafts/reprocessing_updated/dm15_s.pdf', dpi = 300, bbox_inches = 'tight')
    # elif stretch is 'x1':
    #   plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, yerr= new_e_dm15_arr, fmt='o', color='#1b9e77', ms=10)
    #   # plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, fmt='o', color='#1b9e77', ms=10)
    #   plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, yerr=clip_e_dm15_arr, fmt='x', color='black', ms=10, mew=1)
    #   # plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, fmt='x', color='black', ms=30, mew=3)
    #   plt.xlim(-5., 3.)
    #   plt.ylim(.6, 2.5)
    #   plt.xlabel('$x_1$', fontsize = 30)
    #   # plt.savefig('../../../Paper_Drafts/reprocessing_updated/dm15_x1.pdf', dpi = 300, bbox_inches = 'tight')
    # elif stretch is 'delta':
    #   plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, yerr= new_e_dm15_arr, fmt='o', color='#d95f02', ms=10)
    #   # plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, fmt='o', color='#d95f02', ms=10)
    #   plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, yerr=clip_e_dm15_arr, fmt='x', color='black', ms=10, mew=1)
    #   # plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, fmt='x', color='black', ms=30, mew=3)
    #   plt.xlim(-.5, 1.8)
    #   plt.ylim(.6, 2.5)
    #   plt.xlabel('$\Delta$', fontsize = 30)
    #   # plt.savefig('../../../Paper_Drafts/reprocessing_updated/dm15_delta.pdf', dpi = 300, bbox_inches = 'tight')
    # elif stretch is 'delta_lowrv':
    #   print len(new_dm15_arr)
    #   plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, yerr= new_e_dm15_arr, fmt='o', color='#d95f02', ms=10)
    #   # plt.errorbar(new_param_arr, new_dm15_arr, xerr=new_e_param_arr, fmt='o', color='#d95f02', ms=10)
    #   plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, yerr=clip_e_dm15_arr, fmt='x', color='black', ms=10, mew=1)
    #   # plt.errorbar(clip_param_arr, clip_dm15_arr, xerr=clip_e_param_arr, fmt='x', color='black', ms=30, mew=3)
    #   plt.xlim(-.9, 2.0)
    #   plt.ylim(.6, 2.5)
    #   plt.xlabel('$\Delta$', fontsize = 30)
        # plt.savefig('../../../Paper_Drafts/reprocessing_updated/dm15_delta_lowrv.pdf', dpi = 300, bbox_inches = 'tight')
    # plt.show()

    # dm15_interp = interp1d(x, y, bounds_error = True)
    dm15_interp = interp1d(new_x, new_y, bounds_error=False, fill_value=None)

    ###
    # s
    # RMSE:  0.0718342630003
    # Outlier Fraction:  0.0752688172043
    # [ 1.7240621  -4.66995066  3.93727553]

    # x1
    # RMSE:  0.066249618878
    # Outlier Fraction:  0.0504201680672
    # [ 0.01828958 -0.13430543  1.02585001]

    # delta
    # RMSE:  0.070183174903
    # Outlier Fraction:  0.0230769230769
    # [-0.16315158  0.76583297  1.12697325]

    # delta rv=2.5
    # RMSE:  0.07406417245
    # Outlier Fraction:  0.139423076923
    # [-0.04957799  0.56115546  1.11213983]
    return dm15_interp, (new_x, new_y, new_param_arr, new_dm15_arr, new_e_param_arr, 
                        new_e_dm15_arr, clip_param_arr, clip_dm15_arr, clip_e_param_arr, clip_e_dm15_arr)

def plot_dm15_fits(s_results, x1_results, delta_results, delta_lowrv_results):
    gs = gridspec.GridSpec(1, 4, wspace = .0003)
    fig = plt.figure(num = 1, dpi = 100, figsize = [23,7], facecolor = 'w')

    s_plot = plt.subplot(gs[0])
    plt.rc('font', family='serif')
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=10)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=5)
    plt.plot(s_results[0], s_results[1], 'k', linewidth=4)
    plt.ylabel(r'$\Delta m_{15}$ (B) (mag)', fontsize = 30)
    plt.errorbar(s_results[2], s_results[3], xerr=s_results[4], yerr= s_results[5], fmt='o', color='#7570b3', ms=10)
    plt.errorbar(s_results[6], s_results[7], xerr=s_results[8], yerr=s_results[9], fmt='x', color='black', ms=10, mew=1)
    plt.xlim(.4, 1.25)
    plt.ylim(.6, 2.5)
    plt.xlabel('s', fontsize = 30)

    x1_plot = plt.subplot(gs[1], sharey=s_plot)
    plt.rc('font', family='serif')
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=10)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=5)
    plt.plot(x1_results[0], x1_results[1], 'k', linewidth=4)
    plt.errorbar(x1_results[2], x1_results[3], xerr=x1_results[4], yerr= x1_results[5], fmt='o', color='#1b9e77', ms=10)
    plt.errorbar(x1_results[6], x1_results[7], xerr=x1_results[8], yerr=x1_results[9], fmt='x', color='black', ms=10, mew=1)
    plt.xlim(-4.9, 2.9)
    plt.ylim(.6, 2.5)
    plt.xlabel('$x_1$', fontsize = 30)
    plt.setp(x1_plot.get_yticklabels(), visible=False)

    delta_plot = plt.subplot(gs[2], sharey=s_plot)
    plt.rc('font', family='serif')
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=10)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=5)
    plt.plot(delta_results[0], delta_results[1], 'k', linewidth=4)
    plt.errorbar(delta_results[2], delta_results[3], xerr=delta_results[4], yerr= delta_results[5], fmt='o', color='#d95f02', ms=10)
    plt.errorbar(delta_results[6], delta_results[7], xerr=delta_results[8], yerr= delta_results[9], fmt='x', color='black', ms=10, mew=1)
    plt.xlim(-.85, 1.95)
    plt.ylim(.6, 2.5)
    plt.xlabel(r'$\Delta$', fontsize = 30)
    plt.setp(delta_plot.get_yticklabels(), visible=False)

    delta_lowrv_plot = plt.subplot(gs[3], sharey=s_plot)
    plt.rc('font', family='serif')
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=10)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        length=5)
    plt.plot(delta_lowrv_results[0], delta_lowrv_results[1], 'k', linewidth=4)
    plt.errorbar(delta_lowrv_results[2], delta_lowrv_results[3], xerr=delta_lowrv_results[4], yerr= delta_lowrv_results[5], fmt='o', color='#d95f02', ms=10)
    plt.errorbar(delta_lowrv_results[6], delta_lowrv_results[7], xerr=delta_lowrv_results[8], yerr= delta_lowrv_results[9], fmt='x', color='black', ms=10, mew=1)
    plt.xlim(-.85, 2.0)
    plt.ylim(.6, 2.5)
    plt.xlabel(r'$\Delta$', fontsize = 30)
    plt.setp(delta_lowrv_plot.get_yticklabels(), visible=False)
    plt.savefig('../../../Paper_Drafts/reprocessing_updated/dm15_from_fits.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()

def main():
    mn.patch()

    salt = ascii.read("../data/info_files/salt_params_dists.txt", delimiter = r'\s', guess = False)
    salt2 = ascii.read("../data/info_files/salt2_params_dists.txt", delimiter = r'\s', guess = False)
    salt2_hub_res = ascii.read("../data/info_files/SALT2mu_fpan.fitres", delimiter = r'\s', guess = False)
    mlcs31 = ascii.read("../data/info_files/mlcs31_params.txt", delimiter = r'\s', guess = False)
    mlcs17 = ascii.read("../data/info_files/mlcs17_params.txt", delimiter = r'\s', guess = False)
    lcparams = ascii.read("../data/info_files/lc_params.txt", delimiter = r'\s', guess = False)
    host_data = ascii.read("../data/info_files/other_host_data.txt", delimiter = r'\s', guess = False)

    av_dict = build_av_dict('../data/info_files/lowz_rv25_all.fitres')
    delt_dict = build_delt_dict('../data/info_files/lowz_rv25_all.fitres')
    # cfa_dict = build_cfa_dict('..\data\spectra\cfa\cfasnIa_param.dat') #REDO cfa_dict stuff
    cfa_dict, date_dict = read_cfa_info('../data/spectra/cfa/cfasnIa_param.dat', '../data/spectra/cfa/cfasnIa_mjdspec.dat')
    swift_dict = build_swift_dict('../data/spectra/swift_uvspec/swift_uv_log.txt')
    hst_dict = build_hst_dict('../data/spectra/hst_foley/hst_metadata.txt')

    NED_host_dict = build_NED_host_dict('../data/info_files/NED_host_info_simplified.txt')
    vel_dict = build_vel_dict()
    gas_dict = build_gas_dict()
    carbon_dict = build_carbon_dict()
    residual_dict = build_residual_dict()
    tmax_dict = build_tmax_dict(cfa_dict)

    #new redshift function moved from newdb
    bsnip_vals = read_bsnip_data('obj_info_table.txt')
    short_bsnip_dict = create_short_bsnip_dict(bsnip_vals)
    rsd = build_redshift_dict(short_bsnip_dict, cfa_dict)
    NED_red_dict = build_NED_redshift_dict('../data/info_files/NED_redshift_info.txt')

    events = find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,swift_dict,hst_dict,cfa_dict,short_bsnip_dict)

    salt_dict = build_salt_dict(salt)
    salt2_dict = build_salt2_dict(salt2)
    salt2_hub_res_dict = build_salt2_hub_res_dict(salt2_hub_res)
    mlcs31_dict = build_mlcs31_dict(mlcs31)
    mlcs17_dict = build_mlcs17_dict(mlcs17)
    host_dict = build_host_dict(host_data)
    lcparams_dict = build_lcparams_dict(lcparams)

    salt_set = set(salt_dict.keys())
    salt2_set = set(salt2_dict.keys())
    mlcs31_set = set(mlcs31_dict.keys())
    delt_set = set(delt_dict.keys())
    common_SNe = set.intersection(salt_set, salt2_set, delt_set)

    # need this to determine outlier fractions:
    # s: 7/85, x1: 3/85, delta: 2/85, delta25: 6/85
    # RMS from original fit:
    # s: .071, x1: .066, delta: .074
    # salt_dict_new = {}
    # salt2_dict_new = {}
    # mlcs31_dict_new = {}
    # delt_dict_new = {}
    # for cSN in common_SNe:
    #   salt_dict_new[cSN] = salt_dict[cSN]
    #   salt2_dict_new[cSN] = salt2_dict[cSN]
    #   mlcs31_dict_new[cSN] = mlcs31_dict[cSN]
    #   delt_dict_new[cSN] = delt_dict[cSN]

    # dm15_s_interp, s_results = dm15_from_fit_params(events, salt_dict_new, cfa_dict, stretch='s')
    # dm15_x1_interp, x1_results = dm15_from_fit_params(events, salt2_dict_new, cfa_dict, stretch='x1')
    # dm15_delta_interp, delta_results = dm15_from_fit_params(events, mlcs31_dict_new, cfa_dict, stretch='delta')
    # dm15_delta_lowrv_interp, delta_lowrv_fit_results = dm15_from_fit_params(events, delt_dict_new, cfa_dict, stretch='delta_lowrv')


    dm15_s_interp, s_results = dm15_from_fit_params(events, salt_dict, cfa_dict, stretch='s')
    dm15_x1_interp, x1_results = dm15_from_fit_params(events, salt2_dict, cfa_dict, stretch='x1')
    dm15_delta_interp, delta_fit_results = dm15_from_fit_params(events, mlcs31_dict, cfa_dict, stretch='delta')
    dm15_delta_lowrv_interp, delta_lowrv_fit_results = dm15_from_fit_params(events, delt_dict, cfa_dict, stretch='delta_lowrv')
    # plot_dm15_fits(s_results, x1_results, delta_fit_results, delta_lowrv_fit_results)
    ext_dict = MW_extinction_dict()
    # raise TypeError

    # con = sq3.connect('../data/kaepora_v1.db')
    con = sq3.connect('../data/kaepora_v1.db')
    con.execute("""DROP TABLE IF EXISTS Events""")
    con.execute("""CREATE TABLE IF NOT EXISTS Events (SN TEXT, RA TEXT, DEC TEXT, 
                                                          zCMB_salt REAL, e_zCMB_salt REAL, Bmag_salt REAL, e_Bmag_salt REAL, s_salt REAL, e_s_salt REAL, c_salt REAL, e_c_salt REAL, mu_salt REAL, e_mu_salt REAL,
                                                          zCMB_salt2 REAL, e_zCMB_salt2 REAL, Bmag_salt2 REAL, e_Bmag_salt2 REAL, x1_salt2 REAL, e_x1_salt2 REAL, c_salt2 REAL, e_c_salt2 REAL, mu_salt2 REAL, e_mu_salt2 REAL,
                                                          zCMB_mlcs31 REAL, e_zCMB_mlcs31 REAL, mu_mlcs31 REAL, e_mu_mlcs31 REAL, delta_mlcs31 REAL, e_delta_mlcs31 REAL, av_mlcs31 REAL, e_av_mlcs31 REAL,
                                                          zCMB_mlcs17 REAL, e_zCMB_mlcs17 REAL, mu_mlcs17 REAL, e_mu_mlcs17 REAL, delta_mlcs17 REAL, e_delta_mlcs17 REAL, av_mlcs17 REAL, e_av_mlcs17 REAL,
                                                          glon_host REAL, glat_host REAL, cz_host REAL, czLG_host REAL, czCMB_host REAL, mtype_host TEXT, xpos_host REAL, ypos_host REAL, t1_host REAL, filt_host TEXT, Ebv_host REAL,
                                                          zCMB_lc REAL, zhel_lc REAL, mb_lc REAL, e_mb_lc REAL, c_lc REAL, e_c_lc REAL, x1_lc REAL, e_x1_lc REAL, logMst_lc REAL, e_logMst_lc REAL, tmax_lc REAL, e_tmax_lc REAL, cov_mb_s_lc REAL, cov_mb_c_lc REAL, cov_s_c_lc REAL, bias_lc REAL,
                                                          
                                                          Av_MW REAL, Av_25 REAL, MJD_max REAL, Dm15_source REAL, Dm15_from_fits REAL, e_dm15 REAL, separation REAL, NED_host REAL, V_at_max REAl, V_err REAL,
                                                          Redshift REAL, M_b_cfa REAL, M_b_cfa_err REAL, B_minus_V_cfa REAL, B_minus_V_cfa_err REAL, Carbon_presence REAL, Na_presence REAL, Hubble_res REAL,
                                                          Photometry BLOB, csp_photometry BLOB)""")

    salt_none = [None,None,None,None,None,None,None,None,None,None,None,None]
    salt2_none = salt_none
    mlcs31_none = [None,None,None,None,None,None,None,None,None,None]
    mlcs17_none = mlcs31_none
    host_none = [None,None,None,None,None,None,None,None,None,None,None,None,None]
    lc_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
    cfa_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None]
    swift_none = [None,None]
    vel_none = [None,None]
    hst_none = [None,None,None,None,None,None,None]
    
    for event in events:
        print 'Adding data for ' + event + '...'

        #other parameters from metadata files
        ra_salt,dec_salt,zCMB_salt,e_zCMB_salt,Bmag_salt,e_Bmag_salt,s_salt,e_s_salt,c_salt,e_c_salt,mu_salt,e_mu_salt = salt_dict.get(event, salt_none)
        ra_salt2,dec_salt2,zCMB_salt2,e_zCMB_salt2,Bmag_salt2,e_Bmag_salt2,x1_salt2,e_x1_salt2,c_salt2,e_c_salt2,mu_salt2,e_mu_salt2 = salt2_dict.get(event, salt2_none)
        ra_mlcs31,dec_mlcs31,zCMB_mlcs31,e_zCMB_mlcs31,mu_mlcs31,e_mu_mlcs31,delta_mlcs31,e_delta_mlcs31,av_mlcs31,e_av_mlcs31 = mlcs31_dict.get(event, mlcs31_none)
        ra_mlcs17,dec_mlcs17,zCMB_mlcs17,e_zCMB_mlcs17,mu_mlcs17,e_mu_mlcs17,delta_mlcs17,e_delta_mlcs17,av_mlcs17,e_av_mlcs17 = mlcs17_dict.get(event, mlcs17_none)
        ra_host,dec_host,glon_host,glat_host,cz_host,czLG_host,czCMB_host,mtype_host,xpos_host,ypos_host,t1_host,filt_host,Ebv_host = host_dict.get(event, host_none)
        ra_lc,dec_lc,zCMB_lc,zhel_lc,mb_lc,e_mb_lc,c_lc,e_c_lc,x1_lc,e_x1_lc,logMst_lc,e_logMst_lc,tmax_lc,e_tmax_lc,cov_mb_s_lc,cov_mb_c_lc,cov_s_c_lc,bias_lc = lcparams_dict.get(event, lc_none)
        ras = [ra_salt,ra_salt2,ra_mlcs31,ra_mlcs17,ra_host,ra_lc]
        decs = [dec_salt,dec_salt2,dec_mlcs31,dec_mlcs17,dec_host,dec_lc]

        #host, position, and separation
        sep = None
        if xpos_host != None and ypos_host != None:
            sep = (float(xpos_host)**2 + float(ypos_host)**2)**(.5)

        ned_host = NED_host_dict.get(event, None)
        if ned_host is None:
            ned_host = hst_dict.get(event, hst_none)[6]
            if ned_host != None:
                ned_host = float(ned_host)

        ra = None
        for r in ras:
            if r != None:
                ra = r
                break

        dec = None
        for d in decs:
            if d != None:
                dec = d
                break

        #redshift
        redshift = NED_red_dict.get(event, None)
        if redshift == None:
            redshift = rsd.get(event, None)

        #M_B from cfa
        m_b_cfa = cfa_dict.get(event, cfa_none)[7]
        if m_b_cfa != '-99.99' and m_b_cfa != None:
            m_b_cfa = float(m_b_cfa)
            m_b_cfa_err = float(cfa_dict.get(event, cfa_none)[8])
        else:
            m_b_cfa = None
            m_b_cfa_err = None

        #B - V color at max from cfa
        b_minus_v_cfa = cfa_dict.get(event, cfa_none)[9]
        if b_minus_v_cfa != '-9.99' and b_minus_v_cfa != None:
            b_minus_v_cfa = float(b_minus_v_cfa)
            b_minus_v_cfa_err = float(cfa_dict.get(event, cfa_none)[10])
        else:
            b_minus_v_cfa = None
            b_minus_v_cfa_err = None

        #carbon presence
        carbon = carbon_dict.get(event, None)

        #na presence
        na = gas_dict.get(event, None)

        #hubble residual
        # hubble_res = residual_dict.get(event, None)
        hubble_res = salt2_hub_res_dict.get(event, None)

        #MW reddening
        av_mw = ext_dict.get(event, None)

        #MLCS reddening
        av_25 = av_dict.get(event)
        delta_lowrv = delt_dict.get(event)
        if av_mlcs31 is None:
            av_mlcs31 = hst_dict.get(event, hst_none)[4]
            if av_mlcs31 != None:
                av_mlcs31 = float(av_mlcs31)*3.1

        #mjd_max prioritizing mlcs25, cfa, then csp
        mjd_max = tmax_dict.get(event, None)

        #velocity at maximum light
        vel = vel_dict.get(event, vel_none)[0]
        if vel != '-99.0000' and vel != None:
            vel = float(vel)
            e_vel = float(vel_dict.get(event, vel_none)[1])
        else:
            vel = None
            e_vel = None

        if vel is None:
            vel = hst_dict.get(event, hst_none)[2]
            if vel != None:
                vel = float(vel)/1000.

        #dm15 estimation
        dm15_source = cfa_dict.get(event, cfa_none)[4]
        if dm15_source != '9.99' and dm15_source != None:
            dm15_source = float(dm15_source)
            e_dm15 = float(cfa_dict.get(event, cfa_none)[5])
        elif dm15_source == '9.99':
            dm15_source = None
            e_dm15 = None

        if dm15_source is None:
            dm15_source = swift_dict.get(event, swift_none)[1]
            if dm15_source != '-99' and dm15_source != None:
                dm15_source = float(dm15_source)
                e_dm15 = 0.
            elif dm15_source == '-99':
                dm15_source = None
                e_dm15 = None

        if dm15_source is None:
            dm15_source = hst_dict.get(event, hst_none)[0]
            e_dm15 = hst_dict.get(event, hst_none)[1]
            if dm15_source != None:
                dm15_source = float(dm15_source)
                e_dm15 = float(e_dm15)


        if dm15_source is None:
            dm15_from_s = np.nan
            if s_salt != None:
                dm15_from_s = float(dm15_s_interp(s_salt))

            dm15_from_x1 = np.nan
            if x1_salt2 != None:
                dm15_from_x1 = float(dm15_x1_interp(x1_salt2))

            dm15_from_delta = np.nan
            if delta_mlcs31 != None:
                dm15_from_delta = float(dm15_delta_interp(delta_mlcs31))

            dm15_from_delta_lowrv = np.nan
            if delta_lowrv != None:
                dm15_from_delta_lowrv = float(dm15_delta_lowrv_interp(delta_lowrv[0]))

            # dm15_from_fits = np.nanmean([dm15_from_s, dm15_from_x1, dm15_from_delta])
            # dm15s = [dm15_from_delta_lowrv, dm15_from_delta, dm15_from_s, dm15_from_x1] #change order to give fits different priority
            dm15s = [dm15_from_delta, dm15_from_x1, dm15_from_s, dm15_from_delta_lowrv] #change order to give fits different priority
            e_dm15s = [0.070, 0.066, 0.072, 0.074] #errors from fits
            for i, dm in enumerate(dm15s):
                if ~np.isnan(dm):
                    dm15_from_fits = dm
                    e_dm15 = e_dm15s[i]
                    break
                else:
                    dm15_from_fits = np.nan


            if np.isnan(dm15_from_fits):
                dm15_from_fits = None
                e_dm15 = None
        else:
            dm15_from_fits = None

        if event[0:2] == 'sn':
            event = event[2:]

        #photomotry blobs
        sn_phot = phot.get_photometry(event)
        phot_blob = msg.packb(sn_phot)

        csp_sn_phot = phot.get_csp_photometry(event)
        csp_phot_blob = msg.packb(csp_sn_phot)
        con.execute("""INSERT INTO Events(SN, RA, DEC, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
                                                       zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
                                                       zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
                                                       zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
                                                       glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
                                                       zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
                                                       av_MW, av_25, MJD_max, dm15_source, dm15_from_fits, e_dm15, separation, NED_host, v_at_max, v_err, 
                                                       Redshift, M_b_cfa, M_b_cfa_err, B_minus_V_cfa, B_minus_V_cfa_err, Carbon_presence, Na_presence, Hubble_res,
                                                       Photometry, csp_Photometry)
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (event, ra, dec, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
                                  zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
                                  zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
                                  zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
                                  glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
                                  zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
                                  av_mw, av_25, mjd_max, dm15_source, dm15_from_fits, e_dm15, sep, ned_host, vel, e_vel,
                                  redshift, m_b_cfa, m_b_cfa_err, b_minus_v_cfa, b_minus_v_cfa_err, carbon, na, hubble_res,
                                  buffer(phot_blob), buffer(csp_phot_blob))
                    )
    print 'Done'

    con.commit()

if __name__ == "__main__":
    main()


