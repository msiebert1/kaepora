from astropy.io import ascii
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import numpy as np
import composite

def build_salt2_ID_dict(salt2):
    salt2_ID_dict = {}
    for SN in salt2:
        id_num = str(SN['IDSURVEY'])
        salt2_ID_dict[SN['CID'].lower()] = id_num
        
    return salt2_ID_dict

def add_arbitrary_event_column(col_name, data, dtype, db_file):
    con = sq3.connect(db_file)
    cur = con.cursor()
    cur.execute("PRAGMA TABLE_INFO({})".format("Events"))
    names = [tup[1] for tup in cur.fetchall()]

    cur.execute("SELECT SN from Events")
    sn_names = [tup[0] for tup in cur.fetchall()]

    if col_name not in names:
        cur.execute("""ALTER TABLE Events ADD COLUMN """+ col_name + """ """ + dtype)
        for SN in data.keys():
            if SN not in sn_names:
                cur.execute("INSERT INTO Events (SN, "+ col_name + ") VALUES (?,?)", (SN, data[SN]))
            else:
                cur.execute("UPDATE Events SET " + col_name + " = ? where SN = ?", (data[SN], SN))
    else:
        for SN in data.keys():
            if SN not in sn_names:
                cur.execute("INSERT INTO Events (SN, "+ col_name + ") VALUES (?,?)", (SN, data[SN]))
            else:
                print (SN, data[SN], col_name)
                cur.execute("UPDATE Events SET " + col_name + " = ? where SN = ?", (data[SN], SN))
    con.commit()

def add_arbitrary_spectra_column(col_name, data, dtype, db_file):
    con = sq3.connect(db_file)
    cur = con.cursor()
    cur.execute('PRAGMA TABLE_INFO({})'.format("Spectra"))

    names = [tup[1] for tup in cur.fetchall()]
    if col_name not in names:
        cur.execute("""ALTER TABLE Spectra ADD COLUMN """+ col_name + """ """ + dtype)
        for filename in data.keys():
            cur.execute("UPDATE Spectra SET " + col_name + " = ? where filename = ?", (data[filename], filename))
    else:
        for filename in data.keys():
            cur.execute("UPDATE Spectra SET " + col_name + " = ? where filename = ?", (data[filename], filename))
    con.commit()

def add_homogenized_photometry(db_file):
    csp_phot = '../data/info_files/CSP_phot.txt'
    phot_dict = {}
    with open(csp_phot) as csp:
        lines = csp.readlines()
        for line in lines[1:]:
            name = line.split()[0].lower()[2:]
            date = float(line.split()[1])
            band = line.split()[2]
            mag = float(line.split()[3])
            mag_err = float(line.split()[4])
            if name not in phot_dict.keys():
                phot_dict[name] = {}
                phot_dict[name][band] = [[],[],[],'CSP']
                phot_dict[name][band][0].append(date)
                phot_dict[name][band][1].append(mag)
                phot_dict[name][band][2].append(mag_err)
            else:
                if band in phot_dict[name].keys():
                    phot_dict[name][band][0].append(date)
                    phot_dict[name][band][1].append(mag)
                    phot_dict[name][band][2].append(mag_err)
                else:
                    phot_dict[name][band] = [[],[],[],'CSP']
                    phot_dict[name][band][0].append(date)
                    phot_dict[name][band][1].append(mag)
                    phot_dict[name][band][2].append(mag_err)
        for name in phot_dict.keys():
            for band in phot_dict[name].keys():
                mjd_order = np.argsort(np.asarray(phot_dict[name][band][0]))
                phot_dict[name][band][0] = [phot_dict[name][band][0][i] for i in mjd_order]
                phot_dict[name][band][1] = [phot_dict[name][band][1][i] for i in mjd_order]
                phot_dict[name][band][2] = [phot_dict[name][band][2][i] for i in mjd_order]

    con = sq3.connect(db_file)
    cur = con.cursor()
    cur.execute('PRAGMA TABLE_INFO({})'.format("Events"))

    names = [tup[1] for tup in cur.fetchall()]
    # print names
    if 'Homogenized_Photometry' not in names:
        cur.execute("""ALTER TABLE Events ADD COLUMN Homogenized_Photometry TEXT""")

        sql_input = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN"
        # print 'querying'
        SN_Array = composite.grab(sql_input, multi_epoch = False, make_corr = False, selection = 'max_coverage', grab_all=True)
        for SN in SN_Array:
            if SN.name in phot_dict.keys():
                phot_blob = msg.packb(phot_dict[SN.name])
                cur.execute("UPDATE Events SET Homogenized_Photometry = ? where SN = ?", (buffer(phot_blob), SN.name))
            else:
                cur.execute("UPDATE Events SET Homogenized_Photometry = ? where SN = ?", (None, SN.name))
    else:
        sql_input = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN"
        # print 'querying'
        SN_Array = composite.grab(sql_input, multi_epoch = False, make_corr = False, selection = 'max_coverage', grab_all=True)
        for SN in SN_Array:
            if SN.name in phot_dict.keys():
                phot_blob = msg.packb(phot_dict[SN.name])
                cur.execute("UPDATE Events SET Homogenized_Photometry = ? where SN = ?", (buffer(phot_blob), SN.name))
            else:
                cur.execute("UPDATE Events SET Homogenized_Photometry = ? where SN = ?", (None, SN.name))
    con.commit()

def add_flux_cal_scaling(db_file):
    filenames, scales = np.genfromtxt('../data/info_files/flux_cal_data.txt', unpack=True, dtype='string')
    scale_dict = {}
    for i, file in enumerate(filenames):
        scale_dict[file] = float(scales[i])
    add_arbitrary_spectra_column('Flux_Cal_Scale', scale_dict, """REAL""", db_file)

def add_galaxy_metadata(db_file):
    galaxy_data = ascii.read("../data/info_files/localsn_public_cuts.txt", delimiter = r'\s', guess = False)

    local_umg = {}
    global_umg = {}
    local_mass = {}
    global_mass = {}
    localssfr = {}
    globalssfr = {}

    for line in galaxy_data[1:]:
        SN                      = line['ID'].lower()
        if line['local_u-g'] != '-99.0':
            local_umg[SN]           = float(line['local_u-g'])
        if line['global_u-g'] != '-99.0':
            global_umg[SN]          = float(line['global_u-g'])
        if line['localmass'] != '-99.0':
            local_mass[SN]          = float(line['localmass'])
        if line['globalmass'] != '-99.0':
            global_mass[SN]         = float(line['globalmass'])
        if line['localssfr'] != '-99.0':
            localssfr[SN]           = float(line['localssfr'])
        if line['globalssfr'] != '-99.0':
            globalssfr[SN]          = float(line['globalssfr'])

    col_names = ['local_umg', 'global_umg', 'localmass', 'globalmass', 'localssfr', 'globalssfr']
    all_dicts = [local_umg, global_umg, local_mass, global_mass, localssfr, globalssfr]

    for i, sdict in enumerate(all_dicts):
        add_arbitrary_event_column(col_names[i], sdict, """REAL""", db_file)

def add_more_galaxy_metadata(db_file):
    names, masses, sfrs, ssfrs, umgs = np.genfromtxt('../data/info_files/hubblephot_sedfits_global_brief.txt', unpack=True, dtype='string')
    cepheid_data = np.genfromtxt('../data/info_files/cepheidphot_lephare_sedfits.txt', unpack=True, dtype='string')
    ceph_names = cepheid_data[0]
    ceph_ssfrs = cepheid_data[8]

    global_mass = {}
    global_ssfr = {}

    con = sq3.connect(db_file)
    cur = con.cursor()
    spec_table = cur.execute("SELECT SN, globalssfr from EVENTS")
    data = [tup for tup in cur.fetchall()]
    ignore_names = []
    for d in data:
        if d[1] != None:
            ignore_names.append(str(d[0]))

    for i in range(len(names)):
        if ssfrs[i] != '-99.000' and sfrs[i] != '-99.000':
            if names[i][0:2] == 'SN':
                if names[i].lower()[2:] not in ignore_names:
                    global_ssfr[names[i].lower()[2:]] = float(ssfrs[i])
            else:
                if names[i].lower() not in ignore_names:
                    global_ssfr[names[i].lower()] = float(ssfrs[i])

    for i in range(len(ceph_names)):
        if ssfrs[i] != '-99.99':
            if names[i][0:2] == 'SN':
                if ceph_names[i].lower()[2:] not in ignore_names:
                    global_ssfr[ceph_names[i].lower()[2:]] = float(ceph_ssfrs[i])
            else:
                if ceph_names[i].lower() not in ignore_names:
                    global_ssfr[ceph_names[i].lower()] = float(ceph_ssfrs[i])

    add_arbitrary_event_column('globalssfr', global_ssfr, """REAL""", db_file)

def add_all_SALT2_metadata(db_file):
    salt2_data = ascii.read("../data/info_files/SALT2mu_fpan.fitres", delimiter = r'\s', guess = False)

    mu = {}
    mu_err = {}
    mumodel = {}
    mures = {}
    hostmass = {}
    hostmass_err = {}
    x1 = {}
    x1_err = {}
    c = {}
    c_err = {}

    mures_no_mstep = {}
    mures_no_mstep_c = {}
    mures_no_mstep_c_x1 = {}
    mures_no_mstep_x1 = {}

    alpha = 0.14107
    beta = 3.14889
    gamma = 0.04572
    for line in salt2_data:
        SN                      = line['CID'].lower()
        if SN.startswith('at') and not SN.startswith('atlas'):
            SN = SN[2:]
        if '-' in SN:
            SN = SN.replace('-', '')
            # print SN
        mu[SN]                  = float(line['MU'])
        mu_err[SN]              = float(line['MUERR'])
        mumodel[SN]             = float(line['MUMODEL'])
        mures[SN]               = float(line['MURES'])
        hostmass[SN]            = float(line['HOST_LOGMASS'])
        hostmass_err[SN]        = float(line['HOST_LOGMASS_ERR'])
        x1[SN]                  = float(line['x1'])
        x1_err[SN]              = float(line['x1ERR'])
        c[SN]                   = float(line['c'])
        c_err[SN]               = float(line['cERR'])

        if hostmass[SN] >= 10:
            mures_no_mstep[SN] = mures[SN] - gamma
        else:
            mures_no_mstep[SN] = mures[SN] + gamma
        
        mures_no_mstep_c[SN] = mures_no_mstep[SN] + beta*c[SN]
        mures_no_mstep_c_x1[SN] = mures_no_mstep_c[SN] - alpha*x1[SN]
        mures_no_mstep_x1[SN] = mures_no_mstep[SN] - alpha*x1[SN]

    col_names = ['MU', 'MUERR', 'MUMODEL', 'MURES', 'HOST_LOGMASS', 'HOST_LOGMASS_ERR', 'x1', 'x1ERR', 
                 'c', 'cERR', 'MURES_NO_MSTEP', 'MURES_NO_MSTEP_C', 'MURES_NO_MSTEP_C_X1', 'MURES_NO_MSTEP_X1']
    all_dicts = [mu, mu_err, mumodel, mures, hostmass, hostmass_err, x1, x1_err, c, 
                 c_err, mures_no_mstep, mures_no_mstep_c, mures_no_mstep_c_x1, mures_no_mstep_x1]
    for i, sdict in enumerate(all_dicts):
        add_arbitrary_event_column(col_names[i], sdict, """REAL""", db_file)

def add_salt2_survey_ID_column(db_file):
    salt2 = ascii.read("../data/info_files/SALT2mu_fpan.fitres", delimiter = r'\s', guess = False)
    salt2_ID_dict = build_salt2_ID_dict(salt2)

    id_file = "../data/info_files/SURVEY.DEF"
    id_dict = {}
    with open(id_file) as file:
        id_lines = file.readlines()
        for line in id_lines:
            if not line.startswith('#') and len(line.split()) > 2 and line.split()[0] == 'SURVEY:':
                id_dict[line.split()[2]] = line.split()[1]
        sn_id_dict = {}
        for sn in salt2_ID_dict.keys():
            sn_id_dict[sn] = id_dict[salt2_ID_dict[sn]]

    con = sq3.connect(db_file)
    cur = con.cursor()
    cur.execute('PRAGMA TABLE_INFO({})'.format("Events"))

    names = [tup[1] for tup in cur.fetchall()]
    # print names
    if 'salt2_phot_source' not in names:
        cur.execute("""ALTER TABLE Events ADD COLUMN salt2_phot_source TEXT""")

        sql_input = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN"
        # print 'querying'
        SN_Array = composite.grab(sql_input, multi_epoch = False, make_corr = False, selection = 'max_coverage', grab_all=True)
        for SN in SN_Array:
            if SN.name in sn_id_dict.keys():
                # print SN.name, sn_id_dict[SN.name]
                cur.execute("UPDATE Events SET salt2_phot_source = ? where SN = ?", (sn_id_dict[SN.name], SN.name))
            else:
                cur.execute("UPDATE Events SET salt2_phot_source = ? where SN = ?", (None, SN.name))
    else:
        sql_input = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN"
        # print 'querying'
        SN_Array = composite.grab(sql_input, multi_epoch = False, make_corr = False, selection = 'max_coverage', grab_all=True)
        for SN in SN_Array:
            if SN.name in sn_id_dict.keys():
                # print SN.name, sn_id_dict[SN.name]
                cur.execute("UPDATE Events SET salt2_phot_source = ? where SN = ?", (sn_id_dict[SN.name], SN.name))
            else:
                cur.execute("UPDATE Events SET salt2_phot_source = ? where SN = ?", (None, SN.name))
    con.commit()

if __name__ == "__main__":
    # add_salt2_survey_ID_column('../data/kaepora_v1.db')
    add_all_SALT2_metadata('../data/kaepora_v1_DEV.db')
    # add_homogenized_photometry('../data/kaepora_v1.db')
    # add_galaxy_metadata('../data/kaepora_v1_DEV.db')
    # add_more_galaxy_metadata('../data/kaepora_v1_DEV.db')

    





