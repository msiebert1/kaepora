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
    cur.execute('PRAGMA TABLE_INFO({})'.format("Events"))

    names = [tup[1] for tup in cur.fetchall()]
    if col_name not in names:
        cur.execute("""ALTER TABLE Events ADD COLUMN """+ col_name + """ """ + dtype)
        for SN in data.keys():
            cur.execute("UPDATE Events SET " + col_name + " = ? where SN = ?", (data[SN], SN))
    else:
        for SN in data.keys():
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
        for SN in data.keys():
            cur.execute("UPDATE Spectra SET " + col_name + " = ? where filename = ?", (data[filename], filename))
    con.commit()


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

    alpha = 0.14107
    beta = 3.14889
    gamma = 0.04572
    for line in salt2_data:
        SN                      = line['CID'].lower()
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

    col_names = ['MU', 'MUERR', 'MUMODEL', 'MURES', 'HOST_LOGMASS', 'HOST_LOGMASS_ERR', 'x1', 'x1ERR', 
                 'c', 'cERR', 'MURES_NO_MSTEP', 'MURES_NO_MSTEP_C', 'MURES_NO_MSTEP_C_X1']
    all_dicts = [mu, mu_err, mumodel, mures, hostmass, hostmass_err, x1, x1_err, c, 
                 c_err, mures_no_mstep, mures_no_mstep_c, mures_no_mstep_c_x1]
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
    add_all_SALT2_metadata('../data/kaepora_v1.db')

    





