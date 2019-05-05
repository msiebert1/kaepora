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
	add_salt2_survey_ID_column('../data/kaepora_v1.db')

	





