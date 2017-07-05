import build_phot_db as phot_db
import numpy as np
from astropy.io import ascii
import sqlite3 as sq3

# con = sq3.connect('../data/SNe_14_phot_1.db')
# cur = con.cursor()
# sql_input = "SELECT * FROM Photometry where SN = '1997bp'"
# cur.execute(sql_input)
# for row in cur:
# 	print row[67]
salt = ascii.read("..\data\info_files\salt_params_dists.txt", delimiter = '\s', guess = False)
salt2 = ascii.read("..\data\info_files\salt2_params_dists.txt", delimiter = '\s', guess = False)
mlcs31 = ascii.read("..\data\info_files\mlcs31_params.txt", delimiter = '\s', guess = False)
mlcs17 = ascii.read("..\data\info_files\mlcs17_params.txt", delimiter = '\s', guess = False)
lcparams = ascii.read("..\data\info_files\lc_params.txt", delimiter = '\s', guess = False)
host_data = ascii.read("..\data\info_files\other_host_data.txt", delimiter = '\s', guess = False)

av_dict = phot_db.build_av_dict('..\data\info_files\lowz_rv25_all.fitres')
cfa_dict = phot_db.build_cfa_dict('..\data\spectra\cfa\cfasnIa_param.dat')

events = phot_db.find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,cfa_dict)

salt_dict = phot_db.build_salt_dict(salt)
salt2_dict = phot_db.build_salt2_dict(salt2)
mlcs31_dict = phot_db.build_mlcs31_dict(mlcs31)
mlcs17_dict = phot_db.build_mlcs17_dict(mlcs17)
host_dict = phot_db.build_host_dict(host_data)
lcparams_dict = phot_db.build_lcparams_dict(lcparams)

dm15_s_interp = phot_db.dm15_from_fit_params(events, salt_dict, cfa_dict)
dm15_x1_interp = phot_db.dm15_from_fit_params(events, salt2_dict, cfa_dict)
dm15_delta_interp = phot_db.dm15_from_fit_params(events, mlcs31_dict, cfa_dict)


#uncomment plotting code in dm15_from_fit_params in phot_db to see relationships
for event in events:
	dm15_from_s = np.NaN
	dm15_from_x1 = np.NaN
	dm15_from_delta = np.NaN
	if salt_dict.get(event) != None:
		dm15_from_s = float(dm15_s_interp(salt_dict.get(event,np.NaN)[6]))
	if salt2_dict.get(event) != None:
		# print salt2_dict.get(event)[6]
		dm15_from_x1 = float(dm15_x1_interp(salt2_dict.get(event,np.NaN)[6]))
	if mlcs31_dict.get(event) != None:
		dm15_from_delta = float(dm15_delta_interp(mlcs31_dict.get(event,np.NaN)[6]))

	print event, dm15_from_s, dm15_from_x1, dm15_from_delta, np.nanmean([dm15_from_s, dm15_from_x1, dm15_from_delta])