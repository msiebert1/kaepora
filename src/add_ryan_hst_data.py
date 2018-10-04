import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import glob
import matplotlib.pyplot as plt

mn.patch()
data_file_meta = '../../../ryan_hst/hst_metadata.txt'
data_file = '../../../ryan_hst/hst_metadata_spec.txt'
files = glob.glob("../../../ryan_hst/*.flm")
con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
cur = con.cursor()
# cur.execute("DELETE FROM Supernovae where source = 'foley_hst'")
# with open(data_file) as data:
# 	data_dict = {}
# 	for line in data.readlines():
# 		data_dict[line.split()[0]] = line.split()
# 	for spec_file in files:
# 		with open(spec_file) as spec:
# 			spectrum = np.loadtxt(spec)
# 			source = 'foley_hst'
# 			sn_data = data_dict[spec_file.split('\\')[1]]
# 			if sn_data[0][0:2].lower() == 'sn':
# 				sn_name = sn_data[0].split('-')[0][2:]
# 			else:
# 				sn_name = sn_data[0].split('-')[0] + '-' + sn_data[0].split('-')[1]
			
# 			redshift = float(sn_data[2])
# 			phase = float(sn_data[1])
# 			Dm15 = None

# 			interp_spec, sig_noise = prep.compprep(spectrum, sn_name, redshift, source, testing =False)
# 			interped = msg.packb(interp_spec)
# 			name = spec_file.split('\\')[1]
# 			min_wave = spectrum[0][0]
# 			max_wave = spectrum[-1][0]
# 			sig_noise = sig_noise
# 			print name, sn_name, phase, redshift, sig_noise
# 			# plt.plot(interp_spec[0,:],interp_spec[1,:])
# 			# plt.show()
# 			# plt.plot(interp_spec[0,:],interp_spec[2,:])
# 			# plt.show()

# 			m_b = None
# 			bm_vm = None
# 			vel = None
# 			morph = None
# 			carbon = None
# 			gasrich = None
# 			hubble_residual = None
# 			mjd = None
# 			ref = sn_data[3]

# 			sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where filename=" + "'" + name + "'"
# 			cur.execute(sql_input)
# 			in_database = False
# 			for row in cur:
# 				in_database = row[0]
# 			if in_database:
# 				con.execute("""REPLACE INTO Supernovae(Filename, SN, Source,
# 								Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
# 								B_mMinusV_m, Velocity, Morphology, Carbon,
# 								GasRich, snr, Hubble_Res, Interpolated_Spectra, 
# 								MJD, Ref)
# 								VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
# 							(name, sn_name, source, redshift, phase,
# 							 min_wave, max_wave, Dm15, m_b, bm_vm, vel,
# 							 morph, carbon, gasrich, sig_noise, hubble_residual, 
# 							 buffer(interped), mjd, ref))
# 			else:
# 				con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
# 								Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
# 								B_mMinusV_m, Velocity, Morphology, Carbon,
# 								GasRich, snr, Hubble_Res, Interpolated_Spectra, 
# 								MJD, Ref)
# 								VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
# 							(name, sn_name, source, redshift, phase,
# 							 min_wave, max_wave, Dm15, m_b, bm_vm, vel,
# 							 morph, carbon, gasrich, sig_noise, hubble_residual, 
# 							 buffer(interped), mjd, ref))


with open(data_file_meta) as data_meta:
	data_dict = {}
	for line in data_meta.readlines():
		data_dict[line.split()[0]] = line.split()
	for sn in data_dict.keys():
		dm15 = float(data_dict[sn][1])
		dm15_err = float(data_dict[sn][2])
		v_at_max = float(data_dict[sn][3])
		av = float(data_dict[sn][5])*3.1
		ned_host = data_dict[sn][7]
		print sn, dm15, dm15_err, v_at_max, av, ned_host
		cur.execute("UPDATE Photometry SET dm15_source = ? where SN = ?", (dm15, sn))
		cur.execute("UPDATE Photometry SET e_dm15 = ? where SN = ?", (dm15_err, sn))
		cur.execute("UPDATE Photometry SET v_at_max = ? where SN = ?", (v_at_max, sn))
		cur.execute("UPDATE Photometry SET av_mlcs31 = ? where SN = ?", (av, sn))
		cur.execute("UPDATE Photometry SET NED_host = ? where SN = ?", (ned_host, sn))
con.commit()
