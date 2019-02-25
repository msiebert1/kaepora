import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import glob

mn.patch()
obj_file = '../../../maguire_hst/maguire_obj_info.txt'
spec_file = '../../../maguire_hst/maguire_spec_info.txt'
files = glob.glob("../../../maguire_hst/*.asci.txt")
# con = sq3.connect('../../../data/SNe_19_phot_9.db')
# PTF10fps = SN2010cr, PTF10hmv = SN2010dm, PTF10mwb = SN2010gn, PTF10qyx = SN2010gy, PTF10tce = SN 2010ho, PTF10ufj = SN 2010hs, PTF10ygu = SN 2010jn
with open(obj_file) as obj:
	with open(spec_file) as spec:
		obj_dict = {}
		spec_dict = {}
		for line in obj.readlines():
			obj_dict[line.split()[0]] = line.split()[1:]
		for line in spec.readlines():
			spec_dict[line.split()[0]] = line.split()[1:]
		for sp in files:
			with open(sp) as spec:
				spectrum = np.loadtxt(spec)
				source = 'maguire'
				filename = sp.split('\\')[-1]
				sn_obj_data = obj_dict[filename.split('_')[0]]
				sn_spec_data = spec_dict[filename.split('_')[0]]
				if filename.split('_')[0].lower() == 'sn':
					sn_name = filename.split('_')[0][2:]
				else:
					sn_name = filename.split('_')[0]
				print sn_name
				redshift = float(sn_spec_data[4])
				mjd_max = float(sn_obj_data[0])
				mjd_spec = float(sn_spec_data[3])
				phase = (mjd_spec - mjd_max)/(1.+redshift)

				interp_spec, sig_noise = prep.compprep(spectrum, sn_name, redshift, source)
				interped = msg.packb(interp_spec)
				name = spec_file.split('\\')[3]
				min_wave = spectrum[0][0]
				max_wave = spectrum[-1][0]

				m_b = None
				bm_vm = None
				vel = None
				morph = None
				carbon = None
				gasrich = None
				hubble_residual = None
				mjd = None
				ref = 'Maguire et al. (2012)'

				print sn_name, redshift, mjd_max, mjd_spec, phase, min_wave, max_wave, sig_noise
#				con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
							# 	Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
							# 	B_mMinusV_m, Velocity, Morphology, Carbon,
							# 	GasRich, snr, Hubble_Res, Interpolated_Spectra, 
							# 	MJD, Ref)
							# 	VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
							# (name, sn_name, source, redshift, phase,
							#  min_wave, max_wave, Dm15, m_b, bm_vm, vel,
							#  morph, carbon, gasrich, sig_noise, hubble_residual, 
							#  buffer(interped), mjd, ref))
# con.commit()
