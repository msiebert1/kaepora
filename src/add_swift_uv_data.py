import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import glob

mn.patch()
data_file = '../../swift_uvspec/swift_uv_log.txt'
files = glob.glob("..\..\swift_uvspec\*.flm")
con = sq3.connect('../data/SNe_17_phot_1.db')
with open(data_file) as data:
	data_dict = {}
	for line in data.readlines()[1:]:
		data_dict[line.split()[4]] = line.split()[0:4]
	for spec_file in files:
		with open(spec_file) as spec:
			spectrum = np.loadtxt(spec)
			source = 'swift_uv'
			print spec_file	
			sn_data = data_dict[spec_file.split('\\')[3]]
			if sn_data[0][0:2].lower() == 'sn':
				sn_name = sn_data[0][2:]
			else:
				sn_name = sn_name = sn_data[0]
			print sn_name
			redshift = float(sn_data[1])
			phase = float(sn_data[2])
			Dm15 = sn_data[3]
			if Dm15 == '-99':
				Dm15 = None
			else:
				Dm15 = float(Dm15)

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
			ref = None

			con.execute("""INSERT INTO Supernovae(Filename, SN, Source,
							Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
							B_mMinusV_m, Velocity, Morphology, Carbon,
							GasRich, snr, Hubble_Res, Interpolated_Spectra, 
							MJD, Ref)
							VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
						(name, sn_name, source, redshift, phase,
						 min_wave, max_wave, Dm15, m_b, bm_vm, vel,
						 morph, carbon, gasrich, sig_noise, hubble_residual, 
						 buffer(interped), mjd, ref))
con.commit()
