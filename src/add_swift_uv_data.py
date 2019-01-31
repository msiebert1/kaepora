import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep
import os
import glob

mn.patch()

def main():
	data_file = '../data/spectra/swift_uvspec/swift_uv_log.txt'
	files = glob.glob("..\data\spectra\swift_uvspec\*.flm")
	con = sq3.connect('../data/kaepora_v1.db')
	with open(data_file) as data:
		data_dict = {}
		for line in data.readlines()[1:]:
			data_dict[line.split()[4]] = line.split()[0:4]
		for spec_file in files:
			with open(spec_file) as spec:
				spectrum = np.loadtxt(spec)
				source = 'swift_uv'
				print spec_file	
				sn_data = data_dict[spec_file.split('\\')[4]]
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
				name = spec_file.split('\\')[4]
				min_wave = spectrum[0][0]
				max_wave = spectrum[-1][0]

				mjd = None
				ref = None
				con.execute("""INSERT INTO Spectra(Filename, SN, Source,
								Phase, MinWave, MaxWave, snr, Interpolated_Spectra, 
								MJD, Ref)
								VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
							(name, sn_name, source, phase, min_wave, max_wave, 
							 sig_noise, buffer(interped), mjd, ref))
	con.commit()
	
if __name__ == "__main__":
	main()