from __future__ import division
import numpy as np
import pandas as pd
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import prep_osc
import os
import re
import math
import time
import json
from pprint import pprint
import glob
import matplotlib.pyplot as plt

mn.patch()
global c
c = 299792.458

if __name__ == '__main__':
	con = sq3.connect('..\..\osc_data\osc_SNe.db')
	#make sure no prior table in db to avoid doubling/multiple copies of same data
	con.execute("""DROP TABLE IF EXISTS Supernovae""")
	con.execute("""CREATE TABLE IF NOT EXISTS Supernovae (Filename
	                    TEXT PRIMARY KEY, SN Text, Source Text, Redshift REAL,
	                    Phase REAL, MinWave REAL, MaxWave REAL, Dm15 REAL,
	                    M_B REAL, B_mMinusV_m REAL, Velocity REAL,
	                    Morphology INTEGER, Carbon TEXT, GasRich INTEGER, snr REAL,
	                    Interpolated_Spectra BLOB)""")


	files = glob.glob("..\..\osc_data\*.json")
	spectra_key = "spectra"

	for file in files:
		with open (file) as f:
			data = json.load(f)
			filename = file
			file_key = file.split('\\')[-1][:-5]
			data = data.get(file_key)

			if spectra_key in data:
				SN = data["name"]
				source = None
				phase = None
				dm15 = None
				carbon = None
				gasrich = None
				snr = None
				host = data.get("host", [{}])[0].get("value", None)

				redshift = float(data.get("redshift", [{}])[0].get("value", None))
				sn_type = data.get("claimed_type", [{}])[0].get("value", None)
				ebv = float(data.get("ebv", [{}])[0].get("value", None))
				max_absmag = float(data.get("maxabsmag", [{}])[0].get("value", None))
				max_appmag = float(data.get("maxappmag", [{}])[0].get("value", None))
				max_vis_absmag = float(data.get("maxvisualabsmag", [{}])[0].get("value", None))
				max_vis_appmag = float(data.get("maxvisualappmag", [{}])[0].get("value", None))
				velocity = float(data.get("velocity", [{}])[0].get("value", None))
				ra = str(data.get("ra", [{}])[0].get("value", None))
				dec = str(data.get("dec", [{}])[0].get("value", None))

				spec_num = 0
				for spectrum in data[spectra_key]:
					dereddened = "dereddened" in data[spectra_key][spec_num]
					deredshifted = "deredshifted" in data[spectra_key][spec_num]
					if len(data[spectra_key][spec_num]["data"][0]) == 2:
						u_fluxes = data[spectra_key][spec_num]["u_fluxes"]
						waves = np.asarray([float(x[0]) for x in data[spectra_key][spec_num]["data"]])
						fluxes = np.asarray([float(x[1]) for x in data[spectra_key][spec_num]["data"]])
						variances = np.zeros(len(waves))
						minwave = waves[0]
						maxwave = waves[-1]
						interp_spec, sig_noise = prep_osc.compprep(waves, fluxes, variances, redshift, ebv, dereddened, deredshifted, u_fluxes)
						print type(interp_spec)
						#deredden and interpolate
					else:
						u_fluxes = data[spectra_key][spec_num]["u_fluxes"]
						waves = np.asarray([float(x[0]) for x in data[spectra_key][spec_num]["data"]])
						fluxes = np.asarray([float(x[1]) for x in data[spectra_key][spec_num]["data"]])
						variances = np.asarray([float(x[2]) for x in data[spectra_key][spec_num]["data"]])
						minwave = waves[0]
						maxwave = waves[-1]
						interp_spec, sig_noise = prep_osc.compprep(waves, fluxes, variances, redshift, ebv, dereddened, deredshifted, u_fluxes)
						print type(interp_spec)
						#deredden and interpolate
					spec_num += 1
			else:
				pass