import matplotlib.pyplot as plt
import numpy as np
import query_db
import composite
from tabulate import tabulate

def write_table(filename, table):
	with open('../../Paper_Drafts/' + filename, 'w') as file:
		file.write('\\begin{table*}[t]\n')
		file.write('\centering\n')
		file.write(table)
		file.write('\n')
		file.write('\caption{CAPTION}\n')
		file.write('\label{tab:1}\n')
		file.write('\end{table*}')
		file.close()

if __name__ == "__main__":
	# SN_Array = composite.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between -.2 and .2", multi_epoch = True)
	# tab_arr = []
	# for SN in SN_Array:
	# 	wav_range = str(int(np.round(SN.minwave, 0))) + ' - ' + str(int(np.round(SN.maxwave, 0)))
	# 	ref = '...'
	# 	tab_arr.append([SN.name, SN.source, SN.mjd, SN.phase, wav_range, ref])
	# table_1 = tabulate(tab_arr, headers=['SN Name', 'Source', 'mjd', 'Phase', 'Wavelength Range', 'Reference'], tablefmt = 'latex')
	# write_table('table_1.tex', table_1)

	SN_Array = composite.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between -.2 and .2", multi_epoch = True)
	tab_arr = []
	for SN in SN_Array:
		if SN.dm15_cfa != None:
			dm15 = str(np.round(SN.dm15_cfa, 2))
		elif SN.dm15_from_fits != None:
			dm15 = str(np.round(SN.dm15_from_fits, 2)) + '*'
		else:
			dm15 = None

		if SN.av_25 != None:
			av = SN.av_25
		elif SN.av_mlcs31 != None:
			av = SN.av_mlcs31
		elif SN.av_mlcs17 != None:
			av = SN.av_mlcs17
		else:
			av = None

		if SN.SNR != None:
			SN.SNR = np.round(SN.SNR, 2)
		tab_arr.append([SN.name, SN.source, SN.redshift, SN.SNR, dm15, SN.m_b, SN.B_minus_V, SN.velocity, SN.morph, SN.carbon, SN.GasRich, SN.resid, ])
	table_2 = tabulate(tab_arr, headers=['SN Name', 'Source', 'Redshift', 'SNR', '$\Delta m_15 (B)$', 'M_{B}', 'B - V (mag)', 
										 'Velocity (km/s)', 'Host Morphology', 'Carbon Presence', 'Gas Rich', 'Hubble Residual', 'A_{V}'], tablefmt = 'latex')
	write_table('table_2.tex', table_2)

#Photometric metadata:
	# (SN, RA, DEC, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
	#    zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
	#    zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
	#    zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
	#    glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
	#    zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
	#    av_25, dm15_cfa, dm15_from_fits, separation,
	#    Photometry)

#Spectral Metadata:
	# (Filename, SN, Source,
	# 	Redshift, Phase, MinWave, MaxWave, Dm15, M_B,
	# 	B_mMinusV_m, Velocity, Morphology, Carbon,
	# 	GasRich, snr, Hubble_Res, Interpolated_Spectra, MJD)