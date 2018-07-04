import matplotlib.pyplot as plt
import numpy as np
import query_db
import composite
from tabulate import tabulate

def write_table(filename, table):
	with open('../../../Paper_Drafts/' + filename, 'w') as file:
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
	# 	if SN.SNR != None:
	# 		SN.SNR = np.round(SN.SNR, 2)
	# 	ref = '...'
	# 	if SN.ref is not None:
	# 		ref = SN.ref
	# 	else:
	# 		ref = '...'
	# 	tab_arr.append([SN.name, SN.source, SN.mjd, SN.phase, SN.SNR, wav_range, ref])
	# table_1 = tabulate(tab_arr, headers=['SN Name', 'Source', 'mjd', 'Phase', 'SNR', 'Wavelength Range', 'Reference'], tablefmt = 'latex')
	# write_table('table_1.tex', table_1)

	# SN_Array = composite.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between -.2 and .2", multi_epoch = True)
	# tab_arr = []
	# for SN in SN_Array:
	# 	if SN.dm15_cfa != None:
	# 		dm15 = str(np.round(SN.dm15_cfa, 2))
	# 	elif SN.dm15_from_fits != None:
	# 		dm15 = str(np.round(SN.dm15_from_fits, 2)) + '*'
	# 	else:
	# 		dm15 = None

	# 	if SN.av_25 != None:
	# 		av = SN.av_25
	# 	elif SN.av_mlcs31 != None:
	# 		av = SN.av_mlcs31
	# 	elif SN.av_mlcs17 != None:
	# 		av = SN.av_mlcs17
	# 	else:
	# 		av = None

	# 	tab_arr.append([SN.name, SN.source, SN.redshift, dm15, SN.m_b, SN.B_minus_V, SN.velocity, SN.morph, SN.resid, av])
	# table_2 = tabulate(tab_arr, headers=['SN Name', 'Source', 'Redshift', '$\Delta m_{15} (B)$', 'M_{B}', 'B - V (mag)', 
	# 									 'Velocity (km/s)', 'Host Morphology', 'Hubble Residual', 'A_{V}'], tablefmt = 'latex')
	# write_table('table_.tex', table_2)

	tab_arr = []
	SN_Array = composite.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN", multi_epoch = False, make_corr = False)
	num_tmax = 0
	num_red = 0
	num_dm15 = 0
	num_dm15fit = 0
	num_delta = 0
	num_s = 0
	num_x1 = 0
	num_av = 0
	num_host = 0
	num_hostmass = 0
	num_mb = 0
	num_bmv = 0
	num_vel = 0
	num_carbon = 0
	num_gas = 0
	num_hubres=0

	for SN in SN_Array:
		if SN.phase != None:
			num_tmax+=1
		if SN.redshift != None:
			num_red+=1
		if SN.dm15_source != None:
			num_dm15+=1
		if SN.dm15_source is None and SN.dm15_from_fits != None:
			num_dm15fit+=1
		if SN.av_25 != None:
			num_av+=1
		if SN.ned_host != None:
			num_host+=1
		if SN.m_b != None:
			num_mb+=1
		if SN.B_minus_V != None:
			num_bmv+=1
		if SN.v_at_max != None:
			num_vel+=1
		if SN.carbon != None:
			num_carbon+=1
		if SN.GasRich != None:
			num_gas+=1
		if SN.resid != None:
			num_hubres+=1
	# tab_arr.append(['Redshift', num_red])
	# tab_arr.append(['$\Delta m_{15} (B)$', num_dm15])
	# tab_arr.append(['$\Delta m_{15} (B)$ from fit parameter', num_dm15fit])
	# # tab_arr.append(['$\Delta$', num_delta])
	# # tab_arr.append(['s', num_s])
	# # tab_arr.append(['x1', num_x1])
	# tab_arr.append(['$A_V$', num_av])
	# tab_arr.append(['Host Morphology', num_host])
	# tab_arr.append(['Host Stellar Mass', num_hostmass])
	# tab_arr.append(['$M_B$', num_mb])
	# tab_arr.append(['$B-V$', num_bmv])
	# tab_arr.append(['Velocity at Max $m_B$', num_vel])
	# tab_arr.append(['Carbon Presence', num_carbon])
	# tab_arr.append(['Gas Rich', num_gas])
	# tab_arr.append(['Hubble Residual', num_hubres])

	print 'tmax', num_tmax
	print 'Redshift', num_red
	print '$\Delta m_{15} (B)$', num_dm15
	print '$\Delta m_{15} (B)$ from fit parameter', num_dm15fit
	print '$A_V$', num_av
	print 'Host Morphology', num_host
	print 'Host Stellar Mass', num_hostmass
	print '$M_B$', num_mb
	print '$B-V$', num_bmv
	print 'Velocity at Max $m_B$', num_vel
	print 'Carbon Presence', num_carbon
	print 'Gas Rich', num_gas
	print 'Hubble Residual', num_hubres

	# table_3 = tabulate(tab_arr, headers=['Property', 'Number of SNe'], tablefmt = 'latex')
	# write_table('table_3.tex', table_3)

#Photometric metadata:
# (event, ra, dec, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
#                         	      zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
#                         	      zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
#                         	      zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
#                         	      glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
#                         	      zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
#                         	      av_25, dm15_source, dm15_from_fits, e_dm15, sep, ned_host, vel, e_vel,
#                         	      buffer(phot_blob), buffer(csp_phot_blob))

#Spectral Metadata:
# (name, sn_name, source, redshift, phase,
#                          min_wave, max_wave, Dm15, m_b, bm_vm, vel,
#                          morph, carbon, gasrich, sig_noise, hubble_residual, 
#                          buffer(interped), mjd, ref))