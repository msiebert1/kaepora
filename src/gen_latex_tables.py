import matplotlib.pyplot as plt
import numpy as np
import query_db
import composite
from tabulate import tabulate

def write_table(filename, table, caption = None):
	with open('../../../Paper_Drafts/' + filename, 'w') as file:
		file.write('\\begin{table*}[t]\n')
		file.write('\centering\n')
		file.write(table)
		file.write('\n')
		if caption != None:
			file.write('\caption{%s}\n'%caption)
		else:
			file.write('\caption{CAPTION}\n')
		file.write('\label{tab:1}\n')
		file.write('\end{table*}')
		file.close()

if __name__ == "__main__":
	# SN_Array = composite.grab("SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN order by Spectra.SN", multi_epoch = True, make_corr = False, grab_all=True)
	# tab_arr = []
	# refs = []
	# for SN in SN_Array:
	# 	if SN.ref is not None:
	# 		ref = SN.ref
	# 	else:
	# 		ref = 'Unknown'
	# 	refs.append(ref)

	# ref_set = sorted(set(refs), key=refs.index)

	# ref_nums = []
	# for i in range(len(ref_set)):
	# 	ref_nums.append(i+1)
	# for i, ref in enumerate(ref_set):
	# 	print ref, ref_nums[i]

	# i=1
	# for SN in SN_Array:
	# 	wav_range = str(int(np.round(SN.minwave, 0))) + ' - ' + str(int(np.round(SN.maxwave, 0)))
	# 	if SN.SNR != None:
	# 		SNR = SN.SNR
	# 	ref = 'Unknown'
	# 	if SN.ref is not None:
	# 		ref = SN.ref
	# 	else:
	# 		ref = 'Unknown'

	# 	if ref in ref_set:
	# 		ref_num = ref_nums[ref_set.index(ref)]
	# 	# ref_num = ref_dict[ref]

	# 	if SN.phase != None:
	# 		phase = '%.1f'%SN.phase
	# 	else:
	# 		phase = '...'

	# 	if SN.mjd != None:
	# 		mjd = '%.1f'%SN.mjd
	# 	else:
	# 		mjd = '...'

	# 	if len(SN.name) == 5:
	# 		name = SN.name.upper()
	# 	else:
	# 		name = SN.name

	# 	tab_arr.append([name, mjd, phase, SNR, wav_range, ref_num])
	# 	i+=1
	# 	if i > 22:
	# 		break

	# table_1 = tabulate(tab_arr, headers=['SN Name', 'MJD', 'Phase', 'S/N', 'Wavelength Range', 'Reference'], tablefmt = 'latex', floatfmt=".1f")

	# text = " (This table is available in its entirety in a machine-readable form in the online journal. A portion is shown here for guidance regarding its form and content.)"
	# ref_text = 'Full spectral sample. \\textbf{References}: '
	# for i, ref in enumerate(ref_set):
	# 	ref_text = ref_text + '(%d) \\citet{%s}'%(ref_nums[i], ref)
	# 	if i < len(ref_set)-1:
	# 		ref_text = ref_text + '; '

	# caption = ref_text + text
	# write_table('table_1_updated.tex', table_1, caption) #lccrcc
	# raise TypeError

	tab_arr = []
	SN_Array = composite.grab("SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN", multi_epoch = False, make_corr = False)
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
		else:
			print SN.name, SN.filename, SN.source
		if SN.dm15_source != None:
			num_dm15+=1
		if SN.dm15_source is None and SN.dm15_from_fits != None:
			num_dm15fit+=1
		if SN.av_25 != None:
			num_av+=1
		if SN.ned_host != None:
			num_host+=1
		if SN.m_b_cfa != None:
			num_mb+=1
		if SN.b_minus_v_cfa != None:
			num_bmv+=1
		if SN.v_at_max != None:
			num_vel+=1
		if SN.carbon != None:
			num_carbon+=1
		if SN.na != None:
			num_gas+=1
		if SN.hubble_res != None:
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

	print len(SN_Array), 'total events'
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