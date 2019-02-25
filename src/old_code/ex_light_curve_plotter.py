import find_event_data as fed
import composite as comp
import matplotlib.pyplot as plt

def plot_light_curves(light_curves):
	p = light_curves
	times = []
	mags = []
	bands = []
	if p != None:
		for band in p:
			t = []
			m = []
			for i in range(len(p[band][0])):
				t.append(float(p[band][0][i]))
				m.append(float(p[band][1][i][0]))
				# print p[band][1][i][1] #this is the dictionary with additional fields
				# print p[band][1][i][1].get('telescope') #will print telescope field if it exists


			times.append(t)
			mags.append(m)
			bands.append(band)

		for i in range(len(bands)):
			plt.plot(times[i], mags[i], 'o', label = bands[i])
			# if bands[i] == 'I':
			# 	plt.plot(times[i], mags[i], 'o', label = bands[i])
		plt.gca().invert_yaxis()
		plt.legend()
		plt.show()
	else:
		print 'Insufficient Data'

PH_Array = fed.grab_phot_data("SELECT * FROM Photometry")
print len(PH_Array)
print PH_Array[1].name
plot_light_curves(PH_Array[1].light_curves)

SN_Array = comp.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where velocity between -98 and -12", make_corr = False)
print len(SN_Array)
print SN_Array[1].name
plot_light_curves(SN_Array[1].light_curves)

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

# print PH_Array[0].light_curves

