from astropy.io import ascii
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import photometry as phot

def find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,cfa_dict):
	events = []

	for SN in salt:
		events.append(SN['SimbadName'].lower())

	for SN in salt2:
		events.append(SN['SimbadName'].lower())

	for SN in mlcs31:
		events.append(SN['SimbadName'].lower())

	for SN in mlcs17:
		events.append(SN['SimbadName'].lower())

	for SN in lcparams:
		events.append(SN['SimbadName'].lower())

	for SN in host_data:
		events.append(SN['SN'].lower())

	for SN in av_dict.keys():
		events.append(SN.lower())

	for SN in cfa_dict.keys():
		events.append(SN.lower())

	events = set(events)
	return events


def build_salt_dict(salt):
	salt_dict = {}
	for SN in salt:
		ra = SN['RA']
		dec = SN['DEC']
		zCMB = float(SN['zCMB'])
		e_zCMB = float(SN['e_zCMB'])
		Bmag = float(SN['Bmag'])
		e_Bmag = float(SN['e_Bmag'])
		s = float(SN['s'])
		e_s = float(SN['e_s'])
		c = float(SN['c'])	
		e_c = float(SN['e_c'])
		mu = float(SN['mu'])
		e_mu = float(SN['e_mu'])

		salt_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,Bmag,e_Bmag,s,e_s,c,e_c,mu,e_mu]

	return salt_dict

def build_salt2_dict(salt2):
	salt2_dict = {}
	for SN in salt2:
		ra = SN['RA']
		dec = SN['DEC']
		zCMB = float(SN['zCMB'])
		e_zCMB = float(SN['e_zCMB'])
		Bmag = float(SN['Bmag'])
		e_Bmag = float(SN['e_Bmag'])
		c = float(SN['c'])	
		e_c = float(SN['e_c'])
		mu = float(SN['mu'])
		e_mu = float(SN['e_mu'])
		x1 = float(SN['x1'])
		e_x1 = float(SN['e_x1'])

		salt2_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,Bmag,e_Bmag,x1,e_x1,c,e_c,mu,e_mu]
		
	return salt2_dict

def build_mlcs31_dict(mlcs31):
	mlcs31_dict = {}
	for SN in mlcs31:
		ra = SN['RA']
		dec = SN['DEC']
		zCMB = float(SN['zCMB'])
		e_zCMB = float(SN['e_zCMB'])
		mu = float(SN['mu'])
		e_mu = float(SN['e_mu'])
		delta = float(SN['Delta'])
		e_delta = float(SN['e_Delt'])
		av = float(SN['AV'])
		e_av = float(SN['e_AV'])

		mlcs31_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,mu,e_mu,delta,e_delta,av,e_av]
		
	return mlcs31_dict

def build_mlcs17_dict(mlcs17):
	mlcs17_dict = {}
	for SN in mlcs17:
		ra = SN['RA']
		dec = SN['DEC']
		zCMB = float(SN['zCMB'])
		e_zCMB = float(SN['e_zCMB'])
		mu = float(SN['mu'])
		e_mu = float(SN['e_mu'])
		delta = float(SN['Delta'])
		e_delta = float(SN['e_Delt'])
		av = float(SN['AV'])
		e_av = float(SN['e_AV'])

		mlcs17_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,e_zCMB,mu,e_mu,delta,e_delta,av,e_av]
		
	return mlcs17_dict

def build_host_dict(host_data):
	host_dict = {}
	for SN in host_data:
		ra = SN['RA']
		dec = SN['DEC']
		glon = float(SN['GLON'])
		glat = float(SN['GLAT'])
		cz = float(SN['cz'])
		czLG = float(SN['czLG'])
		czCMB = float(SN['czCMB'])
		mtype = SN['MType']
		xpos = float(SN['Xpos'])
		ypos = float(SN['Ypos'])
		t1 = float(SN['t1'])
		filt = SN['Filt']
		Ebv = float(SN['E(B-V)'])
		
		host_dict[SN['SN'].lower()] = [ra,dec,glon,glat,cz,czLG,czCMB,mtype,xpos,ypos,t1,filt,Ebv]
		
	return host_dict

def build_lcparams_dict(lcparams):
	lcparams_dict = {}
	for SN in lcparams:
		ra = SN['RA']
		dec = SN['DEC']
		zCMB = float(SN['zcmb'])
		zhel = float(SN['zhel'])
		mb = float(SN['mb'])
		e_mb = float(SN['e_mb'])
		c = float(SN['c'])	
		e_c = float(SN['e_c'])
		x1 = float(SN['x1'])
		e_x1 = float(SN['e_x1'])
		logMst = float(SN['logMst'])
		e_logMst = float(SN['e_log'])
		tmax = float(SN['tmax'])
		e_tmax = float(SN['e_tmax'])
		cov_mb_s = float(SN['cov_mb_s'])
		cov_mb_c = float(SN['cov_mb_c'])
		cov_s_c = float(SN['cov_s_c'])
		bias = float(SN['bias'])

		lcparams_dict[SN['SimbadName'].lower()] = [ra,dec,zCMB,zhel,mb,e_mb,c,e_c,x1,e_x1,logMst,e_logMst,tmax,e_tmax,cov_mb_s,cov_mb_c,cov_s_c,bias]
		
	return lcparams_dict

def build_av_dict(file):
     with open(file) as f:
        lines = f.readlines()

        av_dict = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                av_dict['sn' + l[1].lower()] = float(l[18])

     return av_dict

def build_cfa_dict(data_file):
	with open(data_file) as data:
		lines = data.readlines()
		cfa_dict = {}
		for line in lines:
			if not line.startswith('#'):
				sndata = line.split()
				# sndict[sndata[0]] = sndata[1:]
				cfa_dict['sn' + sndata[0].lower()] = sndata[1:]

	return cfa_dict

if __name__ == "__main__":
	mn.patch()

	salt = ascii.read("..\data\info_files\salt_params_dists.txt", delimiter = '\s', guess = False)
	salt2 = ascii.read("..\data\info_files\salt2_params_dists.txt", delimiter = '\s', guess = False)
	mlcs31 = ascii.read("..\data\info_files\mlcs31_params.txt", delimiter = '\s', guess = False)
	mlcs17 = ascii.read("..\data\info_files\mlcs17_params.txt", delimiter = '\s', guess = False)
	lcparams = ascii.read("..\data\info_files\lc_params.txt", delimiter = '\s', guess = False)
	host_data = ascii.read("..\data\info_files\other_host_data.txt", delimiter = '\s', guess = False)

	av_dict = build_av_dict('..\data\info_files\lowz_rv25_all.fitres')
	cfa_dict = build_cfa_dict('..\data\spectra\cfa\cfasnIa_param.dat')

	events = find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,cfa_dict)

	salt_dict = build_salt_dict(salt)
	salt2_dict = build_salt2_dict(salt2)
	mlcs31_dict = build_mlcs31_dict(mlcs31)
	mlcs17_dict = build_mlcs17_dict(mlcs17)
	host_dict = build_host_dict(host_data)
	lcparams_dict = build_lcparams_dict(lcparams)

	con = sq3.connect('..\data\SNe_14_phot_1.db')
	con.execute("""DROP TABLE IF EXISTS Photometry""")
	con.execute("""CREATE TABLE IF NOT EXISTS Photometry (SN TEXT, RA TEXT, DEC TEXT, 
														  zCMB_salt REAL, e_zCMB_salt REAL, Bmag_salt REAL, e_Bmag_salt REAL, s_salt REAL, e_s_salt REAL, c_salt REAL, e_c_salt REAL, mu_salt REAL, e_mu_salt REAL,
														  zCMB_salt2 REAL, e_zCMB_salt2 REAL, Bmag_salt2 REAL, e_Bmag_salt2 REAL, x1_salt2 REAL, e_x1_salt2 REAL, c_salt2 REAL, e_c_salt2 REAL, mu_salt2 REAL, e_mu_salt2 REAL,
														  zCMB_mlcs31 REAL, e_zCMB_mlcs31 REAL, mu_mlcs31 REAL, e_mu_mlcs31 REAL, delta_mlcs31 REAL, e_delta_mlcs31 REAL, av_mlcs31 REAL, e_av_mlcs31 REAL,
														  zCMB_mlcs17 REAL, e_zCMB_mlcs17 REAL, mu_mlcs17 REAL, e_mu_mlcs17 REAL, delta_mlcs17 REAL, e_delta_mlcs17 REAL, av_mlcs17 REAL, e_av_mlcs17 REAL,
														  glon_host REAL, glat_host REAL, cz_host REAL, czLG_host REAL, czCMB_host REAL, mtype_host TEXT, xpos_host REAL, ypos_host REAL, t1_host REAL, filt_host TEXT, Ebv_host REAL,
														  zCMB_lc REAL, zhel_lc REAL, mb_lc REAL, e_mb_lc REAL, c_lc REAL, e_c_lc REAL, x1_lc REAL, e_x1_lc REAL, logMst_lc REAL, e_logMst_lc REAL, tmax_lc REAL, e_tmax_lc REAL, cov_mb_s_lc REAL, cov_mb_c_lc REAL, cov_s_c_lc REAL, bias_lc REAL,
														  av_25 REAL, dm15_cfa Real,
														  Photometry BLOB)""")

	salt_none = [None,None,None,None,None,None,None,None,None,None,None,None]
	salt2_none = salt_none
	mlcs31_none = [None,None,None,None,None,None,None,None,None,None]
	mlcs17_none = mlcs31_none
	host_none = [None,None,None,None,None,None,None,None,None,None,None,None,None]
	lc_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	cfa_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	for event in events:
		print 'Adding data for ' + event + '...'
		ra_salt,dec_salt,zCMB_salt,e_zCMB_salt,Bmag_salt,e_Bmag_salt,s_salt,e_s_salt,c_salt,e_c_salt,mu_salt,e_mu_salt = salt_dict.get(event, salt_none)
		ra_salt2,dec_salt2,zCMB_salt2,e_zCMB_salt2,Bmag_salt2,e_Bmag_salt2,x1_salt2,e_x1_salt2,c_salt2,e_c_salt2,mu_salt2,e_mu_salt2 = salt2_dict.get(event, salt2_none)
		ra_mlcs31,dec_mlcs31,zCMB_mlcs31,e_zCMB_mlcs31,mu_mlcs31,e_mu_mlcs31,delta_mlcs31,e_delta_mlcs31,av_mlcs31,e_av_mlcs31 = mlcs31_dict.get(event, mlcs31_none)
		ra_mlcs17,dec_mlcs17,zCMB_mlcs17,e_zCMB_mlcs17,mu_mlcs17,e_mu_mlcs17,delta_mlcs17,e_delta_mlcs17,av_mlcs17,e_av_mlcs17 = mlcs17_dict.get(event, mlcs17_none)
		ra_host,dec_host,glon_host,glat_host,cz_host,czLG_host,czCMB_host,mtype_host,xpos_host,ypos_host,t1_host,filt_host,Ebv_host = host_dict.get(event, host_none)
		ra_lc,dec_lc,zCMB_lc,zhel_lc,mb_lc,e_mb_lc,c_lc,e_c_lc,x1_lc,e_x1_lc,logMst_lc,e_logMst_lc,tmax_lc,e_tmax_lc,cov_mb_s_lc,cov_mb_c_lc,cov_s_c_lc,bias_lc = lcparams_dict.get(event, lc_none)
		ras = [ra_salt,ra_salt2,ra_mlcs31,ra_mlcs17,ra_host,ra_lc]
		decs = [dec_salt,dec_salt2,dec_mlcs31,dec_mlcs17,dec_host,dec_lc]

		ra = None
		for r in ras:
			if r != None:
				ra = r
				break

		dec = None
		for d in decs:
			if d != None:
				dec = d
				break

		av_25 = av_dict.get(event)
		dm15_cfa = cfa_dict.get(event, cfa_none)[4]
		if dm15_cfa != '9.99' and dm15_cfa != None:
			dm15_cfa = float(dm15_cfa)
		elif dm15_cfa == '9.99':
			dm15_cfa = None

		if event[0:2] == 'sn':
			event = event[2:]

		sn_phot = phot.get_photometry(event)
		phot_blob = msg.packb(sn_phot)


		con.execute("""INSERT INTO Photometry(SN, RA, DEC, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
											           zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
											           zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
											           zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
											           glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
											           zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
											           av_25, dm15_cfa,
											           Photometry)
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (event, ra, dec, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
                        	      zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
                        	      zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
                        	      zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
                        	      glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
                        	      zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
                        	      av_25, dm15_cfa,
                        	      buffer(phot_blob))
                    )
		print 'Done'

	con.commit()




