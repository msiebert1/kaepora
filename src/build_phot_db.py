from astropy.io import ascii
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn
import photometry as phot
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

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
		events.append('sn' + SN['SN'].lower())

	for SN in av_dict.keys():
		events.append(SN.lower())

	for SN in cfa_dict.keys():
		events.append(SN.lower())

	other_events = ['sn1989a', 'sn1989m', 'sn1993y', 'sn1993z', 'sn1999dg', 'sn1999do', 'sn2000dr', 'sn2000dx', 'sn2001dl', 'sn2001dt', 'sn2001dw', 'sn2002cv', 'sn2002eb', 
					'sn2002eh', 'sn2002el', 'sn2003gs', 'sn2003gt', 'sn2003he', 'sn2003hs', 'sn2004bl', 'sn2004br', 'sn2004bv', 'sn2004bw', 'sn2004e', 'sn2004go', 'sn2005de', 
					'sn2005di', 'sn2005dm', 'sn2005er', 'sn2005gj', 'sn2006dm', 'sn2007aj', 'sn2007cs', 'sn2007fr', 'sn2007ge', 'sn2007gi', 'sn2007gk', 'sn2007s1', 'sn2008bt', 
					'sn2008cl', 'sn2008dt', 'sn2008dx', 'sn2008ec', 'sn2008ei', 'sn2008ha', 'sn2008hs', 'sn2008s1', 'sn2008s5', 'sn2005A', 'sn2005M', 'sn2005W', 'sn2006dd', 
					'sn2006D', 'sn2006X', 'sn2007A', 'sn2007N', 'sn2007S', 'sn2008C', 'sn2008gl', 'sn2008hu', 'sn2008R', 'sn2009aa', 'sn2009ab', 'sn2009ad', 'sn2009ag', 
					'sn2009D', 'sn2009F', 'sn2009Y','sn1994e', 'sn1994b', 'sn1998dj', 'sn1994j', 'sn2002ec', 'sn2005lt', 'snsnf20080522-011', 'sn1994u', 'sn2008bz', 'sn1998dw', 
					'sn2002ep', 'sn1994x', 'sn2008fr', 'sn2004w', 'sn1997fc', 'sn2007e', 'sn2006es', 'sn2008fg', 'sn2008fj', 'sn2006eb', 'sn1995a', 'sn1995c', 'sn1999gf', 'sn2006ch', 
					'sn2003x', 'sn1998en', 'sn1990g', 'snsnf20080522-000', 'sn1995t', 'sn1993aj', 'sn1993ai', 'sn2004gl', 'sn2005dh', 'sn2002dr', 'sn1993ab', 'sn2005do', 'sn2007v', 
					'sn2006dy', 'sn2001fg', 'sn2006dw', 'sn2006dh', 'sn2006di', 'sn2006do', 'sn1995l', 'sn2008s4', 'sn2007m', 'sn1996o', 'sn2001a', 'sn2003v', 'sn2005bu', 'sn2002gx', 
					'sn2002gg', 'sn2002gf', 'sn2002gc', 'sn2002gb', 'snsnf20080720-001', 'sn1996p', 'sn2005as', 'sn2004gw', 'sn2008dr', 'sn2008ds', 'sn2001cg', 'sn2000ej', 'sn1992ap', 
					'sn2008db', 'sn2002fi', 'snsnf20080623-001', 'sn1998cl', 'sn1998cm', 'sn2006nr', 'sn1998cd', 'sn1997t', 'sn1999aq', 'sn1997fb', 'sn2004fw', 'sn1994ab', 'sn2005x', 
					'sn2008gw', 'sn2007fq', 'sn2005p', 'sn2000dd', 'sn2004fg', 'sn2000df', 'sn2008ge', 'sn2008gh', 'sn2005f', 'sn2004y', 'sn2004bq', 'sn1991bj', 'sn1990m', 'sn1991bd', 
					'sn1991bf', 'sn2003bf', 'sn2000j', 'sn1991bb', 'sn1991bc', 'sn2000q', 'sn1990r', 'snsnf20080514-002', 'sn2007ry', 'sn2001eg', 'sn1999c', 'sn2006ay', 'sn2008r3', 
					'sn2004eq', 'sn2001ew', 'sn2001eu', 'sn2001er', 'sn2008bv', 'sn2008bw', 'sn2002dx', 'sn2003dw', 'sn1991k', 'sn2007sa', 'sn1991ak', 'sn1991b', 'sn1991am', 'sn1991at', 
					'sn1991ay', 'sn2008ez', 'sn2008ey', 'sn2008er', 'sn2003p', 'sn2003ls', 'sn2008s3', 'sn2008gy', 'sn2002hl', 'sn2004di', 'sn2003lb', 'sn2008ee', 'sn2001dn', 'sn2003ax', 
					'sn2003ay', 'sn1998fc', 'sn2003au', 'sn2003av', 'sn2003ij', 'sn2003ik', 'sn2005ec', 'sn2003an', 'sn2003ah', 'sn1999fz', 'sn1992m', 'sn2000al', 'sn2006ct', 'sn2004cv', 
					'sn2003hj', 'sn2006ce', 'sn2004cb', 'sn2008hq', 'sn2008hr', 'sn2007qd', 'sn2008hy', 'sn2008hz', 'sn2002jo', 'sn2002bp', 'sn1993c', 'sn2002bk', 'sn2003bi', 'sn2003bh', 
					'sn2008hj', 'sn2008hk', 'sn2004bz', 'sn2001dd', 'sn2001dm', 'sn2008cd', 'sn2008cf', 'sn2008ca', 'sn2008cb', 'sn2001ds', 'sn2008ct', 'ASASSN-14lp','sn2007sr','sn2008Q',
					'sn2009an','sn2009dc','sn2009ig','sn2009Y','sn2010ev','sn2011aa','sn2011ao','sn2011by','sn2011fe','sn2011iv','sn2012cg','sn2012dn','sn2012fr','sn2012ht','sn2013aa',
					'sn2013cg','sn2014J','sn2016ccz','sn2016coj','sn2016eiy','sn2016ekg','sn2016eoa','sn2016fff','sn2016gsb','sn2016itd','sn2017cbv','sn2017erp']

	for event in other_events:
		events.append(event)

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
		
		host_dict['sn' + SN['SN'].lower()] = [ra,dec,glon,glat,cz,czLG,czCMB,mtype,xpos,ypos,t1,filt,Ebv]
		
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

def build_delt_dict(file):
     with open(file) as f:
        lines = f.readlines()
        delt_dict = {}
        for line in lines:
            l = line.split()    
            if len(l) == 30 and l[0] == 'SN:':
                delt_dict['sn' + l[1].lower()] = [float(l[16]), float(l[17])]
     return delt_dict

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

def build_swift_dict(data_file):
	with open(data_file) as data:
		lines = data.readlines()
		swift_dict = {}
		for line in lines:
			if not line.startswith('#'):
				sndata = line.split()
				# sndict[sndata[0]] = sndata[1:]
				swift_dict[sndata[0].lower()] = [sndata[1],sndata[3]]

	return swift_dict

def build_NED_host_dict(data_file):
	with open(data_file) as data:
		lines = data.readlines()
		NED_host_dict = {}
		for line in lines:
			if not line.startswith('#'):
				sndata = line.split()
				# sndict[sndata[0]] = sndata[1:]
				if sndata[1] != 'None':
					NED_host_dict['sn' + sndata[0].lower()] = sndata[1]
				else:
					NED_host_dict['sn' + sndata[0].lower()] = None
	return NED_host_dict

def dm15_from_fit_params(events, fit_dict, cfa_dict, stretch='N/A'):
	dm15_arr = []
	param_arr = []
	e_param_arr = []
	cfa_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	fit_none = [None,None,None,None,None,None,None,None,None,None]
	lowrv_none = [None,None]
	for event in events:
		if event != 'sn2001da' and event != 'sn2001cp': #outliers
			if event in fit_dict and event in cfa_dict:
				dm15_cfa = cfa_dict.get(event, cfa_none)[4]
				if dm15_cfa != '9.99' and dm15_cfa != None:
					dm15_cfa = float(dm15_cfa)
					dm15_arr.append(dm15_cfa)
					if stretch is 'delta_lowrv':
						param_arr.append(fit_dict.get(event, lowrv_none)[0])
						e_param_arr.append(fit_dict.get(event,lowrv_none)[1])
					else:
						param_arr.append(fit_dict.get(event, fit_none)[6]) #stretch param is 6th index of dicts
						e_param_arr.append(fit_dict.get(event,fit_none)[7])

	# coeffs = np.polyfit(param_arr, dm15_arr, 2)
	weights = 1./np.asarray(e_param_arr)
	coeffs = np.polyfit(dm15_arr, param_arr, 2, w=weights)
	# x = np.linspace(-5, 5.1, 10000)
	x = np.linspace(np.amin(dm15_arr), np.amax(dm15_arr), 10000)
	y = coeffs[0]*x**2. + coeffs[1]*x + coeffs[2]
	print y[0], y[-1]

	# plotted in lc_fitting_relationships.py
	plt.rc('font', family='serif')
	fig, ax = plt.subplots(1,1)
	fig.set_size_inches(10, 8, forward = True)
	plt.minorticks_on()
	plt.xticks(fontsize = 20)
	plt.yticks(fontsize = 20)
	plt.tick_params(
	    which='major', 
	    bottom='on', 
	    top='on',
	    left='on',
	    right='on',
	    length=10)
	plt.tick_params(
		which='minor', 
		bottom='on', 
		top='on',
		left='on',
		right='on',
		length=5)
	plt.plot(y, x, 'k', linewidth=4)
	plt.ylim(.6,2.2)
	plt.ylabel('$\Delta m_{15}$ (B)', fontsize = 30)
	if stretch is 's':
		plt.errorbar(param_arr, dm15_arr, xerr=e_param_arr, fmt='o', color='#7570b3', ms=10)
		plt.xlim(.4, 1.3)
		plt.xlabel('s', fontsize = 30)
		# plt.savefig('../../Paper_Drafts/dm15_s.pdf', dpi = 300, bbox_inches = 'tight')
	elif stretch is 'x1':
		plt.errorbar(param_arr, dm15_arr, xerr=e_param_arr, fmt='o', color='#1b9e77', ms=10)
		plt.xlim(-5., 3.)
		plt.xlabel('$x_1$', fontsize = 30)
		# plt.savefig('../../Paper_Drafts/dm15_x1.pdf', dpi = 300, bbox_inches = 'tight')
	elif stretch is 'delta':
		plt.errorbar(param_arr, dm15_arr, xerr=e_param_arr, fmt='o', color='#d95f02', ms=10)
		plt.xlim(-.5, 1.8)
		plt.xlabel('$\Delta$', fontsize = 30)
		# plt.savefig('../../Paper_Drafts/dm15_delta.pdf', dpi = 300, bbox_inches = 'tight')
	elif stretch is 'delta_lowrv':
		plt.errorbar(param_arr, dm15_arr, xerr=e_param_arr, fmt='o', color='#d95f02', ms=10)
		plt.xlim(-.9, 2.0)
		plt.xlabel('$\Delta$', fontsize = 30)
		# plt.savefig('../../Paper_Drafts/dm15_delta_lowrv.pdf', dpi = 300, bbox_inches = 'tight')
	plt.show()

	# dm15_interp = interp1d(x, y, bounds_error = True)
	dm15_interp = interp1d(y, x, bounds_error=False, fill_value=None)
	return dm15_interp


if __name__ == "__main__":
	mn.patch()

	salt = ascii.read("..\data\info_files\salt_params_dists.txt", delimiter = '\s', guess = False)
	salt2 = ascii.read("..\data\info_files\salt2_params_dists.txt", delimiter = '\s', guess = False)
	mlcs31 = ascii.read("..\data\info_files\mlcs31_params.txt", delimiter = '\s', guess = False)
	mlcs17 = ascii.read("..\data\info_files\mlcs17_params.txt", delimiter = '\s', guess = False)
	lcparams = ascii.read("..\data\info_files\lc_params.txt", delimiter = '\s', guess = False)
	host_data = ascii.read("..\data\info_files\other_host_data.txt", delimiter = '\s', guess = False)

	av_dict = build_av_dict('..\data\info_files\lowz_rv25_all.fitres')
	delt_dict = build_delt_dict('..\data\info_files\lowz_rv25_all.fitres')
	cfa_dict = build_cfa_dict('..\data\spectra\cfa\cfasnIa_param.dat')
	swift_dict = build_swift_dict('..\..\swift_uvspec\swift_uv_log.txt')
	NED_host_dict = build_NED_host_dict('..\data\info_files\NED_host_info_simplified.txt')

	events = find_all_events(salt,salt2,mlcs31,mlcs17,lcparams,host_data,av_dict,cfa_dict)

	salt_dict = build_salt_dict(salt)
	salt2_dict = build_salt2_dict(salt2)
	mlcs31_dict = build_mlcs31_dict(mlcs31)
	mlcs17_dict = build_mlcs17_dict(mlcs17)
	host_dict = build_host_dict(host_data)
	lcparams_dict = build_lcparams_dict(lcparams)

	dm15_s_interp = dm15_from_fit_params(events, salt_dict, cfa_dict, stretch='s')
	dm15_x1_interp = dm15_from_fit_params(events, salt2_dict, cfa_dict, stretch='x1')
	dm15_delta_interp = dm15_from_fit_params(events, mlcs31_dict, cfa_dict, stretch='delta')
	dm15_delta_lowrv_interp = dm15_from_fit_params(events, delt_dict, cfa_dict, stretch='delta_lowrv')

	con = sq3.connect('..\data\SNe_17_phot_3.db')
	con.execute("""DROP TABLE IF EXISTS Photometry""")
	con.execute("""CREATE TABLE IF NOT EXISTS Photometry (SN TEXT, RA TEXT, DEC TEXT, 
														  zCMB_salt REAL, e_zCMB_salt REAL, Bmag_salt REAL, e_Bmag_salt REAL, s_salt REAL, e_s_salt REAL, c_salt REAL, e_c_salt REAL, mu_salt REAL, e_mu_salt REAL,
														  zCMB_salt2 REAL, e_zCMB_salt2 REAL, Bmag_salt2 REAL, e_Bmag_salt2 REAL, x1_salt2 REAL, e_x1_salt2 REAL, c_salt2 REAL, e_c_salt2 REAL, mu_salt2 REAL, e_mu_salt2 REAL,
														  zCMB_mlcs31 REAL, e_zCMB_mlcs31 REAL, mu_mlcs31 REAL, e_mu_mlcs31 REAL, delta_mlcs31 REAL, e_delta_mlcs31 REAL, av_mlcs31 REAL, e_av_mlcs31 REAL,
														  zCMB_mlcs17 REAL, e_zCMB_mlcs17 REAL, mu_mlcs17 REAL, e_mu_mlcs17 REAL, delta_mlcs17 REAL, e_delta_mlcs17 REAL, av_mlcs17 REAL, e_av_mlcs17 REAL,
														  glon_host REAL, glat_host REAL, cz_host REAL, czLG_host REAL, czCMB_host REAL, mtype_host TEXT, xpos_host REAL, ypos_host REAL, t1_host REAL, filt_host TEXT, Ebv_host REAL,
														  zCMB_lc REAL, zhel_lc REAL, mb_lc REAL, e_mb_lc REAL, c_lc REAL, e_c_lc REAL, x1_lc REAL, e_x1_lc REAL, logMst_lc REAL, e_logMst_lc REAL, tmax_lc REAL, e_tmax_lc REAL, cov_mb_s_lc REAL, cov_mb_c_lc REAL, cov_s_c_lc REAL, bias_lc REAL,
														  av_25 REAL, dm15_source Real, dm15_from_fits REAL, separation REAL, NED_host REAL,
														  Photometry BLOB, csp_photometry BLOB)""")

	salt_none = [None,None,None,None,None,None,None,None,None,None,None,None]
	salt2_none = salt_none
	mlcs31_none = [None,None,None,None,None,None,None,None,None,None]
	mlcs17_none = mlcs31_none
	host_none = [None,None,None,None,None,None,None,None,None,None,None,None,None]
	lc_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	cfa_none = [None,None,None,None,None,None,None,None,None,None,None,None,None,None]
	swift_none = [None,None]
	
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

		sep = None
		if xpos_host != None and ypos_host != None:
			sep = (float(xpos_host)**2 + float(ypos_host)**2)**(.5)

		ned_host = NED_host_dict.get(event, None)

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
		delta_lowrv = delt_dict.get(event)

		#dm15 estimation
		dm15_source = cfa_dict.get(event, cfa_none)[4]
		if dm15_source != '9.99' and dm15_source != None:
			dm15_source = float(dm15_source)
		elif dm15_source == '9.99':
			dm15_source = None

		if dm15_source is None:
			dm15_source = swift_dict.get(event, swift_none)[1]
			if dm15_source != '-99' and dm15_source != None:
				dm15_source = float(dm15_source)
			elif dm15_source == '-99':
				dm15_source = None

		if dm15_source is None:
			dm15_from_s = np.NaN
			if s_salt != None:
				dm15_from_s = float(dm15_s_interp(s_salt))

			dm15_from_x1 = np.NaN
			if x1_salt2 != None:
				dm15_from_x1 = float(dm15_x1_interp(x1_salt2))

			dm15_from_delta = np.NaN
			if delta_mlcs31 != None:
				dm15_from_delta = float(dm15_delta_interp(delta_mlcs31))

			dm15_from_delta_lowrv = np.NaN
			if delta_lowrv != None:
				dm15_from_delta_lowrv = float(dm15_delta_lowrv_interp(delta_lowrv[0]))

			# dm15_from_fits = np.nanmean([dm15_from_s, dm15_from_x1, dm15_from_delta])
			dm15s = [dm15_from_delta_lowrv, dm15_from_delta, dm15_from_s, dm15_from_x1] #change order to give fits different priority
			for dm in dm15s:
				if dm != np.NaN:
					dm15_from_fits = dm
					break
				else:
					dm15_from_fits = np.NaN


			if dm15_from_fits == np.NaN:
				dm15_from_fits = None
		else:
			dm15_from_fits = None

		if event[0:2] == 'sn':
			event = event[2:]

		sn_phot = phot.get_photometry(event)
		phot_blob = msg.packb(sn_phot)

		csp_sn_phot = phot.get_csp_photometry(event)
		csp_phot_blob = msg.packb(csp_sn_phot)

		con.execute("""INSERT INTO Photometry(SN, RA, DEC, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
											           zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
											           zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
											           zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
											           glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
											           zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
											           av_25, dm15_source, dm15_from_fits, separation, NED_host,
											           Photometry, csp_Photometry)
                                VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (event, ra, dec, zCMB_salt, e_zCMB_salt, Bmag_salt, e_Bmag_salt, s_salt, e_s_salt, c_salt, e_c_salt, mu_salt, e_mu_salt,
                        	      zCMB_salt2, e_zCMB_salt2, Bmag_salt2, e_Bmag_salt2, x1_salt2, e_x1_salt2, c_salt2, e_c_salt2, mu_salt2, e_mu_salt2,
                        	      zCMB_mlcs31, e_zCMB_mlcs31, mu_mlcs31, e_mu_mlcs31, delta_mlcs31, e_delta_mlcs31, av_mlcs31, e_av_mlcs31,
                        	      zCMB_mlcs17, e_zCMB_mlcs17, mu_mlcs17, e_mu_mlcs17, delta_mlcs17, e_delta_mlcs17, av_mlcs17, e_av_mlcs17,
                        	      glon_host, glat_host, cz_host, czLG_host, czCMB_host, mtype_host, xpos_host, ypos_host, t1_host, filt_host, Ebv_host,
                        	      zCMB_lc, zhel_lc, mb_lc, e_mb_lc, c_lc, e_c_lc, x1_lc, e_x1_lc, logMst_lc, e_logMst_lc, tmax_lc, e_tmax_lc, cov_mb_s_lc, cov_mb_c_lc, cov_s_c_lc, bias_lc,
                        	      av_25, dm15_source, dm15_from_fits, sep, ned_host,
                        	      buffer(phot_blob), buffer(csp_phot_blob))
                    )
		print 'Done'

	con.commit()




