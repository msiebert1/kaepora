import sqlite3 as sq3
import numpy as np 
import datafidelity as df 
import prep as prep
import msgpack as msg
import composite
import matplotlib.pyplot as plt

def fix_2011fe_phases():
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where Supernovae.SN = '2011fe' and source = 'other'"
	print 'querying'
	cur.execute(sql_input)
	for row in cur.fetchall():
		filename  = row[0]
		with open('../data/spectra/other/' + filename) as file:
			lines = file.readlines()
			for line in lines:
				if len(line.split()) > 1 and line.split()[1] == 'TMAX':
					phase = float(line.split()[3])
					print filename, float(line.split()[3])
					cur.execute("UPDATE Supernovae SET phase = ? where filename = ?", (phase, filename))
					break
	con.commit()

def update_bsnip_refs():
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where source = 'bsnip'"
	print 'querying'
	cur.execute(sql_input)
	file = '../data/info_files/bsnip_references.txt'
	with open(file) as f:
		lines = f.readlines()
		for row in cur.fetchall():
			filename  = row[0]
			ref       = row[18]
			for line in lines:
				if filename in line:
					print filename, ref, line.split()[-1]
					new_ref = line.split()[-1]
					cur.execute("UPDATE Supernovae SET Ref = ? where filename = ?", (new_ref, filename))
	con.commit()

def delete_swift_data():
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	cur.execute("DELETE FROM Supernovae where source = 'swift_uv'")
	con.commit()

def build_peak_dict(file):
     with open(file) as f:
        lines = f.readlines()

        peak_mjd_dict = {}
        for line in lines:
            l = line.split()
            if len(l) == 30 and l[0] == 'SN:':
                peak_mjd_dict[l[1].lower()] = float(l[14])
     return peak_mjd_dict

def read_cfa_info(data_file, dates_file):
    """
    Reads Cfa SN information from separate files.
    Output dict format:
    # key,     0         1           2       3        4        5       6     7        8        9
    # SN,  zhel,  tmax(B),  +/-  ref.,  Dm15, +/-  ref.,  M_B   +/-,   B-V,
     10        11       12       13              14
    +/-,  Bm-Vm, +/-,   Phot. ref   Obs Date
    """
    with open(dates_file) as dates:
        lines = dates.readlines()
        date_dict = {}
        for line in lines:
            if not line.startswith('#'):
                cleandate = line.split()

                if cleandate:
                    # date_dict[cleandate[0]] = cleandate[1]
                    date_dict[cleandate[0]] = cleandate[1]

    with open(data_file) as data:
        lines = data.readlines()
        sndict = {}
        for line in lines:
            if not line.startswith('#'):
                sndata = line.split()
                # sndict[sndata[0]] = sndata[1:]
                sndict[sndata[0].lower()] = sndata[1:]

    return sndict, date_dict

def update_phases_from_mlcs():
	#updates phases using time of max from David's mlcs fits
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where Supernovae.SN"
	print 'querying'
	peak_mjd_dict = build_peak_dict('..\data\info_files\lowz_rv25_all.fitres')
	sndict, date_dict = read_cfa_info('../data/spectra/cfa/cfasnIa_param.dat',
                                  '../data/spectra/cfa/cfasnIa_mjdspec.dat')
	cur.execute(sql_input)
	for row in cur.fetchall():
		filename  = row[0]
		name      = row[1]
		source    = row[2]
		redshift  = row[3]
		phase     = row[4]
		mjd         = row[17]
		# if mjd is None and source == 'cfa':
		# 	print name, source, float(date_dict[filename])
		# 	mjd = float(date_dict[filename])
		# 	cur.execute("UPDATE Supernovae SET mjd = ? where filename = ?", (mjd, filename))
		if name in peak_mjd_dict and mjd is not None and redshift is not None:
			print name, filename, source, mjd, peak_mjd_dict[name], phase
			phase = (mjd - peak_mjd_dict[name])/(1.+redshift)
			cur.execute("UPDATE Supernovae SET phase = ? where filename = ?", (phase, filename))
	con.commit()
def fix_snr_measurements():
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN"
	print 'querying'
	bad_ivars = []
	SN_Array = composite.grab(sql_input, multi_epoch = True, make_corr = False, selection = 'max_coverage', grab_all=True)
	for SN in SN_Array:
		nan_bool_flux = np.isnan(SN.flux)
		non_nan_data = np.where(nan_bool_flux == False)
		non_nan_data = np.array(non_nan_data[0])
		if len(non_nan_data) > 0:
			x1 = non_nan_data[0]
			x2 = non_nan_data[-1]
			error = 1./np.sqrt(SN.ivar[x1:x2])
			snr = np.nanmedian(SN.flux[x1:x2] / error)
			print SN.filename, SN.SNR, snr
			if not np.isnan(snr):
				cur.execute("UPDATE Supernovae SET snr = ? where filename = ?", (snr, SN.filename))
	con.commit()
	print 'done'

def repair_bad_variance_spectra():
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	sql_input = "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN"
	print 'querying'
	bad_ivars = []
	SN_Array = composite.grab(sql_input, multi_epoch = True, make_corr = False, selection = 'max_coverage', grab_all=True)
	for SN in SN_Array:
		if len(np.where(np.isnan(SN.ivar))[0] == True) == 5500:
			bad_ivars.append(SN)

	for SN in bad_ivars:
		print SN.source
		if SN.source == 'bsnip':
			path = '../data/spectra/bsnip/' + SN.filename
		elif SN.source == 'cfa':
			path = '../data/spectra/cfa/' + SN.filename.split('-')[0] + '/' + SN.filename
		elif SN.source == 'other':
			path = '../data/spectra/other/' + SN.filename
		elif SN.source == 'uv':
			path = '../data/spectra/uv/' + SN.filename
		spectrum = np.loadtxt(path)
		newdata, snr = prep.compprep(spectrum, SN.name, SN.redshift, SN.source, use_old_error=False, testing = True)
		try:
			interped = msg.packb(newdata)
			cur.execute("UPDATE Supernovae SET snr = ? where filename = ?", (snr.value, SN.filename))
			cur.execute("UPDATE Supernovae SET Interpolated_Spectra = ? where filename = ?", (buffer(interped), SN.filename))
			print "Added: ", SN.filename
		except Exception, e:
			print "Interp failed: ", SN.filename
		# raise TypeError
		# old_wave = spectrum[:, 0]
		# old_flux = spectrum[:, 1]
		# real_error = spectrum[:, 2]
		# old_error = np.zeros(len(old_wave), float)
		# new_ivar = df.genivar(old_wave, old_flux, old_error)
		# new_var = 1./new_ivar
		# real_var = real_error*real_error
		# new_var = new_var*2.02
		# newdata = Interpo(new_wave, new_flux, new_ivar)
	con.commit()
	print 'done'

def add_more_host_data():
	datafile = '../data/info_files/more_host_info.txt'
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	with open(datafile) as data:
		lines = data.readlines()
		for line in lines:
			sndata = line.split()
			print sndata[0].lower(), sndata[3]
			cur.execute("UPDATE Photometry SET NED_host = ? where SN = ?", (sndata[3], sndata[0].lower()))
	con.commit()

def update_references():
	datafile = '../data/info_files/more_references.txt'
	con = sq3.connect('../data/SNIaDB_Spec_v20_phot_v10.db')
	cur = con.cursor()
	with open(datafile) as data:
		lines = data.readlines()
		for line in lines:
			filename, ref = line.split()
			if ref == "Unknown":
				ref = None
			print filename, ref
			cur.execute("UPDATE Supernovae SET Ref = ? where filename = ?", (ref, filename))
	con.commit()

if __name__ == "__main__":
	# fix_2011fe_phases()
	# delete_swift_data()
	# update_bsnip_refs()
	# repair_bad_variance_spectra()
	# update_phases_from_mlcs()
	# fix_snr_measurements()
	# add_more_host_data()
	update_references()