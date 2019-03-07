import composite
import numpy as np
import sys
import argparse
import spectral_analysis as sa
import kaepora_plot as kplot
import warnings
from tabulate import tabulate
warnings.filterwarnings("ignore")

def set_min_num_spec(composites, num):

	for comp in composites:
		comp.spec_bin = np.array(comp.spec_bin)
		valid_range = np.where(comp.spec_bin >= num)[0]
		comp.x1, comp.x2 = valid_range[0], valid_range[-1]


def normalize_comps(composites, scale=1., w1=3500, w2=9000.):
	for comp in composites:
		scale_range = np.where((comp.wavelength > w1) & (comp.wavelength < w2))
		# norm = scale/np.amax(comp.flux[comp.x1:comp.x2])
		norm = scale/np.nanmax(comp.flux[scale_range])
		comp.flux = norm*comp.flux

		# if len(comp.RMSE) > 0:
		# 	comp.RMSE = comp.RMSE*(norm)

		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			comp.low_conf = norm*comp.low_conf
			comp.up_conf = norm*comp.up_conf
		comp.ivar /= (norm)**2
	return composites


def normalize_comp(comp):
	norm = 1./np.nanmax(comp.flux[comp.x1:comp.x2])
	comp.flux = norm*comp.flux

	if len(comp.RMSE) > 0:
		comp.RMSE = comp.RMSE*(norm)

	if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
		comp.low_conf = norm*comp.low_conf	
		comp.up_conf = norm*comp.up_conf
	comp.ivar /= (norm)**2
	return comp, norm

def grab(query, multi_epoch = False, make_corr = False, grab_all=True):
	spec_array = composite.grab(query, multi_epoch = multi_epoch, make_corr = make_corr, grab_all=grab_all)
	return spec_array


def make_composite(num_queries, query_strings, boot=False, medmean = 1, selection = 'max_coverage', gini_balance=False, verbose=True, 
		 multi_epoch=True, combine=True, low_av_test=None, measure_vs = False, og_arr=False):
	# num_queries = int(sys.argv[1])
	# query_strings = sys.argv[2:]

	composites = []
	sn_arrays = []
	og_sn_arrays = []
	boot_sn_arrays = []
	store_boots = True
	for n in range(num_queries):
		if og_arr:
			comp, arr, og_arr, boots = composite.main(query_strings[n],boot=boot, medmean = medmean, 
											selection = selection, gini_balance=gini_balance, combine=combine,
											verbose=verbose, multi_epoch=multi_epoch, low_av_test=low_av_test, og_arr=og_arr)
			og_sn_arrays.append(og_arr)
		else:
			comp, arr, boots = composite.main(query_strings[n],boot=boot, medmean = medmean, 
											selection = selection, gini_balance=gini_balance, combine=combine,
											verbose=verbose, multi_epoch=multi_epoch, low_av_test=low_av_test, og_arr=og_arr)
		if store_boots:
			boot_sn_arrays.append(boots)
		composites.append(comp)
		sn_arrays.append(arr)

	composite.optimize_scales(composites, composites[0], True)
	composites = normalize_comps(composites)
	# if measure_vs:
	# 	for comp in composites:
	# 		dm15 = np.round(np.nanmean(comp.dm15_array[comp.x1:comp.x2]),2)
	# 		# r = sa.measure_si_ratio(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], vexp = .001, dm15=dm15)
	# 		v_strong, si_min_wave = sa.measure_velocity(comp.wavelength[comp.x1:comp.x2],comp.flux[comp.x1:comp.x2], 5900., 6300.)
	# 		print 'v = ', v_strong

	# 	for boot in boot_sn_arrays:
	# 		vs = []
	# 		for b in boot:
	# 			v_strong, si_min_wave = sa.measure_velocity(b.wavelength[b.x1:b.x2],b.flux[b.x1:b.x2], 5900., 6300.)
	# 			vs.append(v_strong)
	# 		v_err = np.nanstd(vs)
	# 		print 'v_err = ', v_err

	if og_arr:
		return composites, sn_arrays, og_sn_arrays, boot_sn_arrays
	else:
		return composites, sn_arrays, boot_sn_arrays


def save_comps_to_files(composites, prefix):
	#save to file
	for SN in composites:
		set_min_num_spec(composites, 5)
		phase = np.round(np.average(SN.phase_array[SN.x1:SN.x2]), 2)
		dm15 = np.round(np.average(SN.dm15_array[SN.x1:SN.x2]), 2)
		z = np.round(np.average(SN.red_array[SN.x1:SN.x2]), 3)
		num = np.amax(np.array(SN.spec_bin[SN.x1:SN.x2]))
		print phase, dm15, z
		set_min_num_spec(composites, 1)

		if phase >= 0.:
			sign = 'p'
		else:
			sign = 'm'
		abs_phase = np.absolute(phase)
		phase_str = str(abs_phase)
		dm15_str = str(dm15)
		z_str = str(z)
		num_str = str(SN.num_sne)
		num_spec_str = str(SN.num_spec)

		file_path = '../data/S19_Composite_Spectra/' + prefix + '_N=' + num_str + '_Nspec=' + num_spec_str + '_phase='+ sign + phase_str + '_dm15=' + dm15_str + '_z=' + z_str+'.txt'
		print file_path
		with open(file_path, 'w') as file:
			file.write('# SQL Query: ' + SN.query + '\n')
			wave = np.array(SN.wavelength[SN.x1:SN.x2])
			flux = np.array(SN.flux[SN.x1:SN.x2])
			up_conf = np.array(SN.up_conf[SN.x1:SN.x2]) - flux
			low_conf = flux - np.array(SN.low_conf[SN.x1:SN.x2])
			phase_arr = np.array(SN.phase_array[SN.x1:SN.x2])
			dm15_arr = np.array(SN.dm15_array[SN.x1:SN.x2])
			red_arr = np.array(SN.red_array[SN.x1:SN.x2])
			spec_per_sn_arr = np.array(SN.spec_bin[SN.x1:SN.x2])
			data = np.c_[wave,flux,low_conf,up_conf,phase_arr,dm15_arr, red_arr, spec_per_sn_arr]
			table = tabulate(data, headers=['Wavelength', 'Flux', '1-Sigma Lower', '1-Sigma Upper', 'Phase', 'Dm15', 'Redshift', 'SNe per Bin'], 
												tablefmt = 'ascii')
			file.write(table)

if __name__ == "__main__":
	composites = []
	SN_Arrays = []
	boot_sn_arrays = []
	store_boots = True

	boot = sys.argv[1]
	if boot == 'nb':
		boot = False
	else:
		boot = True

	query_strings = sys.argv[2:]

	num_queries = len(query_strings)

	for n in range(num_queries):
		c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = False, make_corr=False, multi_epoch=True, combine=False)
		# c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = False, multi_epoch=True, combine=True) 
		# composites.append(c)
		# SN_Arrays.append(sn_arr)
		# if store_boots:
		# 	boot_sn_arrays.append(boots)

		#use this for good composites
		# c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = True, verbose=True, multi_epoch=True, combine=True)
		# composites.append(c)
		# SN_Arrays.append(sn_arr)
		# if store_boots:
		# 	boot_sn_arrays.append(boots)
			
		# c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = True, verbose=True, multi_epoch=True, combine=True)
		composites.append(c)
		SN_Arrays.append(sn_arr)
		if store_boots:
			boot_sn_arrays.append(boots)

	for i, comp in enumerate(composites):
		dm15s = []
		kplot.plot_comp_and_all_spectra(comp, SN_Arrays[i], show_ivar=True)
		# for SN in SN_Arrays[i]:
		# 	if SN.dm15_source is not None:
		# 		dm15s.append(SN.dm15_source)
		# 	else:
		# 		dm15s.append(SN.dm15_from_fits)
		# plt.hist(dm15s)
		# plt.show()
	# comp2 = composites[1]
	# for sn in sn_arr:
	# 	r = sa.measure_si_ratio(sn.wavelength[sn.x1:sn.x2], sn.flux[sn.x1:sn.x2])
	# 	print r

	bad_vels = ['1995bd','2012cg','1991t','1998ab','2006oa','2009ig','2012fr','1995ac', '2001ah',
            '1997br','2008ia','2012ht','2007sr','2002do','2007ax','2003fa','1998aq','2011fe']
	# for SN in SN_Arrays[0]:
	# 	# if SN.name not in bad_vels and SN.source != 'swift_uv':
	# 	if SN.name and SN.source != 'swift_uv':
	# 		var = 1./SN.ivar
	# 		vexp, SNR = sa.autosmooth(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], var_y = var[SN.x1:SN.x2])
	# 		v, si_min_wave = sa.measure_velocity(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], 5900., 6300., vexp=vexp, plot=True)
	# 		print SN.name, v


	# for i, comp in enumerate(composites):
	# 	comp.name = "Comp" + str(i)
	# 	dm15 = np.round(np.nanmean(comp.dm15_array[comp.x1:comp.x2]),2)
	# 	# r = sa.measure_si_ratio(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], vexp = .001, dm15=dm15)
	# 	v_strong, si_min_wave = sa.measure_velocity(comp.wavelength[comp.x1:comp.x2],comp.flux[comp.x1:comp.x2], 5900., 6300.)
	# 	print 'v = ', v_strong


	# vs = []
	# for b in boot_sn_arrays[0]:
	# 	v_strong, si_min_wave = sa.measure_velocity(b.wavelength[b.x1:b.x2],b.flux[b.x1:b.x2], 5900., 6300.)
	# 	vs.append(v_strong)
	# v_err = np.nanstd(vs)
	# print 'v_err = ', v_err


	# r = sa.measure_si_ratio(comp2.wavelength[comp2.x1:comp2.x2], comp2.flux[comp2.x1:comp2.x2], vexp = .001)
	# print r
	# set_min_num_spec(composites, 10)
	set_min_num_spec(composites, 2)
	# set_min_num_spec(composites, 20)
	normalize_comps(composites)
	# for comp in composites:
	# 	sa.measure_si_velocity(comp)
	
	#plotting rutines to visualize composites
	# set_min_num_spec(composites, 5)
	kplot.comparison_plot(composites)
	# scaled_plot(composites, min_num_show=2, min_num_scale=2)
	# stacked_plot(composites)

	# save_comps_to_files(composites)
	
		