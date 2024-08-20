import matplotlib.pyplot as plt
import numpy as np
import kaepora as kpora
import kaepora_plot as kplot
import flux_calibration as fc
import composite
import numpy as np
import pysynphot
from astropy import units as u
from scipy.interpolate import interp2d
import glob
from mpl_toolkits.mplot3d import Axes3D
import gp2d_george as gp2d
import prep 
import pysynphot
import matplotlib
import spectral_analysis as sa

def make_colorbar(params):
    norm = matplotlib.colors.Normalize(vmin=np.min(params),vmax=np.max(params))
    c_m = matplotlib.cm.gist_rainbow
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    return s_m

def create_query_list(phase_range, dm15_range, velocity_range=None, phase_step=1, phase_bin_size=3, phase_step_bin_scale = True):
	queries= kpora.create_query_strings_for_sequence(phase_range, dm15_range, velocity_range=velocity_range, 
															 phase_step = phase_step, phase_bin_size=phase_bin_size, 
															 phase_step_bin_scale = phase_step_bin_scale)
	return queries


def create_composites(queries):
	composites, sn_arrays, og_SN_Arrays, boot_sn_arrays = kpora.make_composite(queries, db_file = '../data/kaepora_v1.2.db', 
																					  shape_param='dm15', boot=False, make_corr=True, 
																					  av_corr=True, medmean=1, verbose=False, gini_balance=True, 
																					  combine=True, scale_region=[4000,9000], get_og_arr = True)
	return composites


def calc_x1_from_comp_list(composites):
	dm15s = []
	composites_custom_dict = {}
	for comp in composites:
		if comp is not None:
		    phase = np.average(comp.phase_array[comp.x1:comp.x2])
		    dm15s.append(np.average(comp.dm15_array[comp.x1:comp.x2]))
		    composites_custom_dict[phase] = comp
	print (np.median(dm15s))
	dm15_med = np.median(dm15s)
	x_1 = sa.calc_x1_from_dm15(dm15_med)
	print (x_1)
	return composites_custom_dict, x_1

def scale_composites_to_photometry(lc_template, composites_custom_dict, phase_range=[-20,45], x_1=0):
	if lc_template == 'salt':
		phot = fc.create_salt_phot(x_1)
	elif lc_template == 'hsiao':
		phot = fc.create_hsiao_phot()

	min_phase = phase_range[0]
	max_phase = phase_range[1]
	del_phases = []
	for phase in composites_custom_dict.keys():
	    composites_custom_dict[phase].event_data = {}
	    composites_custom_dict[phase].event_data['Homogenized_Photometry'] = phot #assign template here
	    if phase < min_phase or phase > max_phase:
	        del_phases.append(phase)
	for p in del_phases:
	    print ("spectrum outside phase range")
	    del composites_custom_dict[p]

	fc.plot_light_curves(phot, lc_template, spec_dates = composites_custom_dict.keys(), fit=True)
	print (composites_custom_dict)
	calibrated_spec_array = fc.scale_composites_to_photometry(composites_custom_dict, verbose=True, template=lc_template)

	fc.plot_spectra(calibrated_spec_array)

	return calibrated_spec_array


def prepare_spec(wrange, calibrated_spec_array):
	spec_arr_new = []
	minwaves = []
	maxwaves = []

	for spec in calibrated_spec_array:
	    if spec.wavelength[spec.x1] < wrange[0] and spec.wavelength[spec.x2] > wrange[1]:
	        spec_arr_new.append(spec)

	specnew, phases, waves = fc.prepare_spec(spec_arr_new,waverange=wrange,plot=True)

	return specnew, phases, waves

def mklecho_tophat(X2, Y2, Z2, p_int = [-10,45]):
	inds = []
	plow = p_int[0]
	phigh = p_int[1]
	phases_interp = []
	for i, y in enumerate(Y2):
	    if y[0] > plow and y[0] < phigh:
	        inds.append(i)
	    phases_interp.append(y[0])
	    
	phase_diff = np.diff(phases_interp)[0]

	lecho = np.nansum(Z2[inds]*phase_diff,axis=0)

	return lecho

def make_le_template(phase_range, dm15_range, lc_template, wrange=[3500, 9000], plot=False):

	queries = create_query_list(phase_range, dm15_range, velocity_range=None, 
								phase_step=1, phase_bin_size=3, phase_step_bin_scale = True)
	
	composites = create_composites(queries)

	composites_custom_dict, x_1 = calc_x1_from_comp_list(composites)

	if lc_template == 'salt':
		phase_range = [-20,90]
	else:
		phase_range = [-20,90]

	calibrated_spec_array = scale_composites_to_photometry(lc_template, composites_custom_dict, phase_range=phase_range, x_1=x_1)

	specnew, phases, waves = prepare_spec(wrange, calibrated_spec_array)

	X2, Y2, Z2 = fc.linear_interpolation(specnew, waves, phases)

	if plot:
		plt.figure(figsize=[10,10])
		plt.imshow(Z2,aspect=20, origin='lower', interpolation='none')
		plt.show()

		for i,z in enumerate(Z2):
		    print (Y2[i][0])
		    plt.plot(X2[i], z)
		    plt.show()

		fig, ax, = kplot.basic_format()

		s_m = make_colorbar(phases)
		for i,z in enumerate(Z2):
		    plt.plot(X2[i], z, color = s_m.to_rgba(Y2[i][0]))
		plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize = 35)
		plt.ylabel('Flux', fontsize = 35)
		cax=ax.inset_axes([1, 0., 0.05, 1])
		plt.colorbar(s_m, cax=cax,label='Phase')
		plt.show()

	return X2, Y2, Z2










