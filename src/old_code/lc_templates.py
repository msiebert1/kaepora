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

def create_query_list(phase_range, dm15_range, velocity_range=None, phase_step=1, phase_bin_size=3, phase_step_bin_scale = True):
	queries= kpora.create_query_strings_for_sequence(phase_range, dm15_range, velocity_range=velocity_range, 
															 phase_step = phase_step, phase_bin_size=phase_bin_size, 
															 phase_step_bin_scale = phase_step_bin_scale)
	return queries


def create_composites(queries)
	composites, sn_arrays, og_SN_Arrays, boot_sn_arrays = kpora.make_composite(queries, db_file = '../data/kaepora_v1.2.db', 
																					  shape_param='dm15', boot=False, make_corr=True, 
																					  av_corr=True, medmean=1, verbose=False, gini_balance=True, 
																					  combine=True, scale_region=[4000,9000], get_og_arr = True)
	return composites


def calc_x1_from_comp_list(composites):
	dm15s = []
	composites_custom_dict = {}
	for comp in composites:
	    phase = np.average(comp.phase_array[comp.x1:comp.x2])
	    dm15s.append(np.average(comp.dm15_array[comp.x1:comp.x2]))
	    composites_custom_dict[phase] = comp
	print (np.median(dm15s))
	dm15_med = np.median(dm15s)
	x_1 = sa.calc_x1_from_dm15(dm15_med)
	print (x_1)
	return x_1

def make_lc_template(lc_template):
	if lc_template == 'salt':
		phot = fc.create_salt_phot(x_1)
	elif lc_template == 'hsiao':
		phot = fc.create_hsiao_phot()

	min_phase = -20
	max_phase = 90
	del_phases = []
	for phase in composites_custom_dict.keys():
	    composites_custom_dict[phase].event_data = {}
	    composites_custom_dict[phase].event_data['Homogenized_Photometry'] = phot #assign template here
	    if phase < min_phase or phase > max_phase:
	        del_phases.append(phase)
	for p in del_phases:
	    print ("spectrum outside phase range")
	    del composites_custom_dict[p]
	return phot, composites_custom_dict




def make_lc_template():
	queries = create_query_list(phase_range, dm15_range, velocity_range=None, 
								phase_step=1, phase_bin_size=3, phase_step_bin_scale = True)
	
	composites = create_composites(queries)

	x_1 = calc_x1_from_comp_list(composites)










