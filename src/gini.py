import numpy as np
def make_gini_ranges():
    gini_ranges = []
    
#     x = np.linspace(1000, 12000, num = 45)
    # x = np.linspace(1000, 12000, num = 23)
    x = np.linspace(1000, 12000, num = 12)
    # x = np.linspace(1000, 12000, num = 2)
    for i in range(len(x)-1):
        gini_ranges.append((x[i], x[i+1]))
    # print gini_ranges
    return gini_ranges

def set_SN_gweights(sn_arr, gini_ranges):
    for SN in sn_arr:
        g_weights = []
        for g_r in gini_ranges:
            g_locs = np.where((SN.wavelength >= g_r[0]) & (SN.wavelength < g_r[1]))[0]
            g_w = np.sum(SN.ivar[g_locs])
            if g_w == 0.:
                g_w = np.nan
            g_weights.append(g_w)
        SN.g_array = g_weights

def calc_gini_coeffs(sn_arr):
    num_specs = []
    for i in range(len(sn_arr[0].g_array)):
        num_specs.append(len(sn_arr))
        for SN in sn_arr:
            if np.isnan(SN.g_array[i]):
                num_specs[i] -= 1
    # print num_specs
    
    gini_coeffs = []
    for i in range(len(sn_arr[0].g_array)):
        gini_num = 0.
        gini_denom = 0.
        for SNi in sn_arr:
            g_wi = SNi.g_array[i]
            jsum = 0.
            for SNj in sn_arr:
                g_wj = SNj.g_array[i]
                g_diff = np.absolute(g_wi - g_wj)
                jsum = np.nansum([jsum, g_diff])
            gini_num = np.nansum([gini_num,jsum])
            gini_denom = np.nansum([gini_denom, g_wi])
        gini_denom *= 2.*num_specs[i]
        gini_coeffs.append(gini_num/gini_denom)
    return gini_coeffs, num_specs

def gini_coeffs(sn_arr):
	gini_ranges = make_gini_ranges()
	set_SN_gweights(sn_arr, gini_ranges)
	gini_coeffs, num_specs = calc_gini_coeffs(sn_arr)
	return gini_coeffs, num_specs, gini_ranges

def calc_deweight_ranges(sn_arr, gini_coeffs, gini_ranges, gini_range_meds, tol=.5):
	deweight_spec = []
	deweight_ranges = []
	replace_ivar_sns = []
	swaps = []
	# max_weight_list = []

	# g_list = []
	# sum_gs = []
	# for j, SN in enumerate(sn_arr):
	#     g_list.append(SN.g_array[i])
	#     sum_gs.append(np.nansum(SN.g_array[i]))
	# g_T = np.transpose(g_list)

	# print len(g_T), len(gini_coeffs)

	# for i in range(len(gini_coeffs)):
	#     if gini_coeffs[i] >= tol:
	#     	max_weight_SN_ind = np.argmax(g_T[i])
	#     	print sn_arr[max_weight_SN_ind].name
	# raise TypeError
	biasing_SNs = []
	biasing_tuples = []

	#first find where coeffs are greater than tolerance
	for i in range(len(gini_coeffs)):
	    if gini_coeffs[i] >= tol:
	        g_list = []
	        #populate g_list with individual SN weights in this wave range
	        for j, SN in enumerate(sn_arr):
	            g_list.append(SN.g_array[i])
	        #if SN has no weight, set to nan
	        for k, g in enumerate(g_list):
	            if np.isnan(g_list[k]):
	                g_list[k] = 0.
	        #sort list to pick out SN with the largest weight in this wave range
	        g_sort = np.argsort(g_list)
	        g_sort_vals = np.sort(g_list)
	        biasing_tuples.append((i, sn_arr[g_sort[-1]])) #store index of wave range and biasing SN


	#biasing_tuples contains most weighted SN for each gini coeff > tol

	# print g_sort
	scale_dict = {}
	scale_ref_dict = {}
	deweight_SNs = []
	for tup in biasing_tuples:
		# print tup[1].g_array[tup[0]]
		g = gini_ranges[tup[0]] # wavelength range
		g_locs = np.where((tup[1].wavelength >= g[0]) & (tup[1].wavelength < g[1]))[0] #indices valid for this gini range
		non_zero_locs = np.where(tup[1].ivar[g_locs] != 0.)
		# print tup[0], tup[1].name, gini_range_meds[tup[0]], np.nansum(tup[1].ivar[g_locs][non_zero_locs]), gini_range_meds[tup[0]]/np.nansum(tup[1].ivar[g_locs][non_zero_locs])
		if tup[1].name in scale_dict:
			med_in_range = np.nansum(tup[1].ivar[g_locs][non_zero_locs])
			scale_dict[tup[1].name] = np.amax([scale_dict[tup[1].name],med_in_range]) #store max of biasing SN in the unbalanced region
			scale_ref_dict[tup[1].name] = np.amax([scale_ref_dict[tup[1].name],gini_range_meds[tup[0]]]) #store max of all SNe in the unbalanced region
		else:
			scale_dict[tup[1].name] = np.nansum(tup[1].ivar[g_locs][non_zero_locs])
			scale_ref_dict[tup[1].name] = gini_range_meds[tup[0]]
		deweight_SNs.append(tup[1])
	# for el in scale_dict.keys():
	# 	print el, scale_dict[el], 'Ref: ', scale_ref_dict[el], 'Fraction: ', scale_dict[el]/scale_ref_dict[el]

	deweight_SNs = set(deweight_SNs)
	# for SN in deweight_SNs:
	# 	print SN.name
	# scale_to = np.median(tot_weights)
	# print scale_to
	# tot_sort = np.argsort(tot_weights)
	# print 
	# for SN in deweight_SNs:
	# 	print SN.name, np.nanmedian(SN.ivar[SN.x1:SN.x2])
	# raise TypeError

	# s = 0
	# found_scale_spec = False
	# while not found_scale_spec:
	# 	s+=1
	# 	med_weights = []
	# 	for SN in deweight_SNs:
	# 		med_weights.append(np.nanmedian(SN.ivar[SN.x1:SN.x2]))

	# 	temp = sn_arr[tot_sort[-s]]
	# 	if sn_arr[tot_sort[-s]] not in deweight_SNs and all(np.nanmedian(temp.ivar[temp.x1:temp.x2]) < w for w in med_weights):
	# 		found_scale_spec = True
	# 		scale_ivar_sn = sn_arr[tot_sort[-s]]

	# for SN in deweight_SNs:
		# print 'Scaling', SN.name, 'to', sn_arr[tot_sort[-s]].name, 'Fraction: ', np.nanmedian(SN.ivar[SN.x1:SN.x2])/np.nanmedian(scale_ivar_sn.ivar[scale_ivar_sn.x1:scale_ivar_sn.x2]), np.nanmedian(SN.ivar[SN.x1:SN.x2])
		
		# if len(prev_swaps) < 1:
		# 	found_scale_spec = True
		# elif (sn_arr[g_sort[-1]].name == prev_swaps[i][0] and sn_arr[g_sort[-s]].name != prev_swaps[i][1]) and g_list[g_sort[-1]]/g_list[g_sort[-s]] > 1.01:
		# 	found_scale_spec = True
		# elif (sn_arr[g_sort[-1]].name != prev_swaps[i][0] and sn_arr[g_sort[-s]].name != prev_swaps[i][1]) and g_list[g_sort[-1]]/g_list[g_sort[-s]] > 1.01:
		# 	found_scale_spec = True

	#         print 'Scale ', sn_arr[g_sort[-1]].name, ' with ', sn_arr[g_sort[-s]].name, 'Fraction: ', g_list[g_sort[-1]]/g_list[g_sort[-s]], 'in range', gini_ranges[i], s
	#         # print 'Deweighting ', sn_arr[g_sort[-1]].name, '  by factor of 1/2 in range ', gini_ranges[i]
	#         for g in g_list:
	#         	tot_weight = np.nansum(g)
	#         max_weight_list.append(sn_arr[g_sort[-1]])

	#         swaps.append((sn_arr[g_sort[-1]].name, sn_arr[g_sort[-s]].name))
	#         scale_ivar_sns.append(sn_arr[g_sort[-s]])
	#         max_ind = np.nanargmax(g_list)
	#         deweight_spec.append(sn_arr[max_ind])
	#         deweight_ranges.append(gini_ranges[i])
	#     else:
	#     	swaps.append((None,None))
	# return deweight_SNs, scale_ivar_sn, deweight_ranges, scale_to
	return deweight_SNs, scale_dict, scale_ref_dict

def deweight_biasing_SNe(deweight_SNs, scale_dict, scale_ref_dict):
	for SN in deweight_SNs:
		# scale_median = np.nanmedian(scale_ivar_sn.ivar[scale_ivar_sn.x1:scale_ivar_sn.x2])
		# sn_median = np.nanmedian(SN.ivar[SN.x1:SN.x2])
		# print sn_median/scale_to
		# SN.ivar = SN.ivar*(scale_median/sn_median)

		# print SN.name, scale_ref_dict[SN.name], scale_dict[SN.name], scale_ref_dict[SN.name]/scale_dict[SN.name]
		scale = scale_ref_dict[SN.name]/scale_dict[SN.name]
		SN.ivar = SN.ivar*scale
	# print 