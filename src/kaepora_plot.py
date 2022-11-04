import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import numpy as np

import composite

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
		comp.RMSE = norm*comp.RMSE

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

def make_colorbar(composites, cmap_kind='diff'):
	params = []

	if cmap_kind == 'dm15':
		for comp in composites:
			params.append(np.average(comp.dm15_array[comp.x1:comp.x2]))
		norm = matplotlib.colors.Normalize(vmin=np.min(params),vmax=np.max(params))
		c_m = matplotlib.cm.plasma
	elif cmap_kind == 'diff':
		norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites))
		c_m = matplotlib.cm.coolwarm
	elif cmap_kind == 'other':
		norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites))
		c_m = matplotlib.cm.viridis
	# norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites) + .5)
	# norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites) + 1.)
	s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
	s_m.set_array([])

	return s_m

def basic_format(size=[16,12]):
	plt.rc('font', family='serif')
	fig, ax = plt.subplots(1,1)
	fig.set_size_inches(size[0], size[1], forward = True)
	plt.minorticks_on()
	plt.xticks(fontsize = 20)
	plt.yticks(fontsize = 20)
	plt.tick_params(
	    which='major', 
	    bottom='on', 
	    top='on',
	    left='on',
	    right='on',
	    direction='in',
	    length=20)
	plt.tick_params(
	    which='minor', 
	    bottom='on', 
	    top='on',
	    left='on',
	    right='on',
	    direction='in',
	    length=10)
	return fig, ax
	
def scaled_plot(composites, min_num_show = 5, min_num_scale = 5, include_spec_bin = False, scaleto=10., scale_region=None, zoom=True, fs=[15,12], ticks=[13,8], xlim=None, include_phase_dm15=False,
				 legend_labels = None, rm_last_label=False, expand_ratio=False, text = '', dashes = None, savename = None):
	color_dict = {"Comp": (0.2298057, 0.298717966, 0.753683153, 1.0), "Comp0": "#000080", "Comp1": (0.705673158, 0.01555616, 0.150232812, 1.0), "Comp2": "limegreen", "Comp3": "orange", "Hsiao": "orange", "Foley": "crimson", "Nugent": "turquoise", "SALT2": "black"}
	# order_dict = {"Comp": 2, "Comp0": 2, "Comp1": 2, "Comp2": 2,"Hsiao": 3, "Foley08": 1, "Nugent": 4, "SALT2": 5}
	order_dict = {"Comp": 2, "Comp0": 2, "Comp1": 2, "Comp2": 2, "Comp3": 2, "Hsiao": 2.1, "Foley": 1.9, "Nugent": 2.2, "SALT2": 2.3}
	# if not include_spec_bin:
	# 	h = [3,1]
	# 	gs = gridspec.GridSpec(2, 1, height_ratios=h, hspace = .003)
	# 	fig = plt.figure(num = 1, dpi = 100, figsize = [10,5], facecolor = 'w')
	# else:
	# 	h = [3,1,1]
	# 	gs = gridspec.GridSpec(3, 1, height_ratios=h, hspace = .003)
	# 	fig = plt.figure(num = 1, dpi = 100, figsize = [10,6.25], facecolor = 'w')
	# if not zoom:
	# 	h = [3,1,1]
	# 	gs = gridspec.GridSpec(3, 1, height_ratios=h, hspace = .003)
	# 	fig = plt.figure(num = 1, dpi = 100, figsize = [10,6.25], facecolor = 'w')
	# else:
	# 	h = [3,1,1,1]
	# 	gs = gridspec.GridSpec(4, 1, height_ratios=h, hspace = .003)
	# 	fig = plt.figure(num = 1, dpi = 100, figsize = [10,7.5], facecolor = 'w')

	# h = [3,1,1]
	# gs = gridspec.GridSpec(3, 1, height_ratios=h, hspace = .003)
	# fig = plt.figure(num = 1, dpi = 100, figsize = [10,8], facecolor = 'w')

	if include_phase_dm15:
		i = 0
		k = 1

		red_minmax = []
		phase_minmax = []
		spec_minmax = []
		dm15_minmax = []
		ranges = []
		max_range = composites[0]
		min_range = composites[0]
		for comp in composites:
			red_minmax.append([np.nanmin(comp.red_array[comp.x1:comp.x2]), np.nanmax(comp.red_array[comp.x1:comp.x2])])
			phase_minmax.append([np.nanmin(comp.phase_array[comp.x1:comp.x2]), np.nanmax(comp.phase_array[comp.x1:comp.x2])])
			spec_minmax.append([np.nanmin(comp.spec_bin[comp.x1:comp.x2]), np.nanmax(comp.spec_bin[comp.x1:comp.x2])])
			dm15_minmax.append([np.nanmin(comp.dm15_array[comp.x1:comp.x2]), np.nanmax(comp.dm15_array[comp.x1:comp.x2])])
			therange = comp.wavelength[comp.x2] - comp.wavelength[comp.x1]
			if therange > (max_range.wavelength[max_range.x2] - max_range.wavelength[max_range.x1]):
				max_range = comp
			if therange < (min_range.wavelength[min_range.x2] - min_range.wavelength[min_range.x1]):
				min_range = comp

		red_bounds = [None,None]
		phase_bounds = [None,None]
		spec_bounds = [None,None]
		dm15_bounds = [None,None]
		red_bounds[0] = np.nanmin(red_minmax)
		red_bounds[1] = np.nanmax(red_minmax)
		phase_bounds[0] = np.nanmin(phase_minmax)
		phase_bounds[1] = np.nanmax(phase_minmax)
		spec_bounds[0] = np.nanmin(spec_minmax)
		spec_bounds[1] = np.nanmax(spec_minmax)
		dm15_bounds[0] = np.nanmin(dm15_minmax)
		dm15_bounds[1] = np.nanmax(dm15_minmax)

	h = [3,1,1,1]
	gs = gridspec.GridSpec(4, 1, height_ratios=h, hspace = .003)
	fig = plt.figure(num = 1, dpi = 100, figsize = [10,7.5], facecolor = 'w')

	s_m = make_colorbar(composites)

	# composites, scales = composite.optimize_scales(composites, composites[0], True)
	set_min_num_spec([composites[0]], min_num_scale)
	if scale_region:
		normalize_comps(composites,w1=scale_region[0], w2=scale_region[1])
	else:
		normalize_comps(composites)
	composites, scales = composite.optimize_scales(composites, composites[0], True)
	if scale_region:
		normalize_comps(composites,w1=scale_region[0], w2=scale_region[1])
	else:
		normalize_comps(composites)

	set_min_num_spec([composites[0]], min_num_show)

	for i, comp in enumerate(composites):
		comp.flux *= scaleto
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			comp.up_conf *= scaleto
			comp.low_conf *= scaleto

		if i == 0:
			comp.name = "Comp"
		if i == 1:
			comp.name = "Comp1"
		c = color_dict[comp.name]

		if comp.name != "Comp":
			lw = 2
		else:
			lw = 2
		# phase = np.average(comp.phase_array[comp.x1:comp.x2])
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color = s_m.to_rgba(phase))
		# if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
		# 	plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
		# 					 comp.up_conf[comp.x1:comp.x2], color = s_m.to_rgba(phase), alpha = 0.5)
		rel_flux = plt.subplot(gs[0])
		plt.rc('font', family='serif')
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)
		plt.ylabel('Relative Flux', fontsize=fs[0])
		# plt.ylim([-.05*scaleto, scaleto + .15*scaleto])
		# plt.ylim([-.05*scaleto, scaleto + .3*scaleto])
		if dashes:
			plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color=c, linewidth=lw, zorder=order_dict[comp.name], dashes=dashes)
		else:
			plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color=c, linewidth=lw, zorder=order_dict[comp.name])
		# plt.ylim([0.,1.3])
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
							 comp.up_conf[comp.x1:comp.x2], color=c, alpha = 0.5, zorder=order_dict[comp.name])

		rel_flux.set_ylim([0., 1.25*max(composites[0].flux[composites[0].x1:composites[0].x2])])
		plt.setp(rel_flux.get_xticklabels(), visible=False)
		res = plt.subplot(gs[1], sharex = rel_flux)
		plt.rc('font', family='serif')
		res.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		plt.ylabel('Ratio', fontsize=fs[0])
		if comp.name != "Comp2" and comp.name != "SALT2":
			plt.plot(comp.wavelength[composites[0].x1:composites[0].x2], np.clip(comp.flux[composites[0].x1:composites[0].x2]/composites[0].flux[composites[0].x1:composites[0].x2], -1, 10), color=c, linewidth=lw, zorder=order_dict[comp.name])
		else:
			plt.plot(comp.wavelength[comp.x1:comp.x2], np.clip(comp.flux[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2], -1, 10), color=c, linewidth=lw, zorder=order_dict[comp.name])
		# plt.plot(comp.wavelength[composites[0].x1:composites[0].x2], comp.flux[composites[0].x1:composites[0].x2]/composites[0].flux[composites[0].x1:composites[0].x2])
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			low_resid = comp.low_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
			up_resid = comp.up_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
			# low_resid = comp.low_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
			# up_resid = comp.up_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], np.clip(low_resid, -1, 10), np.clip(up_resid, -1, 10), color=c, alpha = 0.5, zorder=order_dict[comp.name])

		if 'comp' in comp.name.lower():
			print ('Phase: ', np.average(comp.phase_array[comp.x1:comp.x2]))
			print ('dm15: ', np.average(comp.dm15_array[comp.x1:comp.x2]))
			print ('Redshift: ', np.nanmean(comp.red_array[comp.x1:comp.x2]))

		if zoom:
			res2 = plt.subplot(gs[2], sharex = rel_flux)
			plt.ylabel('Ratio', fontsize=fs[0])
			if comp.name != "Comp2" and comp.name != "SALT2":
				plt.plot(comp.wavelength[composites[0].x1:composites[0].x2], np.clip(comp.flux[composites[0].x1:composites[0].x2]/composites[0].flux[composites[0].x1:composites[0].x2], -1, 10), color=c, linewidth=lw, zorder=order_dict[comp.name])
			else:
				plt.plot(comp.wavelength[comp.x1:comp.x2], np.clip(comp.flux[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2], -1, 10), color=c, linewidth=lw, zorder=order_dict[comp.name])
			# plt.plot(comp.wavelength[composites[0].x1:composites[0].x2], comp.flux[composites[0].x1:composites[0].x2]/composites[0].flux[composites[0].x1:composites[0].x2])
			if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
				low_resid = comp.low_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
				up_resid = comp.up_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
				# low_resid = comp.low_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
				# up_resid = comp.up_conf[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
				plt.fill_between(comp.wavelength[comp.x1:comp.x2], np.clip(low_resid, -1, 10), np.clip(up_resid, -1, 10), color=c, alpha = 0.5, zorder=order_dict[comp.name])


	majorLocator = MultipleLocator(1)
	majorFormatter = FormatStrFormatter('%d')
	minorLocator = MultipleLocator(.5)

	res.axes.yaxis.set_major_locator(majorLocator)
	res.axes.yaxis.set_major_formatter(majorFormatter)

	res.axes.yaxis.set_minor_locator(minorLocator)

	if zoom:
		res2.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		majorLocator = MultipleLocator(.2)
		majorFormatter = FormatStrFormatter('%.1f')
		minorLocator = MultipleLocator(.1)

		res2.axes.yaxis.set_major_locator(majorLocator)
		res2.axes.yaxis.set_major_formatter(majorFormatter)

		res2.axes.yaxis.set_minor_locator(minorLocator)
	else:
		res.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		majorLocator = MultipleLocator(1)
		majorFormatter = FormatStrFormatter('%.1f')
		minorLocator = MultipleLocator(.5)

		res.axes.yaxis.set_major_locator(majorLocator)
		res.axes.yaxis.set_major_formatter(majorFormatter)

		res.axes.yaxis.set_minor_locator(minorLocator)
		# print res2.axes.get_yticks().tolist()

	# plt.yticks(np.arange(-.25*scaleto, .25*scaleto, .5))
	# plt.ylim([-.25*scaleto, .25*scaleto])
	# labels=res.axes.get_yticks().tolist()
	# labels[0] = ''
	# labels[2] = ''
	# labels[4] = ''
	# labels[6] = ''
	# labels[8] = ''
	# print labels
	# res.set_yticklabels(labels)
	# ticks=res.axes.get_yticks()

	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('Phase', fontdict = font)
	
	if include_spec_bin:
		plt.setp(res.get_xticklabels(), visible=False)
		if zoom:
			plt.setp(res2.get_xticklabels(), visible=False)

	if include_spec_bin:
		if not zoom:
			spec = plt.subplot(gs[2], sharex = rel_flux)
		else:
			spec = plt.subplot(gs[3], sharex = rel_flux)
		spec.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		if max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) > 90: 
			majorLocator = MultipleLocator(30)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(10)
		elif max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) <= 90 and max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) > 45:
			majorLocator = MultipleLocator(20)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(10)
		elif max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) <= 45 and max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) > 23:
			majorLocator = MultipleLocator(10)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(5)
		elif max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) <= 23 and max(composites[0].spec_bin[composites[0].x1:composites[0].x2]) > 10:
			majorLocator = MultipleLocator(5)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(1)
		else:
			majorLocator = MultipleLocator(2)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(1)


		spec.axes.yaxis.set_major_locator(majorLocator)
		spec.axes.yaxis.set_major_formatter(majorFormatter)
		spec.axes.yaxis.set_minor_locator(minorLocator)

		plt.ylim([0., 1.25*max(composites[0].spec_bin[composites[0].x1:composites[0].x2])])
		plt.ylabel('SNe/Bin', fontsize=fs[1])
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(composites[0].wavelength[composites[0].x1:composites[0].x2], composites[0].spec_bin[composites[0].x1:composites[0].x2], color="#000080",linewidth=lw)
		if composites[1].name == "Comp2":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].spec_bin[composites[1].x1:composites[1].x2], color="limegreen",linewidth=lw)
		if composites[1].name == "Comp3":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].spec_bin[composites[1].x1:composites[1].x2], color="orange",linewidth=lw)

		if rm_last_label:
			labels=spec.axes.get_yticks().tolist()
			labels[-2]=''
			for i, l in enumerate(labels):
				if l != '':
					labels[i] = str(int(float(labels[i])))
			spec.set_yticklabels(labels)
		else:
			labels=spec.axes.get_yticks().tolist()
			for i, l in enumerate(labels):
				if l != '':
					labels[i] = str(int(float(labels[i])))
			spec.set_yticklabels(labels)

	if include_phase_dm15:
		plt.setp(res.get_xticklabels(), visible=False)
		if zoom:
			plt.setp(res2.get_xticklabels(), visible=False)

	if include_phase_dm15:
		phase = plt.subplot(gs[2], sharex = rel_flux)
		phase.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		plt.setp(phase.get_xticklabels(), visible=False)
		plt.ylabel('Phase (d)')
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(composites[0].wavelength[composites[0].x1:composites[0].x2], composites[0].phase_array[composites[0].x1:composites[0].x2], color="#000080",linewidth=lw)
		if composites[1].name == "Comp2":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].phase_array[composites[1].x1:composites[1].x2], color="limegreen",linewidth=lw)
		if composites[1].name == "Comp3":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].phase_array[composites[1].x1:composites[1].x2], color="orange",linewidth=lw)

		p = np.average(comp.phase_array[comp.x1:comp.x2])
		majorLocator = MultipleLocator(1.0)
		majorFormatter = FormatStrFormatter('%d')
		minorLocator = MultipleLocator(.5)
		phase.axes.yaxis.set_major_locator(majorLocator)
		phase.axes.yaxis.set_major_formatter(majorFormatter)
		phase.axes.yaxis.set_minor_locator(minorLocator)
		phase.axes.set_ylim([phase_bounds[0] - 1.2, phase_bounds[1] + 1.2])


		delt = plt.subplot(gs[3], sharex = rel_flux)
		delt.axes.spines['top'].set_visible(False)
		plt.minorticks_on()
		plt.xticks(fontsize = 10)
		plt.yticks(fontsize = 9)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[0],
			direction='in',
			zorder=20)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=ticks[1],
			direction='in',
			zorder=20)

		plt.setp(delt.get_xticklabels(), visible=True)
		plt.ylabel('$\Delta$m$_{15}$ (mag)')
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(composites[0].wavelength[composites[0].x1:composites[0].x2], composites[0].dm15_array[composites[0].x1:composites[0].x2], color="#000080",linewidth=lw)
		if composites[1].name == "Comp2":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].dm15_array[composites[1].x1:composites[1].x2], color="limegreen",linewidth=lw)
		if composites[1].name == "Comp3":
			plt.plot(composites[1].wavelength[composites[1].x1:composites[1].x2], composites[1].dm15_array[composites[1].x1:composites[1].x2], color="orange",linewidth=lw)

		dm15 = np.average(comp.dm15_array[comp.x1:comp.x2])
		majorLocator = MultipleLocator(.1)
		majorFormatter = FormatStrFormatter('%0.1f')
		minorLocator = MultipleLocator(.05)
		delt.axes.yaxis.set_major_locator(majorLocator)
		delt.axes.yaxis.set_major_formatter(majorFormatter)
		delt.axes.yaxis.set_minor_locator(minorLocator)
		delt.axes.set_ylim([dm15_bounds[0] - .12, dm15_bounds[1] + .12])

		# if rm_last_label:
		# 	labels=spec.axes.get_yticks().tolist()
		# 	labels[-2]=''
		# 	for i, l in enumerate(labels):
		# 		if l != '':
		# 			labels[i] = str(int(float(labels[i])))
		# 	spec.set_yticklabels(labels)
		# else:
		# 	labels=spec.axes.get_yticks().tolist()
		# 	for i, l in enumerate(labels):
		# 		if l != '':
		# 			labels[i] = str(int(float(labels[i])))
		# 	spec.set_yticklabels(labels)

	if xlim:
		plt.xlim(xlim)
		rel_flux.text(4500, 9.5, '4 days after peak brightness', fontsize=15, horizontalalignment='left')
		rel_flux.text(5330, 5.0, 'S II', fontsize=10, horizontalalignment='left')
		rel_flux.text(3700, 2.8, 'Ca II H&K', fontsize=10, horizontalalignment='left')
		rel_flux.text(6080, 3.0, 'Si II', fontsize=10, horizontalalignment='left')
	else:
		plt.xlim([composites[0].wavelength[composites[0].x1]-100,composites[0].wavelength[composites[0].x2]+100])
	plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize=fs[0])

	if zoom:
		res.set_ylim([-.25,2.25])
	else:
		# res.set_ylim([0.,4.])
		res.set_ylim([0.5,1.5])
	if savename is not None and '91bg' in savename:
		res.set_ylim([0,5.5])
	if savename is not None and 'middm15' in savename:
		rel_flux.set_ylim([-.5,12.5])
	if zoom:
		res2.set_ylim([0.,4.])
	if legend_labels:
		# rel_flux.legend(legend_labels)
		if 'Month' in text:
			rel_flux.legend(legend_labels, bbox_to_anchor=(0.5, 0.45, 0.48, 0.5), fontsize=10, title=text)
		else:
			rel_flux.legend(legend_labels, loc='best', bbox_to_anchor=(0.5, 0.45, 0.48, 0.5), title=text)
	if savename is not None:
		plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
		plt.savefig(savename+'.png', dpi = 300, bbox_inches = 'tight')
	plt.show()

def comparison_plot(composites, scale_type = False, min_num_show = 1, min_num_scale = 2, template=None, scaleto=10., scale_region = None, compare_ind = 0, rm_last_label = False,
					legend_labels = None, title=None, cmap_kind='diff', extra= False, zoom_ratio=False, label_features=False, savename=None, verbose=True):

	# plt.style.use('ggplot')
	colors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]
	h = [3,1,1,1,1,1]

	gs = gridspec.GridSpec(6, 1, height_ratios=h, hspace = .001)
	fig = plt.figure(num = 1, dpi = 100, figsize = [10,15])
	plt.rc('font', family='serif')
	s_m = make_colorbar(composites, cmap_kind = cmap_kind)
	lw = 2

	if scale_region:
		normalize_comps(composites,w1=scale_region[0], w2=scale_region[1])
	else:
		normalize_comps(composites)
	set_min_num_spec(composites, min_num_scale)
	composites, scales = composite.optimize_scales(composites, composites[0], scale_type)
	if scale_region:
		normalize_comps(composites,w1=scale_region[0], w2=scale_region[1])
	else:
		normalize_comps(composites)
	set_min_num_spec(composites, min_num_show)


	i = 0
	k = 1

	red_minmax = []
	phase_minmax = []
	spec_minmax = []
	dm15_minmax = []
	ranges = []
	max_range = composites[0]
	min_range = composites[0]
	for comp in composites:
		red_minmax.append([np.nanmin(comp.red_array[comp.x1:comp.x2]), np.nanmax(comp.red_array[comp.x1:comp.x2])])
		phase_minmax.append([np.nanmin(comp.phase_array[comp.x1:comp.x2]), np.nanmax(comp.phase_array[comp.x1:comp.x2])])
		spec_minmax.append([np.nanmin(comp.spec_bin[comp.x1:comp.x2]), np.nanmax(comp.spec_bin[comp.x1:comp.x2])])
		dm15_minmax.append([np.nanmin(comp.dm15_array[comp.x1:comp.x2]), np.nanmax(comp.dm15_array[comp.x1:comp.x2])])
		therange = comp.wavelength[comp.x2] - comp.wavelength[comp.x1]
		if therange > (max_range.wavelength[max_range.x2] - max_range.wavelength[max_range.x1]):
			max_range = comp
		if therange < (min_range.wavelength[min_range.x2] - min_range.wavelength[min_range.x1]):
			min_range = comp

	red_bounds = [None,None]
	phase_bounds = [None,None]
	spec_bounds = [None,None]
	dm15_bounds = [None,None]
	red_bounds[0] = np.nanmin(red_minmax)
	red_bounds[1] = np.nanmax(red_minmax)
	phase_bounds[0] = np.nanmin(phase_minmax)
	phase_bounds[1] = np.nanmax(phase_minmax)
	spec_bounds[0] = np.nanmin(spec_minmax)
	spec_bounds[1] = np.nanmax(spec_minmax)
	dm15_bounds[0] = np.nanmin(dm15_minmax)
	dm15_bounds[1] = np.nanmax(dm15_minmax)

	for comp in composites:
		comp.flux *= scaleto
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			comp.up_conf *= scaleto
			comp.low_conf *= scaleto

	for comp in composites:
		# param = np.average(comp.dm15_array[comp.x1:comp.x2])
		param = k
		if (comp.x2 - comp.x2) > (max_range.x2 - max_range.x1):
			max_range = comp
		rel_flux = plt.subplot(gs[0])
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		plt.setp(rel_flux.get_xticklabels(), visible=False)
		plt.ylabel('Relative Flux', fontsize=20)
		rel_flux.axes.set_ylim([-.05*scaleto, 1.1*scaleto])
		# plt.gca().axes.yaxis.set_ticklabels([])
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
							 comp.up_conf[comp.x1:comp.x2], color = s_m.to_rgba(param), alpha = 0.5)

		# if comp.RMSE != None:
		# 	up_err = comp.flux + comp.RMSE
		# 	low_err = comp.flux - comp.RMSE
		# 	plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_err[comp.x1:comp.x2], 
		# 						 up_err[comp.x1:comp.x2], color = 'gold', alpha = 0.7)

			# plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
			# 	             comp.up_conf[comp.x1:comp.x2], color = colors[i%len(colors)], alpha = 0.5)

		# var = plt.subplot(gs[1], sharex = rel_flux)
		# plt.setp(var.get_xticklabels(), visible=False)
		# plt.ylabel('Variance', fontdict = font)
		# # plt.plot(comp.wavelength[comp.x1:comp.x2], comp.ivar[comp.x1:comp.x2], color = s_m.to_rgba(param))
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.ivar[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))

		res = plt.subplot(gs[1], sharex = rel_flux)
		# res.yaxis.set_ticks(np.arange(-.11, .23, .07))
		# res.axes.set_ylim([-.18, .24])
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		plt.setp(res.get_xticklabels(), visible=False)
		plt.ylabel('Ratio', fontsize=15)
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], (comp.flux[comp.x1:comp.x2]/composites[compare_ind].flux[comp.x1:comp.x2]), linewidth = lw, color = s_m.to_rgba(param))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			low_resid = comp.low_conf[comp.x1:comp.x2]/composites[compare_ind].flux[comp.x1:comp.x2]
			up_resid = comp.up_conf[comp.x1:comp.x2]/composites[compare_ind].flux[comp.x1:comp.x2]
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, up_resid, color = s_m.to_rgba(param), alpha = 0.5)
			# plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, up_resid, color = colors[i%len(colors)], alpha = 0.5)
		# if comp.RMSE != None:
		# 	low_resid = low_err[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
		# 	up_resid = up_err[comp.x1:comp.x2]/composites[0].flux[comp.x1:comp.x2]
		# 	plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, 
		# 						 up_resid, color = 'gold', alpha = 0.7)

		spec = plt.subplot(gs[2], sharex = rel_flux)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		plt.setp(spec.get_xticklabels(), visible=False)
		plt.ylabel('SNe/Bin', fontsize=15)
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))

		phase = plt.subplot(gs[3], sharex = rel_flux)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		# avg_phase = np.nanmean(comp.phase_array[comp.x1:comp.x2])
		# avg_phase_rnd = np.round(avg_phase,1)
		# phase.yaxis.set_ticks(np.arange(avg_phase_rnd-1.5, avg_phase_rnd+1.51, .5))
		plt.setp(phase.get_xticklabels(), visible=False)
		plt.ylabel('Phase (d)', fontsize=15)
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.phase_array[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.phase_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		# phase.axes.set_ylim([avg_phase - 2., avg_phase + 2.])

		delt = plt.subplot(gs[4], sharex = rel_flux)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		# avg_dm15 = np.nanmean(comp.dm15_array[comp.x1:comp.x2])
		# avg_dm15_rnd = np.round(avg_dm15,2)
		# delt.yaxis.set_ticks(np.arange(avg_dm15_rnd-.1, avg_dm15_rnd+.11, .05))
		plt.setp(delt.get_xticklabels(), visible=False)

		if 'shape_param' in dir(comp):
			plt.ylabel(comp.shape_param, fontsize=15)
		else:
			plt.ylabel('$\Delta$m$_{15}$ (mag)', fontsize=15)
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		# delt.axes.set_ylim([avg_dm15 - .2, avg_dm15 + .2])

		z = plt.subplot(gs[5], sharex = rel_flux)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		plt.yticks(fontsize = 10)
		plt.tick_params(
			which='major', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			direction="in",
			length=5)
		avg_z = np.nanmean(comp.red_array[comp.x1:comp.x2])
		# print avg_z
		# avg_z_rnd = np.round(avg_z,3)
		# z.yaxis.set_ticks(np.arange(avg_z_rnd-.006, avg_z_rnd+.0061, .005))
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.red_array[comp.x1:comp.x2], color = s_m.to_rgba(param))
		if not extra:
			plt.ylabel('Redshift', fontsize=15)
			plt.plot(comp.wavelength[comp.x1:comp.x2], comp.red_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
			if ~np.isnan(avg_z):
				z.axes.set_ylim([avg_z - .015, avg_z + .01])
		else:
			# plt.ylabel('HR (mag)',fontsize=15)
			# avg_m1 = np.nanmean(composites[0].hr_array[comp.x1:comp.x2])
			# avg_m2 = np.nanmean(composites[1].hr_array[comp.x1:comp.x2])
			# avg_m = (avg_m1 + avg_m2)/2.
			# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.hr_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
			# ['E', 'E0', 'E1', 'E2', 'E23', 'E3', 'E6', 'S0', 'S0a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm', 'Sm', 'Sp']
			# z.axes.set_ylim([avg_m - 7, avg_m + 7])
			# y = [5,10,15]
			# labels = ['E','Sa', 'Scd']
			# plt.yticks(y, labels)

			plt.ylabel('Velocity',fontsize=15)
			avg_m1 = np.nanmean(composites[0].vel[comp.x1:comp.x2])
			avg_m2 = np.nanmean(composites[1].vel[comp.x1:comp.x2])
			avg_m = (avg_m1 + avg_m2)/2.
			plt.plot(comp.wavelength[comp.x1:comp.x2], comp.vel[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
			z.axes.set_ylim([avg_m1 - 5, avg_m2 + 5])
		

		z.axes.set_xlim([comp.wavelength[comp.x1]-200., comp.wavelength[comp.x2]+200.])
		# labels=z.axes.get_yticks().tolist()
		# labels[-1]=''
		# print labels
		

		# red = np.average(comp.red_array[comp.x1:comp.x2])
		# # z.set_yticklabels(labels)
		# majorLocator = MultipleLocator(.005)
		# majorFormatter = FormatStrFormatter('%0.3f')
		# minorLocator = MultipleLocator(.0025)
		# z.axes.yaxis.set_major_locator(majorLocator)
		# z.axes.yaxis.set_major_formatter(majorFormatter)
		# z.axes.yaxis.set_minor_locator(minorLocator)
		# z.axes.set_ylim([0., red_bounds[1]+0.003])

		# labels=res.axes.get_yticks().tolist()
		# labels[0]=''
		# labels[-1]=''
		# res.set_yticklabels(labels)
		if not zoom_ratio:
			majorLocator = MultipleLocator(.5)
			majorFormatter = FormatStrFormatter('%.1f')
			minorLocator = MultipleLocator(.25)
			res.axes.yaxis.set_major_locator(majorLocator)
			res.axes.yaxis.set_major_formatter(majorFormatter)
			res.axes.yaxis.set_minor_locator(minorLocator)
			res.axes.set_ylim([0., 2.])
		else:
			majorLocator = MultipleLocator(.5)
			majorFormatter = FormatStrFormatter('%.1f')
			minorLocator = MultipleLocator(.25)
			res.axes.yaxis.set_major_locator(majorLocator)
			res.axes.yaxis.set_major_formatter(majorFormatter)
			res.axes.yaxis.set_minor_locator(minorLocator)
			res.axes.set_ylim([0.4, 1.6])


		# res.axes.set_ylim([-2., 2.])

		if spec_bounds[1] > 45:
			majorLocator = MultipleLocator(20)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(10)
		elif spec_bounds[1] > 20 and spec_bounds[1] <= 45:
			majorLocator = MultipleLocator(10)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(5)
		else:
			majorLocator = MultipleLocator(5)
			majorFormatter = FormatStrFormatter('%d')
			minorLocator = MultipleLocator(2.5)
		spec.axes.yaxis.set_major_locator(majorLocator)
		spec.axes.yaxis.set_major_formatter(majorFormatter)
		spec.axes.yaxis.set_minor_locator(minorLocator)
		spec.axes.spines['top'].set_visible(False)
		spec.axes.set_ylim([0, spec_bounds[1]*1.23])
		if rm_last_label:
			labels=spec.axes.get_yticks().tolist()
			labels[-2]=''
			for i, l in enumerate(labels):
				if l != '':
					labels[i] = str(int(float(labels[i])))
			spec.set_yticklabels(labels)

		p = np.average(comp.phase_array[comp.x1:comp.x2])
		# labels=phase.axes.get_yticks().tolist()
		# labels[0]=''
		# labels[-1]=''
		# phase.set_yticklabels(labels)
		majorLocator = MultipleLocator(1.0)
		majorFormatter = FormatStrFormatter('%d')
		minorLocator = MultipleLocator(.5)
		phase.axes.yaxis.set_major_locator(majorLocator)
		phase.axes.yaxis.set_major_formatter(majorFormatter)
		phase.axes.yaxis.set_minor_locator(minorLocator)
		phase.axes.set_ylim([phase_bounds[0] - 1.2, phase_bounds[1] + 1.2])

		dm15 = np.nanmean(comp.dm15_array[comp.x1:comp.x2])
		# labels=delt.axes.get_yticks().tolist()
		# print labels
		# labels[0]=''
		# labels[-1]=''
		# delt.set_yticklabels(labels)
		majorLocator = MultipleLocator(.1)
		majorFormatter = FormatStrFormatter('%0.1f')
		minorLocator = MultipleLocator(.05)
		delt.axes.yaxis.set_major_locator(majorLocator)
		delt.axes.yaxis.set_major_formatter(majorFormatter)
		delt.axes.yaxis.set_minor_locator(minorLocator)
		delt.axes.set_ylim([dm15_bounds[0] - .1, dm15_bounds[1] + .1])

		# labels=res.axes.get_yticks().tolist()
		# labels[0]=''
		# labels[-1]=''
		# res.set_yticklabels(labels)

		# labels=rel_flux.axes.get_yticks().tolist()
		# labels[0]=''
		# rel_flux.set_yticklabels(labels)

		# labels=spec.axes.get_yticks().tolist()
		# labels[0]=''
		# labels[-1]=''
		# spec.set_yticklabels(labels)

		# res.yaxis.set_ticks(np.arange(-.11, .23, .07))
		# res.axes.set_ylim([-.5, .5])

		# avg_phase = np.nanmean(comp.phase_array[comp.x1:comp.x2])
		# avg_phase_rnd = np.round(avg_phase,1)
		# phase.yaxis.set_ticks(np.arange(avg_phase_rnd-1.5, avg_phase_rnd+1.51, .5))
		# phase.axes.set_ylim([avg_phase - 2., avg_phase + 2.])

		# avg_dm15 = np.nanmean(comp.dm15_array[comp.x1:comp.x2])
		# avg_dm15_rnd = np.round(avg_dm15,2)
		# delt.yaxis.set_ticks(np.arange(avg_dm15_rnd-.1, avg_dm15_rnd+.11, .05))
		# delt.axes.set_ylim([avg_dm15 - .2, avg_dm15 + .2])

		# avg_z = np.nanmean(comp.red_array[comp.x1:comp.x2])
		# avg_z_rnd = np.round(avg_z,3)
		# z.yaxis.set_ticks(np.arange(avg_z_rnd-.006, avg_z_rnd+.0061, .005))
		# z.axes.set_ylim([avg_z - .01, avg_z + .01])

		i+=1
		k+=1
		if verbose:
			print ('Phase: ', np.average(comp.phase_array[comp.x1:comp.x2]))
			print ('dm15: ', np.average(comp.dm15_array[comp.x1:comp.x2]))
			print ('Redshift: ', np.nanmean(comp.red_array[comp.x1:comp.x2]))
			print ('HR: ', np.average(comp.hr_array[comp.x1:comp.x2]))

	z.axes.set_xlim([min_range.wavelength[min_range.x1]-200., min_range.wavelength[min_range.x2]+200.])

	plt.xlabel('Rest Wavelength ($\mathrm{\AA}$)', fontsize=20)
	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('Phase', fontdict = font)
	# for ax in fig.axes:
	# 	ax.set_axis_bgcolor('white')
	# plt.savefig('../../Paper_Drafts/multi_panel_composite.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/multi_panel_host.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../Paper_Drafts/multi_panel_host.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../Paper_Drafts/split_host_lowdm15.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/split_host_lowdm15.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../Paper_Drafts/split_host_highdm15.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/split_host_highdm15.png', dpi = 300, bbox_inches = 'tight')
	rel_flux.set_ylim([0., 1.25*np.nanmax(composites[0].flux[composites[0].x1:composites[0].x2])])
	if legend_labels:
		rel_flux.legend(legend_labels, bbox_to_anchor=(0.48, 0.45, 0.48, 0.5), fontsize=15)
	if title:
		rel_flux.set_title(title)
	if label_features:
		rel_flux.annotate('Si II', xy=(6150, 3.5), xytext=(6550, 4.5),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('S II', xy=(5400, 5), xytext=(5550, 6.5),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('Ca II', xy=(8300, 1), xytext=(8500, 2),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('Ca II', xy=(3800, 3.5), xytext=(3900, 1.5),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('O I', xy=(7450, 1), xytext=(7550, 2),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('Mg II', xy=(4350, 8.5), xytext=(4550, 9.2),
            arrowprops=dict(facecolor='black'), fontsize=15)
		rel_flux.annotate('Si II', xy=(4000, 9.6), xytext=(4250, 10.),
            arrowprops=dict(facecolor='black'), fontsize=15)
	if savename is not None:
		plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
		plt.savefig(savename+'.png', dpi = 300, bbox_inches = 'tight')

	plt.show()

def feature_plot(composites, scale_region = [5600, 6500], boot=False, min_num_show = 1, min_num_scale = 5, 
				 scale=10., cmap_kind = 'dm15', labels = None, savename = None):
	wave1=scale_region[0]
	wave2=scale_region[1]
	set_min_num_spec([composites[0]], min_num_scale)
	normalize_comps(composites)
	composites, scales = composite.optimize_scales(composites, composites[0], True, scale_region=scale_region)
	set_min_num_spec([composites[0]], min_num_show)

	s_m = make_colorbar(composites, cmap_kind=cmap_kind)
	basic_format()
# 	plt.rc('font', family='serif')
# 	fig, ax = plt.subplots(1,1)
# #     ax.get_yaxis().set_ticks([])
# #     plt.ylim([-311,12])
# 	plt.ylabel('Scaled Flux', fontsize = 20)
# 	plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 20)
# 	plt.minorticks_on()
# 	plt.xticks(fontsize = 10)
# 	plt.yticks(fontsize = 10)
# 	plt.tick_params(
# 		which='major', 
# 		bottom='on', 
# 		top='on',
# 		left='on',
# 		right='on',
# 		length=10)
# 	plt.tick_params(
# 		which='minor', 
# 		bottom='on', 
# 		top='on',
# 		left='on',
# 		right='on',
# 		length=5)
# 	fig.set_size_inches(6, 6, forward = True)
	i = 0
	color = '#3F5D7D'
	lw=4
	for i, comp in enumerate(composites):
		roi = np.where((comp.wavelength >= wave1) & (comp.wavelength < wave2))[0]
		r1 = roi[0]
		r2 = roi[-1]
		# if i==1:
		# 	color = 'darkred'
		dm15 = np.average(comp.dm15_array[r1:r2])
		print (dm15)
#         buff = 200*np.log10(phase+20)
		plt.plot(comp.wavelength[r1:r2], scale*comp.flux[r1:r2], color = s_m.to_rgba(i+1), linewidth = lw, label = labels[i])

		if boot:
			low_conf = comp.low_conf[r1:r2]
			up_conf = comp.up_conf[r1:r2]
			plt.fill_between(comp.wavelength[r1:r2], scale*low_conf, scale*up_conf, color = s_m.to_rgba(i+1), alpha = 0.3)
#         plt.title('All Phase Composite Spectra', fontdict = font1, fontsize = 40)
		i+=1
	# plt.xlim([5500.,6500.])
	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('$\Delta m_{15}$ (B) (mag)', fontsize=15)
	plt.xlim([wave1+20, wave2-20])
	plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
	plt.ylabel('Relative Flux', fontsize = 30)
	plt.gca().axvline(x=5690., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)
	plt.gca().axvline(x=5360., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)
	# plt.ylim([0, .6])
	# plt.ylim([-1*len(composites)+3.5,.6])
	# plt.savefig('../../Paper_Drafts/dm15_split_max.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack_zoom.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/individual_spectra.png', dpi = 300, bbox_inches = 'tight')
	plt.legend(bbox_to_anchor=(0.48, 0.45, 0.48, 0.5), fontsize = 30)
	if savename is not None:
		plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
	plt.show()

def feature_panel(composites, wave_range = [5600, 6500], boot=False, min_num_show = 1, min_num_scale = 5, buff = None, pre_buff=False, text=False,
				 scale=10., cmap_kind = 'dm15', labels = None, axarr = None, kj = None, velocity = False, rest_wave = None, savename = None):
	wave1=wave_range[0]
	wave2=wave_range[1]
	set_min_num_spec([composites[0]], min_num_scale)
	normalize_comps(composites)
	composites, scales = composite.optimize_scales(composites, composites[0], True, scale_range=True, wave1=wave1, wave2=wave2)
	set_min_num_spec([composites[0]], min_num_show)


	s_m = make_colorbar(composites, cmap_kind=cmap_kind)

	i = 0
	color = '#3F5D7D'
	lw=2
	k = kj[0]
	j = kj[1]
	axarr[j].axes.yaxis.set_ticklabels([])
	roi = np.where((composites[0].wavelength >= wave1) & (composites[0].wavelength < wave2))[0]
	r1 = roi[0]
	r2 = roi[-1]
	scale = 10*(1./np.amax(composites[0].flux[r1:r2]))
	c = 299792.458
	for i, comp in enumerate(composites):
		roi = np.where((comp.wavelength >= wave1) & (comp.wavelength < wave2))[0]
		r1 = roi[0]
		r2 = roi[-1]

		# if i==1:
		# 	color = 'darkred'
		dm15 = np.average(comp.dm15_array[r1:r2])
#         buff = 200*np.log10(phase+20)
		if buff != None:
			if velocity:
				wave_arr = np.asarray(comp.wavelength[r1:r2])
				vel_arr = -1.*c*((rest_wave/wave_arr)**2. - 1)/(1+((rest_wave/wave_arr)**2.))/10.**3
				axarr[j].plot(vel_arr, scale*comp.flux[r1:r2] - buff, color = s_m.to_rgba(i+1), linewidth = lw, label = labels[i])
			else:
				axarr[j].plot(comp.wavelength[r1:r2], scale*comp.flux[r1:r2] - buff, color = s_m.to_rgba(i+1), linewidth = lw, label = labels[i])
		else:
			if velocity:
				wave_arr = np.asarray(comp.wavelength[r1:r2])
				vel_arr = -1.*c*((rest_wave/wave_arr)**2. - 1)/(1+((rest_wave/wave_arr)**2.))/10.**3
				axarr[k,j].plot(vel_arr, scale*comp.flux[r1:r2], color = s_m.to_rgba(i+1), linewidth = lw, label = labels[i])
			else:
				axarr[k,j].plot(comp.wavelength[r1:r2], scale*comp.flux[r1:r2], color = s_m.to_rgba(i+1), linewidth = lw, label = labels[i])
		if boot:
			low_conf = comp.low_conf[r1:r2]
			up_conf = comp.up_conf[r1:r2]
			if buff != None:
				if velocity:
					wave_arr = np.asarray(comp.wavelength[r1:r2])
					vel_arr = -1.*c*((rest_wave/wave_arr)**2. - 1)/(1+((rest_wave/wave_arr)**2.))/10.**3
					axarr[j].fill_between(vel_arr, scale*low_conf - buff, scale*up_conf - buff, color = s_m.to_rgba(i+1), alpha = 0.3)
				else:
					axarr[j].fill_between(comp.wavelength[r1:r2], scale*low_conf - buff, scale*up_conf - buff, color = s_m.to_rgba(i+1), alpha = 0.3)
			else:
				if velocity:
					wave_arr = np.asarray(comp.wavelength[r1:r2])
					vel_arr = -1.*c*((rest_wave/wave_arr)**2. - 1)/(1+((rest_wave/wave_arr)**2.))/10.**3
					axarr[k,j].fill_between(vel_arr, scale*low_conf, scale*up_conf, color = s_m.to_rgba(i+1), alpha = 0.3)
				else:
					axarr[k,j].fill_between(comp.wavelength[r1:r2], scale*low_conf, scale*up_conf, color = s_m.to_rgba(i+1), alpha = 0.3)
#         plt.title('All Phase Composite Spectra', fontdict = font1, fontsize = 40)
	if text:
		plt.text(4080, np.mean(composites[0].flux[r1:r2])-1.3*buff+9, text, fontsize=15, horizontalalignment='right')
		plt.text(8680, np.mean(composites[0].flux[r1:r2])-1.3*buff+9, text, fontsize=15, horizontalalignment='left')
		# plt.text(5800, np.mean(composites[0].flux[r1:r2])-1.3*buff+9, text, fontsize=15, horizontalalignment='right')
		# plt.text(3.5, np.mean(composites[0].flux[r1:r2])-1.3*buff+9, text, fontsize=15, horizontalalignment='left')
		i+=1
	else:
		plt.gca().axvline(x=4610., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)
		plt.gca().axvline(x=4260., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)

	if pre_buff:
		plt.ylim([-16,12])
	# plt.xlim([5500.,6500.])
	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('$\Delta m_{15}$ (B) (mag)', fontsize=15)
	# axarr[i].set_ylim([wave1+20, wave2-20])
	# plt.ylim([0, .6])
	# plt.ylim([-1*len(composites)+3.5,.6])
	# plt.savefig('../../Paper_Drafts/dm15_split_max.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack_zoom.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/individual_spectra.png', dpi = 300, bbox_inches = 'tight')

def stacked_plot(composites, boot=False, savename = None):
	plt.rc('font', family='serif')
	fig, ax = plt.subplots(1,1)
#     ax.get_yaxis().set_ticks([])
#     plt.ylim([-311,12])
	plt.ylabel('Relative Flux + Constant', fontsize = 20)
	plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 20)
	plt.minorticks_on()
	plt.xticks(fontsize = 10)
	plt.yticks(fontsize = 10)
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
	fig.set_size_inches(3.5, 13.5, forward = True)
	fig.set_size_inches(10.5, 13.5, forward = True)
	plt.gca().axes.yaxis.set_ticklabels([])

	i = 0
	scale = 0.5*.8
	color = '#3F5D7D'
	if boot:
		lw=2
	else:
		lw=2
	for comp in composites:
		if i ==5:
			scale = .5*.8
		# if i==1:
		# 	color = 'darkred'
		dm15 = np.average(comp.dm15_array[comp.x1:comp.x2])
#         buff = 200*np.log10(phase+20)
		buff  = i*scale
		ax.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - buff, color = color, linewidth = lw)
		if boot:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2]- buff, comp.up_conf[comp.x1:comp.x2] - buff, alpha = 0.5)
		plt.text(8500, comp.flux[comp.x2-1130] - i*scale + .1*scale , r'$\langle$$\mathbf{\Delta m_{15}} $(B)$\rangle$ = ' + str(round(dm15, 2)), fontsize=13)
#         plt.title('All Phase Composite Spectra', fontdict = font1, fontsize = 40)
		i+=1
	# plt.xlim([5500.,6500.])
	plt.xlim([composites[0].wavelength[composites[0].x1]-200., composites[0].wavelength[composites[0].x2]+200.])
	# plt.ylim([-1*len(composites)+3.5,.6])
	plt.ylim([(-1*len(composites)+1)*scale,1.1])
	# plt.savefig('../../Paper_Drafts/dm15_split_max.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack_zoom.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/individual_spectra.png', dpi = 300, bbox_inches = 'tight')
	if savename is not None:
		plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
	plt.show()

def make_animation(composites):
	fig, ax = plt.subplots()
	plt.xlim([3000,10000])
	plt.ylim([0.,5.e-15])
	s_m = make_colorbar(composites)

	phase = np.average(composites[0].phase_array[composites[0].x1:composites[0].x2])
	line, = ax.plot(composites[0].wavelength[composites[0].x1:composites[0].x2], 
					composites[0].flux[composites[0].x1:composites[0].x2], color = s_m.to_rgba(phase))

	def animate(i):
		if i == 0:
			ax.cla()
		phase = np.average(composites[i].phase_array[composites[i].x1:composites[i].x2])
		ax.plot(composites[i].wavelength[composites[i].x1:composites[i].x2], 
					composites[i].flux[composites[i].x1:composites[i].x2], color = s_m.to_rgba(phase))
		plt.xlim([3000,10000])
		plt.ylim([0.,5.e-15])
		# line.set_xdata(composites[i].wavelength[composites[i].x1:composites[i].x2])
		# line.set_ydata(composites[i].flux[composites[i].x1:composites[i].x2])  # update the data
		return line,

	animation.FuncAnimation(fig, animate, np.arange(0, len(composites)), repeat = True)
	plt.show()

def plot_comp_and_all_spectra(comp, SN_Array, show_ivar = False, comp2=None, comp3=None, one_color=False, xlim = None, ylim=None, scale_region=None,
							  boots=False, scaleto=10, dm15=False, comp_name = "Inverse-Variance Weighted \nComposite Spectrum", savename=None):
	# norm = 1./np.amax(comp.flux)
	# comp.flux = comp.flux*norm
	# for SN in SN_Array:
	# 	SN.flux = SN.flux*norm
	composite.optimize_scales(SN_Array,comp, True, scale_region=scale_region)
	# for comp in composites:
	# 	comp.flux *= scaleto
	# 	if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
	# 		comp.up_conf *= scaleto
	# 		comp.low_conf *= scaleto
	# for SN in SN_Array:
	# 	SN.flux *= scaleto
	# 	SN.ivar /= (scaleto)**2.

	if not show_ivar:
		if dm15:
			h = [4,1]
			gs = gridspec.GridSpec(2, 1, height_ratios=h, hspace = .003)
			fig = plt.figure(num = 1, dpi = 100, figsize = [12,10], facecolor = 'w')

		rel_flux = plt.subplot(gs[0])
		plt.rc('font', family='serif')
		# fig, ax = plt.subplots(1,1)
		# fig.set_size_inches(10, 8, forward = True)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		# ax.xaxis.set_ticks(np.arange(np.round(comp.wavelength[comp.x1:comp.x2][0],-3), np.round(comp.wavelength[comp.x1:comp.x2][-1],-3),1000))
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
		max_snr = SN_Array[0]
		max_snr_num = np.sum(SN_Array[0].ivar[np.where((SN_Array[0].wavelength >= 2000) & (SN_Array[0].wavelength < 3000))[0]])
		for i in range(len(SN_Array)):
			# if SN_Array[i].SNR > max_snr.SNR:
			# 	max_snr = SN_Array[i]
			weight = np.sum(SN_Array[i].ivar[np.where((SN_Array[i].wavelength >= 2000) & (SN_Array[i].wavelength < 3000))[0]])
			if weight > max_snr_num:
				max_snr = SN_Array[i]
				max_snr_num = weight
		print (max_snr.name)
		for i in range(len(SN_Array)):
			# plt.plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].flux[SN_Array[i].x1:SN_Array[i].x2], color = '#7570b3', alpha = .5)
			# plt.plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].flux[SN_Array[i].x1:SN_Array[i].x2])
			if one_color:
				# plt.plot(SN_Array[i].wavelength, SN_Array[i].flux, color='#7570b3', alpha = .7)
				if SN_Array[i].name == max_snr.name:
					plt.plot(SN_Array[i].wavelength, 10*SN_Array[i].flux, color='crimson', linewidth=1,dashes=[4, 2], zorder=5, label="Highest Weighted Spectrum: \n" + max_snr.name)
				else:
					if i==0:
						plt.plot(SN_Array[i].wavelength, 10*SN_Array[i].flux, color='black', alpha = .2, zorder=0, label='Individual Spectra')
					else:
						plt.plot(SN_Array[i].wavelength, 10*SN_Array[i].flux, color='black', alpha = .2, zorder=0)
			else:
				plt.plot(SN_Array[i].wavelength, 10*SN_Array[i].flux, drawstyle='steps-mid', label=SN_Array[i].name)
		plt.plot(comp.wavelength[comp.x1:comp.x2], 10*comp.flux[comp.x1:comp.x2], 'orange', linewidth = 3, dashes=[12, 2], zorder=3, label=comp_name)
		# if boots:
		# 	plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
		# 				 comp.up_conf[comp.x1:comp.x2], color='orange', alpha = 0.3, zorder=3)
		if comp2:
			plt.plot(comp2.wavelength[comp2.x1:comp2.x2], 10*comp2.flux[comp2.x1:comp2.x2], color='purple', linewidth = 3, dashes=[12, 2], zorder=4, label="Gini Re-weighted \nComposite Spectrum")
			# if boots:
			# 	plt.fill_between(comp2.wavelength[comp2.x1:comp2.x2], comp2.low_conf[comp2.x1:comp2.x2], 
			# 				 comp2.up_conf[comp2.x1:comp2.x2], color='purple', alpha = 0.3,zorder=4)
		if comp3:
			plt.plot(comp3.wavelength[comp3.x1:comp3.x2], 10*comp3.flux[comp3.x1:comp3.x2], color='turquoise', linewidth = 3, dashes=[12, 2], zorder=2, label="Median \nComposite Spectrum")
			# if boots:
			# 	plt.fill_between(comp3.wavelength[comp3.x1:comp3.x2], comp3.low_conf[comp3.x1:comp3.x2], 
			# 				 comp3.up_conf[comp3.x1:comp3.x2], color='turquoise', alpha = 0.3,zorder=2)

		# plt.ylim(0,80)
		plt.rc('font', family='serif')
		plt.ylabel('Relative Flux', fontsize = 30)
		if xlim:
			plt.xlim(xlim)
		if ylim:
			plt.ylim(ylim)
		# "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase >= -3 and phase <= 3 and morphology >= 9"
		# plt.savefig('../../Paper_Drafts/scaled.png', dpi = 300, bbox_inches = 'tight')
		
		plt.legend(fontsize=15)

		# dm15_avg = plt.subplot(gs[1], sharex = rel_flux)
		# plt.minorticks_on()
		# plt.ylabel('$\Delta$m$_{15}$ (mag)', fontsize = 20)
		# plt.xticks(fontsize = 20)
		# # ax.xaxis.set_ticks(np.arange(np.round(comp.wavelength[comp.x1:comp.x2][0],-3), np.round(comp.wavelength[comp.x1:comp.x2][-1],-3),1000))
		# plt.yticks(fontsize = 20)
		# plt.tick_params(
		# 	which='major', 
		# 	bottom='on', 
		# 	top='on',
		# 	left='on',
		# 	right='on',
		# 	length=10)
		# plt.tick_params(
		# 	which='minor', 
		# 	bottom='on', 
		# 	top='on',
		# 	left='on',
		# 	right='on',
		# 	length=5)

		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15_array[comp.x1:comp.x2], 'orange', linewidth = 3, zorder=3)
		# if comp2:
		# 	plt.plot(comp2.wavelength[comp2.x1:comp2.x2], comp2.dm15_array[comp2.x1:comp2.x2], color='purple', linewidth = 3, zorder=4)
		# if comp3:
		# 	plt.plot(comp3.wavelength[comp3.x1:comp3.x2], comp3.dm15_array[comp3.x1:comp3.x2], color='turquoise', linewidth = 3, zorder=2)
		# if xlim:
		# 	plt.xlim(xlim)

		# majorLocator = MultipleLocator(.1)
		# majorFormatter = FormatStrFormatter('%0.1f')
		# minorLocator = MultipleLocator(.05)
		# dm15_avg.axes.yaxis.set_major_locator(majorLocator)
		# dm15_avg.axes.yaxis.set_major_formatter(majorFormatter)
		# dm15_avg.axes.yaxis.set_minor_locator(minorLocator)
		# dm15_avg.axes.set_ylim([.8, 1.2])
		plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
		plt.gca().get_yaxis().set_ticks([])
		# plt.xlim(6500,8000)

		if savename is not None:
			plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
		plt.show()
	else:

		plt.rc('font', family='serif')
		fig, ax = plt.subplots(2,1)
		fig.set_size_inches(10, 8, forward = True)
		plt.minorticks_on()
		plt.xticks(fontsize = 20)
		# ax[1].xaxis.set_ticks(np.arange(np.round(comp.wavelength[comp.x1:comp.x2][0],-3), np.round(comp.wavelength[comp.x1:comp.x2][-1],-3),1000))
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
		for i in range(len(SN_Array)):
			# plt.plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].flux[SN_Array[i].x1:SN_Array[i].x2], color = '#7570b3', alpha = .5)
			# ax[0].plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].flux[SN_Array[i].x1:SN_Array[i].x2], drawstyle='steps-mid')
			# ax[1].plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].ivar[SN_Array[i].x1:SN_Array[i].x2], drawstyle='steps-mid')
			ax[0].plot(SN_Array[i].wavelength, SN_Array[i].flux, drawstyle='steps-mid', linewidth = .3)
			ax[1].plot(SN_Array[i].wavelength, SN_Array[i].ivar, drawstyle='steps-mid', linewidth = .3)
		ax[0].plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], 'k', linewidth = 2, drawstyle='steps-mid')
		ax[0].axvline(x=5890.)
		ax[0].set_ylabel('Relative Flux', fontsize = 30)
		ax[1].set_xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
		# "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase >= -3 and phase <= 3 and morphology >= 9"
		plt.show()