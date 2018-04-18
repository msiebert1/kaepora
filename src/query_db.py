import composite
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import argparse
import spectral_analysis as sa


def make_colorbar(composites):
	params = []
	# for comp in composites:
	# 	params.append(np.average(comp.dm15_array[comp.x1:comp.x2]))
	# norm = matplotlib.colors.Normalize(vmin=np.min(params),vmax=np.max(params))
	norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites))
	# norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites) + .5)
	# norm = matplotlib.colors.Normalize(vmin=1.,vmax=len(composites) + 1.)
	c_m = matplotlib.cm.coolwarm
	# c_m = matplotlib.cm.plasma
	s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
	s_m.set_array([])

	return s_m

def scaled_plot(composites):

	font = {'family' : 'serif',
		'color'  : 'black',
		'weight' : 'bold',
		'size'   : 10,
		}

	fig = plt.figure(num = 1, dpi = 100, figsize = [10,10], facecolor = 'w')
	s_m = make_colorbar(composites)

	composites, scales = composite.optimize_scales(composites, composites[0], True)

	for comp in composites:
		phase = np.average(comp.phase_array[comp.x1:comp.x2])

		plt.ylabel('Relative Flux', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color = s_m.to_rgba(phase))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
							 comp.up_conf[comp.x1:comp.x2], color = s_m.to_rgba(phase), alpha = 0.5)

	plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
	cb = plt.colorbar(s_m, ax = fig.axes)
	cb.set_label('Phase', fontdict = font)
	plt.show()

def set_min_num_spec(composites, num):

	for comp in composites:
		comp.spec_bin = np.array(comp.spec_bin)
		valid_range = np.where(comp.spec_bin >= num)[0]
		comp.x1, comp.x2 = valid_range[0], valid_range[-1]


def comparison_plot(composites, scale_type = False, min_num_show = 1, min_num_scale = 5):

	# plt.style.use('ggplot')
	colors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]
	h = [3,1,1,1,1,1,1]

	gs = gridspec.GridSpec(6, 1, height_ratios=h, hspace = .001)
	fig = plt.figure(num = 1, dpi = 100, figsize = [10,15])
	plt.rc('font', family='serif')
	s_m = make_colorbar(composites)
	lw = 3

	set_min_num_spec(composites, min_num_scale)
	composites, scales = composite.optimize_scales(composites, composites[0], scale_type)
	set_min_num_spec(composites, min_num_show)


	i = 0
	k = 1
	for comp in composites:
		# param = np.average(comp.dm15_array[comp.x1:comp.x2])
		param = k

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
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=5)
		plt.setp(rel_flux.get_xticklabels(), visible=False)
		plt.ylabel('Relative Flux')
		rel_flux.axes.set_ylim([0, 1.4])
		# plt.gca().axes.yaxis.set_ticklabels([])
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
							 comp.up_conf[comp.x1:comp.x2], color = s_m.to_rgba(param), alpha = 0.5)
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
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=5)
		plt.setp(res.get_xticklabels(), visible=False)
		plt.ylabel('Residuals')
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			low_resid = comp.low_conf[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2]
			up_resid = comp.up_conf[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2]
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, up_resid, color = s_m.to_rgba(param), alpha = 0.5)
			# plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, up_resid, color = colors[i%len(colors)], alpha = 0.5)

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
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=5)
		plt.setp(spec.get_xticklabels(), visible=False)
		plt.ylabel('Spectra/Bin')
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
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=5)
		# avg_phase = np.nanmean(comp.phase_array[comp.x1:comp.x2])
		# avg_phase_rnd = np.round(avg_phase,1)
		# phase.yaxis.set_ticks(np.arange(avg_phase_rnd-1.5, avg_phase_rnd+1.51, .5))
		plt.setp(phase.get_xticklabels(), visible=False)
		plt.ylabel('Phase (d)')
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
			length=10)
		plt.tick_params(
			which='minor', 
			bottom='on', 
			top='on',
			left='on',
			right='on',
			length=5)
		# avg_dm15 = np.nanmean(comp.dm15_array[comp.x1:comp.x2])
		# avg_dm15_rnd = np.round(avg_dm15,2)
		# delt.yaxis.set_ticks(np.arange(avg_dm15_rnd-.1, avg_dm15_rnd+.11, .05))
		plt.setp(delt.get_xticklabels(), visible=False)
		plt.ylabel('$\Delta$m$_{15}$')
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		# delt.axes.set_ylim([avg_dm15 - .2, avg_dm15 + .2])

		z = plt.subplot(gs[5], sharex = rel_flux)
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
		# avg_z = np.nanmean(comp.red_array[comp.x1:comp.x2])
		# avg_z_rnd = np.round(avg_z,3)
		# z.yaxis.set_ticks(np.arange(avg_z_rnd-.006, avg_z_rnd+.0061, .005))
		plt.ylabel('Redshift')
		# plt.plot(comp.wavelength[comp.x1:comp.x2], comp.red_array[comp.x1:comp.x2], color = s_m.to_rgba(param))
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.red_array[comp.x1:comp.x2], linewidth = lw, color = s_m.to_rgba(param))
		# z.axes.set_ylim([avg_z - .01, avg_z + .01])
		z.axes.set_xlim([comp.wavelength[comp.x1]-200., comp.wavelength[comp.x2]+200.])
		labels=z.axes.get_yticks().tolist()
		labels[-1]=''
		# print labels
		z.set_yticklabels(labels)

		labels=res.axes.get_yticks().tolist()
		labels[0]=''
		labels[-1]=''
		res.set_yticklabels(labels)

		labels=phase.axes.get_yticks().tolist()
		labels[0]=''
		labels[-1]=''
		phase.set_yticklabels(labels)

		labels=delt.axes.get_yticks().tolist()
		labels[0]=''
		labels[-1]=''
		delt.set_yticklabels(labels)

		labels=res.axes.get_yticks().tolist()
		labels[0]=''
		labels[-1]=''
		res.set_yticklabels(labels)

		labels=rel_flux.axes.get_yticks().tolist()
		labels[0]=''
		rel_flux.set_yticklabels(labels)

		labels=spec.axes.get_yticks().tolist()
		labels[0]=''
		labels[-1]=''
		spec.set_yticklabels(labels)

		# res.yaxis.set_ticks(np.arange(-.11, .23, .07))
		res.axes.set_ylim([-.5, .5])
		rel_flux.axes.set_ylim([0., 1.1])

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

		z.axes.set_xlim([comp.wavelength[comp.x1]-200., comp.wavelength[comp.x2]+200.])

		i+=1
		k+=1

	plt.xlabel('Rest Wavelength ($\AA$)')
	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('Phase', fontdict = font)
	# for ax in fig.axes:
	# 	ax.set_axis_bgcolor('white')
	# plt.savefig('../../Paper_Drafts/multi_panel_composite.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/multi_panel_host.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../Paper_Drafts/multi_panel_host.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../Paper_Drafts/split_host_lowdm15.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/split_host_lowdm15.png', dpi = 300, bbox_inches = 'tight')
	plt.savefig('../../Paper_Drafts/split_host_highdm15.pdf', dpi = 300, bbox_inches = 'tight')
	plt.savefig('../../FLASH/split_host_highdm15.png', dpi = 300, bbox_inches = 'tight')
	plt.show()

def stacked_plot(composites, boot=False):
	plt.rc('font', family='serif')
	fig, ax = plt.subplots(1,1)
#     ax.get_yaxis().set_ticks([])
#     plt.ylim([-311,12])
	plt.ylabel('Relative Flux + Const.', fontsize = 20)
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
	plt.gca().axes.yaxis.set_ticklabels([])

	i = 0
	color = '#3F5D7D'
	if boot:
		lw=2
	else:
		lw=4
	for comp in composites:
		# if i==1:
		# 	color = 'darkred'
		dm15 = np.average(comp.dm15_array[comp.x1:comp.x2])
#         buff = 200*np.log10(phase+20)
		buff  = .5*i
		ax.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - buff, color = color, linewidth = lw)
		if boot:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2]- buff, comp.up_conf[comp.x1:comp.x2] - buff, alpha = 0.5)
		# plt.text(9000, comp.flux[comp.x2] - i + .1, '$\mathbf{\Delta m_{15}}$ = ' + str(round(dm15, 2)), fontsize=20)
#         plt.title('All Phase Composite Spectra', fontdict = font1, fontsize = 40)
		i+=1
	plt.xlim([5500.,6500.])
	plt.ylim([-1*len(composites)+3.5,.6])
	# plt.savefig('../../Paper_Drafts/dm15_split_max.pdf', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack.png', dpi = 300, bbox_inches = 'tight')
	# plt.savefig('../../FLASH/host_stack_zoom.png', dpi = 300, bbox_inches = 'tight')
	plt.show()

def normalize_comps(composites):
	for comp in composites:
		norm = 1./np.amax(comp.flux[comp.x1:comp.x2])
		comp.flux = norm*comp.flux
		comp.low_conf = norm*comp.low_conf
		comp.up_conf = norm*comp.up_conf
		comp.ivar *= (norm)**2
	return composites

def normalize_comp(comp):
	norm = 1./np.amax(comp.flux[comp.x1:comp.x2])
	comp.flux = norm*comp.flux
	comp.low_conf = norm*comp.low_conf	
	comp.up_conf = norm*comp.up_conf
	comp.ivar *= (norm)**2
	return comp, norm

def main(num_queries, query_strings, boot='nb', medmean = 1, selection = 'max_coverage'):
	# num_queries = int(sys.argv[1])
	# query_strings = sys.argv[2:]

	composites = []
	sn_arrays = []
	boot_sn_arrays = []
	store_boots = True
	for n in range(num_queries):
		comp, arr, boots = composite.main(query_strings[n],boot=boot, medmean = medmean, selection = selection)
		if store_boots:
			boot_sn_arrays.append(boots)
		composites.append(comp)
		sn_arrays.append(arr)

	# composite.optimize_scales(composites, composites[0], True)
	# composites = normalize_comps(composites)

	return composites, sn_arrays, boot_sn_arrays


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

def plot_comp_and_all_spectra(comp, SN_Array):
	# norm = 1./np.amax(comp.flux)
	# comp.flux = comp.flux*norm
	# for SN in SN_Array:
	# 	SN.flux = SN.flux*norm
	# composite.optimize_scales(SN_Array,comp, True)
	plt.rc('font', family='serif')
	fig, ax = plt.subplots(1,1)
	fig.set_size_inches(10, 8, forward = True)
	plt.minorticks_on()
	plt.xticks(fontsize = 20)
	ax.xaxis.set_ticks(np.arange(np.round(comp.wavelength[comp.x1:comp.x2][0],-3), np.round(comp.wavelength[comp.x1:comp.x2][-1],-3),1000))
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
		plt.plot(SN_Array[i].wavelength[SN_Array[i].x1:SN_Array[i].x2], SN_Array[i].flux[SN_Array[i].x1:SN_Array[i].x2])
	plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], 'k', linewidth = 6)
	plt.ylabel('Relative Flux', fontsize = 30)
	plt.xlabel('Rest Wavelength ' + "($\mathrm{\AA}$)", fontsize = 30)
	# "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase >= -3 and phase <= 3 and morphology >= 9"
	# plt.savefig('../../Paper_Drafts/scaled.png', dpi = 300, bbox_inches = 'tight')
	plt.show()


def save_comps_to_files(composites):
	for SN in composites:
		phase = np.round(np.average(SN.phase_array[SN.x1:SN.x2]), 1)
		vel = np.round(np.average(SN.vel[SN.x1:SN.x2]), 1)
		print phase, vel

	#save to file
	for SN in composites:
		phase = np.round(np.average(SN.phase_array[SN.x1:SN.x2]), 1)
		vel = np.round(np.average(SN.vel[SN.x1:SN.x2]), 1)

		if phase >= 0.:
			sign = 'p'
		else:
			sign = 'm'
		abs_phase = np.absolute(phase)
		phase_str = str(abs_phase)

		abs_vel = np.absolute(vel)
		vel_str = str(abs_vel)

		with open('../../' + sign + phase_str + 'days' + '.flm', 'w') as file:
			wave = np.array(SN.wavelength[SN.x1:SN.x2])
			flux = np.array(SN.flux[SN.x1:SN.x2])
			data = np.array([wave,flux])
			data = data.T
			np.savetxt(file, data)

		plt.plot(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2])
	plt.show()

if __name__ == "__main__":
	composites = []
	SN_Arrays = []
	boot_sn_arrays = []
	store_boots = True

	boot = sys.argv[1]
	query_strings = sys.argv[2:]

	num_queries = len(query_strings)

	for n in range(num_queries):
		c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, make_corr=False, multi_epoch=True, selection='max_coverage')
		# c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, selection='max_coverage_choose_uv')
		composites.append(c)
		SN_Arrays.append(sn_arr)
		if store_boots:
			boot_sn_arrays.append(boots)
	# composite.optimize_scales(composites, composites[0], True)

	for i, comp in enumerate(composites):
		plot_comp_and_all_spectra(comp, SN_Arrays[i])
	# comp2 = composites[1]
	# for sn in sn_arr:
	# 	r = sa.measure_si_ratio(sn.wavelength[sn.x1:sn.x2], sn.flux[sn.x1:sn.x2])
	# 	print r

	for comp in composites:
		r = sa.measure_si_ratio(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], vexp = .001)
		print r
	# r = sa.measure_si_ratio(comp2.wavelength[comp2.x1:comp2.x2], comp2.flux[comp2.x1:comp2.x2], vexp = .001)
	# print r
	# set_min_num_spec(composites, 10)
	# set_min_num_spec(composites, 5)
	set_min_num_spec(composites, 1)
	normalize_comps(composites)
	# for comp in composites:
	# 	sa.measure_si_velocity(comp)
	
	#plotting rutines to visualize composites
	# set_min_num_spec(composites, 5)
	comparison_plot(composites)
	# # scaled_plot(composites)
	# stacked_plot(composites)

	# save_comps_to_files(composites)
	
		