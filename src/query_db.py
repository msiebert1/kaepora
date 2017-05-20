import composite
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation


def make_colorbar(composites):
	phases = []
	for comp in composites:
		phases.append(np.average(comp.phase_array[comp.x1:comp.x2]))

	norm = matplotlib.colors.Normalize(vmin=np.min(phases),vmax=np.max(phases))
	# c_m = matplotlib.cm.viridis
	c_m = matplotlib.cm.winter
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


def comparison_plot(composites):

	h = [3,1,1,1,1,1,1]
	font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 10,
        }

	gs = gridspec.GridSpec(7, 1, height_ratios=h, hspace = .001)
	fig = plt.figure(num = 1, dpi = 100, figsize = [10,10], facecolor = 'w')
	s_m = make_colorbar(composites)

	composites, scales = composite.optimize_scales(composites, composites[0], True)

	for comp in composites:
		phase = np.average(comp.phase_array[comp.x1:comp.x2])

		rel_flux = plt.subplot(gs[0])
		plt.ylabel('Relative Flux', fontdict = font)
		rel_flux.axes.get_yaxis().set_ticks([])
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], color = s_m.to_rgba(phase))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], comp.low_conf[comp.x1:comp.x2], 
				             comp.up_conf[comp.x1:comp.x2], color = s_m.to_rgba(phase), alpha = 0.5)

		plt.subplot(gs[1], sharex = rel_flux)
		plt.ylabel('Variance', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.ivar[comp.x1:comp.x2], color = s_m.to_rgba(phase))

		plt.subplot(gs[2], sharex = rel_flux)
		plt.ylabel('Residuals', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2], color = s_m.to_rgba(phase))
		if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
			low_resid = comp.low_conf[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2]
			up_resid = comp.up_conf[comp.x1:comp.x2] - composites[0].flux[comp.x1:comp.x2]
			plt.fill_between(comp.wavelength[comp.x1:comp.x2], low_resid, up_resid, color = s_m.to_rgba(phase), alpha = 0.5)

		plt.subplot(gs[3], sharex = rel_flux)
		plt.ylabel('Spectra/Bin', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.spec_bin[comp.x1:comp.x2], color = s_m.to_rgba(phase))

		plt.subplot(gs[4], sharex = rel_flux)
		plt.ylabel('Phase [d]', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.phase_array[comp.x1:comp.x2], color = s_m.to_rgba(phase))

		plt.subplot(gs[5], sharex = rel_flux)
		plt.ylabel('$\Delta$m$_{15}$', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.dm15[comp.x1:comp.x2], color = s_m.to_rgba(phase))

		plt.subplot(gs[6], sharex = rel_flux)
		plt.ylabel('Redshift', fontdict = font)
		plt.plot(comp.wavelength[comp.x1:comp.x2], comp.red_array[comp.x1:comp.x2], color = s_m.to_rgba(phase))

	plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
	# cb = plt.colorbar(s_m, ax = fig.axes)
	# cb.set_label('Phase', fontdict = font)
	for ax in fig.axes:
		ax.set_axis_bgcolor('white')

	plt.show()

def stacked_plot(composites):
	fig, ax = plt.subplots(1,1)
	fig.set_size_inches(20, 12., forward = True)
	ax.get_yaxis().set_ticks([])
	composites, scales = composite.optimize_scales(composites, composites[0], True)

	i = 0
	for comp in composites:
		dm15 = np.average(comp.dm15[comp.x1:comp.x2])
		comp.flux = 1.e15*comp.flux 
		ax.plot(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2] - 2.*i, color = '#3F5D7D', linewidth = 5)
		plt.text(9700, comp.flux[comp.x2] - 2*i + 1, dm15)
		i+=1

	plt.ylabel('Relative Flux + const.')
	plt.xlabel('Wavelength (A)')
	plt.show()



def main(num_queries, query_strings):
	# num_queries = int(sys.argv[1])
	# query_strings = sys.argv[2:]

	composites = []
	for n in range(num_queries):
		composites.append(composite.main(query_strings[n]))

	composite.optimize_scales(composites, composites[0], True)
	return composites


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

if __name__ == "__main__":
	composites = []

	num_queries = int(sys.argv[1])
	query_strings = sys.argv[2:]

	# query_strings = [
	# 				 "SELECT * FROM Supernovae where phase < -10",
	#                  "SELECT * FROM Supernovae where phase >= -10 and phase < -9",
	#                  "SELECT * FROM Supernovae where phase >= -9 and phase < -8",
	#                  "SELECT * FROM Supernovae where phase >= -8 and phase < -7",
	#                  "SELECT * FROM Supernovae where phase >= -7 and phase < -6",
	#                  "SELECT * FROM Supernovae where phase >= -6 and phase < -5",
	#                  "SELECT * FROM Supernovae where phase >= -5 and phase < -4",
	#                  "SELECT * FROM Supernovae where phase >= -4 and phase < -3",
	#                  "SELECT * FROM Supernovae where phase >= -3 and phase < -2",
	#                  "SELECT * FROM Supernovae where phase >= -2 and phase < -1",
	#                  "SELECT * FROM Supernovae where phase >= -1 and phase < 0",
	#                  "SELECT * FROM Supernovae where phase >= 0 and phase < 1",
	#                  "SELECT * FROM Supernovae where phase >= 1 and phase < 2",
	#                  "SELECT * FROM Supernovae where phase >= 2 and phase < 3",
	#                  "SELECT * FROM Supernovae where phase >= 3 and phase < 4",
	#                  "SELECT * FROM Supernovae where phase >= 4 and phase < 5",
	#                  "SELECT * FROM Supernovae where phase >= 5 and phase < 6",
	#                  "SELECT * FROM Supernovae where phase >= 6 and phase < 7",
	#                  "SELECT * FROM Supernovae where phase >= 7 and phase < 8",
	#                  "SELECT * FROM Supernovae where phase >= 8 and phase < 9",
	#                  "SELECT * FROM Supernovae where phase >= 9 and phase < 10",
	#                  "SELECT * FROM Supernovae where phase >= 10 and phase < 11",
	#                  "SELECT * FROM Supernovae where phase >= 11 and phase < 12",
	#                  "SELECT * FROM Supernovae where phase >= 12 and phase < 13",
	#                  "SELECT * FROM Supernovae where phase >= 13 and phase < 14",
	#                  "SELECT * FROM Supernovae where phase >= 14 and phase < 15",
	#                  "SELECT * FROM Supernovae where phase >= 15 and phase < 16",
	#                  "SELECT * FROM Supernovae where phase >= 16 and phase < 17",
	#                  "SELECT * FROM Supernovae where phase >= 17 and phase < 18",
	#                  "SELECT * FROM Supernovae where phase >= 18 and phase < 19",
	#                  "SELECT * FROM Supernovae where phase >= 19 and phase < 20",
	#                  "SELECT * FROM Supernovae where phase >= 20 and phase < 21",
	#                  "SELECT * FROM Supernovae where phase >= 21 and phase < 22",
	#                  "SELECT * FROM Supernovae where phase >= 22 and phase < 23",
	#                  "SELECT * FROM Supernovae where phase >= 23 and phase < 24",
	#                  "SELECT * FROM Supernovae where phase >= 24 and phase < 25",
	#                  "SELECT * FROM Supernovae where phase >= 25 and phase < 26",
	#                  "SELECT * FROM Supernovae where phase >= 26 and phase < 27",
	#                  "SELECT * FROM Supernovae where phase >= 27 and phase < 28",
	#                  "SELECT * FROM Supernovae where phase >= 28 and phase < 29",
	#                  "SELECT * FROM Supernovae where phase >= 29 and phase < 30",
	#                  "SELECT * FROM Supernovae where phase >= 30 and phase < 31",
	#                  "SELECT * FROM Supernovae where phase >= 31 and phase < 32",
	#                  "SELECT * FROM Supernovae where phase >= 32 and phase < 33",
	#                  "SELECT * FROM Supernovae where phase >= 33 and phase < 34",
	#                  "SELECT * FROM Supernovae where phase >= 34 and phase < 35",
	#                  "SELECT * FROM Supernovae where phase >= 35 and phase < 36",
	#                  "SELECT * FROM Supernovae where phase >= 36 and phase < 37",
	#                  "SELECT * FROM Supernovae where phase >= 37 and phase < 38",
	#                  "SELECT * FROM Supernovae where phase >= 38 and phase < 39",
	#                  "SELECT * FROM Supernovae where phase >= 39 and phase < 40",
	#                  "SELECT * FROM Supernovae where phase >= 40 and phase < 42",
	#                  "SELECT * FROM Supernovae where phase >= 42 and phase < 44",
	#                  "SELECT * FROM Supernovae where phase >= 44 and phase < 46",
	#                  "SELECT * FROM Supernovae where phase >= 46 and phase < 49",
	#                  "SELECT * FROM Supernovae where phase >= 49 and phase < 52",
	#                  "SELECT * FROM Supernovae where phase >= 52 and phase < 55",
	#                  "SELECT * FROM Supernovae where phase >= 55 and phase < 58",
	#                  "SELECT * FROM Supernovae where phase >= 58 and phase < 61",
	#                  "SELECT * FROM Supernovae where phase >= 61 and phase < 65",
	#                  "SELECT * FROM Supernovae where phase >= 65 and phase < 70",
	#                  "SELECT * FROM Supernovae where phase >= 70 and phase < 78",
	#                  "SELECT * FROM Supernovae where phase >= 78 and phase < 88",
	#                  "SELECT * FROM Supernovae where phase >= 88 and phase < 98",
	#                  "SELECT * FROM Supernovae where phase >= 98 and phase < 125",
	#                  "SELECT * FROM Supernovae where phase >= 125 and phase < 180",
	#                  "SELECT * FROM Supernovae where phase >= 180"
	#                  ]
	num_queries = len(query_strings)

	for n in range(num_queries):
		composites.append(composite.main(query_strings[n]))

	composite.optimize_scales(composites, composites[0], True)
	
	set_min_num_spec(composites, 5)
	comparison_plot(composites)
	# scaled_plot(composites)
	stacked_plot(composites)