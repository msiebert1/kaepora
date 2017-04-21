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
	c_m = matplotlib.cm.viridis
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
	cb = plt.colorbar(s_m, ax = fig.axes)
	cb.set_label('Phase', fontdict = font)
	for ax in fig.axes:
		ax.set_axis_bgcolor('white')

	plt.show()

def stacked_plot(composites):
	"finish"

# set_min_num_spec(composites, 5)

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

	for n in range(num_queries):
		composites.append(composite.main(query_strings[n]))

	composite.optimize_scales(composites, composites[0], True)
	
	comparison_plot(composites)
	scaled_plot(composites)
