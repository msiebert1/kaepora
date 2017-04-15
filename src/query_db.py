import composite
import Plotting
from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import matplotlib.gridspec as gridspec

num_queries = int(sys.argv[1])
query_strings = sys.argv[2:]

composites = []
for n in range(num_queries):
	composites.append(composite.main(query_strings[n]))

def make_colorbar(composites):
	phases = []
	for comp in composites:
		phases.append(np.average(comp.phase_array[comp.x1:comp.x2]))

	norm = matplotlib.colors.Normalize(vmin=np.min(phases),vmax=np.max(phases))
	c_m = matplotlib.cm.plasma
	s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
	s_m.set_array([])

	return s_m

#set x1 and x2 to where spec_bin > X
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
		plt.ylabel('Residuals', fontdict = font)
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

comparison_plot(composites)