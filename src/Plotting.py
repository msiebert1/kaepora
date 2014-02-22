import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter

###########################################
#             HOW TO USE
## At the top of YOUR code write the following
# import ../../src/Plotting
# 
## Place after you have calculated your data and fill in information
# 
# plot_data    = [wavelengths, avg_flux, plus, minus]       # Expects [X, Y, Pos_Y Err, Neg_Y Err]
# title        = "Super Cool"                               # Within the quotes place plot title
# image_title  = "../personal/YOUR_REP/WHOA.png"            # Name the image
# xlabel       = "Wavelength ($\AA$)"                       # Label the X axis
# ylabel       = "COOLLLLLL"				    # Label the Y axis
# legend1      = "AVG Spectrum"                             # Name the title you wish to see in the legend
# 
## The following line will plot the data
#
# Plotting.main(plot_data,title, image_title , xlabel, ylabel, legend1)
#
#
###########################################
def main(plot_data, title, image_title, xlabel, ylabel, legend1):

# re-name variables
	xaxis = plot_data[0] 
	yaxis = plot_data[1]
	err_p = plot_data[2] 
	err_n = plot_data[3]

# Legend parameters, must come before plotting
	params = {'legend.fontsize': 15, 'legend.linewidth':2}
	plt.rcParams.update(params)

# Create the plot
	ax1 = subplot(1,1,1)

# Read legend labels from files
#	ax1.plot(xaxis[:],yaxis[:],'k',label = name[5] ) # (X_comp, Y_comp, 'color', curve name)
#	ax1.plot(xaxis[:],yaxis[:],'k',label = legend1 ) # still reading in test_data
	ax1.plot(xaxis,yaxis,'k',label = legend1 )

# Error surrounds data
#	plt.fill_between(xaxis[:],err_p[:],err_n[:],alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')
	plt.fill_between(xaxis,err_p,err_n,alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')

# Remove legend box frame 
	l = legend()
	l.draw_frame(False)
	draw()

#Set the visual range. Automatic range is ugly. 
	xmin = int(float(xaxis[0]))
	xmax = int(float(xaxis[-1]))
	plt.xlim((xmin,xmax))

#Label the figure and show
	plt.title( title )
	plt.xlabel( xlabel )
	plt.ylabel( ylabel )
	plt.savefig( image_title )

	plt.show()
	

