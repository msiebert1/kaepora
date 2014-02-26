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
## Place the following after you have calculated your data 
## fill in information
# 
# plot_data_1  = [wavelength_1, Spectra_1, plus, minus]     # Expects [X, Y, Pos_Y Err, Neg_Y Err]
# plot_data_2  = [wavelength_2, Spectra_2, plus, minus]     # Expects [X, Y, Pos_Y Err, Neg_Y Err]
# title        = "Super Cool"                               # Within the quotes place plot title
# image_title  = "../personal/YOUR_REP/WHOA.png"            # Name the image
# xlabel       = "Wavelength ($\AA$)"                       # Label the X axis
# ylabel       = "Relative Flux"			    # Label the Y axis
# legend_1     = "First Spectra"                            # Name the title to see in the legend
# legend_2     = "Second Spectra"                           # Name second data array
#
## The following line will plot the data
#
# Plotting.main(plot_data,title, image_title , xlabel, ylabel, legend_1,legend_2)
#
#
###########################################
def main(plot_data_1, plot_data_2, title, image_title, xlabel, ylabel, legend_1,legend_2):

# re-name variables
	xaxis_1 = plot_data[0] 
	yaxis_1 = plot_data[1]
	err_p_1 = plot_data[2] 
	err_n_1 = plot_data[3]
	xaxis_2 = plot_data[0] 
	yaxis_2 = plot_data[1]
	err_p_2 = plot_data[2] 
	err_n_2 = plot_data[3]

# Legend parameters, must come before plotting
	params = {'legend.fontsize': 15, 'legend.linewidth':2}
	plt.rcParams.update(params)

# Create the plot
	ax1 = subplot(1,1,1)

# Read legend labels from files
	ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 ) # (X_comp, Y_comp, 'color', curve name)
	ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

# Error surrounds data
	plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
	plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

# Remove legend box frame 
	l = legend()
	l.draw_frame(False)
	draw()

#Set the visual range. Automatic range is ugly. 
	xmin = int(float(xaxis_1[0]))
	xmax = int(float(xaxis_2[-1]))
	plt.xlim((xmin,xmax))

#Label the figure and show
	plt.title( title )
	plt.xlabel( xlabel )
	plt.ylabel( ylabel )
	plt.savefig( image_title )

	plt.show()

