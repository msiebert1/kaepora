import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter

"""
### Tasks
[X] what figures to plot 
[X] atm just main light curve plot 
[X] flexible data label 
[X] need to account for spaces in labels
[X] pipe functions
[X] personalizing the plots from the data
[] histogram
[] zoom

[] Write using latex font or similar
[] stacked data (different y axis) 

"""
def main(fig_type,plot_data_1,plot_data_2,plot_data_3,image_title,plot_labels,legend_names):
# re-name variables
	xaxis_1   = plot_data_1[0] 
	yaxis_1   = plot_data_1[1]
	err_p_1   = plot_data_1[2] 
	err_n_1   = plot_data_1[3]

	xaxis_2   = plot_data_2[0] 
	yaxis_2   = plot_data_2[1]
	err_p_2   = plot_data_2[2] 
	err_n_2   = plot_data_2[3]

	xaxis_3   = plot_data_3[0] 
	yaxis_3_1 = plot_data_3[1]
	yaxis_3_2 = plot_data_3[2] 

	title     = plot_labels[0]
	xlabel    = plot_labels[1]
	ylabel_1  = plot_labels[2]
	ylabel_2  = plot_labels[3]

	legend_1  = legend_names[0]
	legend_2  = legend_names[1]
	legend_3  = legend_names[2]
	legend_4  = legend_names[3]
	
# Legend parameters, must come before plotting
	params = {'legend.fontsize': 15, 'legend.linewidth':2}
	plt.rcParams.update(params)

	f = figure()
	subplots_adjust(hspace=0.001)


	if fig_type == 2.3: # Plotting details for two figures with a total of three sets of data
		ax1 = subplot(2,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

		# Plotting details for two figures with one set of data
		ax2 = subplot(2,1,2, sharex=ax1)
		plt.ylabel( ylabel_2 )
		ax2.plot(xaxis_3,yaxis_3_1,'b',label = legend_3 )
		ax2.plot(xaxis_3,yaxis_3_2,'g',label = legend_4 )
		
#Removes x-axis labels for ax1. May need later. 	
	# 	Commented code is to remove labels for both figures
	#	xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
	#	setp(xticklabels, visible=False)
		setp(ax1.get_xticklabels(), visible = False)
	
	elif fig_type == 1.2: # Plotting details for single figure with two sets of data
		ax1 = subplot(1,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

	else : # Plotting single plot with single set of data
		ax1 = subplot(1,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'k',label = legend_1 )
		plt.fill_between(xaxis_1,err_p_1,err_n_2,alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')


# Remove legend box frame 
	l = legend()
	l.draw_frame(False)
	draw()

# Should change this so that its not always the first set of data -> only works for Carbon data
#Set the visual range. Automatic range is ugly. 
	xmin = int(float(xaxis_1[0]))
	xmax = int(float(xaxis_2[-1]))
	plt.xlim((xmin,xmax))

#Label the figure and show
	plt.xlabel( xlabel )
	plt.savefig( image_title )
	plt.show()


"""	if num_of_plots == 1:
		ax1 = subplot(1,1,1) # Create the plot

		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 ) # Read legend labels from files
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

# Error surrounds data
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
"""
#	if num_of_plots == 2:

