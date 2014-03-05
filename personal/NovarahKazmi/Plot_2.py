import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *

#import matplotlib.gridspec as gridspec # Doesn't exist on EWS Linux

#import matplotlib.font_manage
#from matplotlib.font_manager import FontProperties
#from matplotlib.ticker import FuncFormatter

# Show_Data = [Light_curve, Residual, Age, Spectra_Bin, Age, Delta, Redshift] 

#def main(fig_type,height,plot_data_1,plot_data_2,plot_data_3,image_title,plot_labels,legend_names):
def main(fig_type, Show_Data, Plots, num_data_set, image_title , plot_labels , legend)):
# Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
	#                    0              1          2            3    4      5
	#Height =           [8,             2,         3,           2,   2,     2]

	#Plots = [0,1,2,3,4,5] # Plots to generate

	#h = []

	#for m in Plots:
	#    h.append(Height[m])

# re-name variables

#	for i in range(len(Plots)):
#		xaxis[] = Show_Data[:][0][0] 

	
	if len(Show_Data[:][0]) > 0 :
		xaxis_1   = Show_Data[:][0][0] 
		yaxis_1   = Show_Data[:][0][1]
		err_p_1   = Show_Data[:][0][2] # Maybe change this to second x and y array
		err_n_1   = Show_Data[:][0][3] 
	if len(Show_Data[:][1]) > 0 :
		xaxis_2   = Show_Data[:][1][0]
		yaxis_2   = Show_Data[:][1][1]
		err_p_2   = Show_Data[:][1][2]
		err_n_2   = Show_Data[:][1][3]
	if len(Show_Data[:][2]) > 0:
		xaxis_3   = Show_Data[:][2][0] 
		yaxis_3   = Show_Data[:][2][1]
		xaxis_3   = Show_Data[:][2][2] # Will change this to accept potential second x-axis
		yaxis_3   = Show_Data[:][2][3]
	if len(Show_Data[:][3) > 0 :
		xaxis_4   = Show_Data[:][3][0]
		yaxis_4   = Show_Data[:][3][1]
		err_p_4   = Show_Data[:][3][2]
		err_n_4   = Show_Data[:][3][3]
	if len(Show_Data[:][4]) > 0 :
		xaxis_5   = Show_Data[:][4][0]
		yaxis_5   = Show_Data[:][4][1]
		err_p_5   = Show_Data[:][4][2]
		err_n_5   = Show_Data[:][4][3]
	if len(Show_Data[:][5]) > 0 :
		xaxis_6   = Show_Data[:][5][0]
		yaxis_6   = Show_Data[:][5][1]
		err_p_6   = Show_Data[:][5][2]
		err_n_6   = Show_Data[:][5][3]

# What should I do with legend

	title     = plot_labels[0]	
	xlabel    = plot_labels[1] # Dont need to label x and y
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
	#plt.figure(num = 2, dpi = 100, figsize = [8, np.sum(h)], facecolor = 'w')
	#gs = gridspec.GridSpec(len(Plots), 1, height_ratios = h, hspace = 0.001)

# Changing font parameters. 
	params = {'legend.fontsize': 8,  
		  'legend.linewidth': 2,
	          'legend.font': 'serif',
	          'xtick.labelsize': 8, 
	          'ytick.labelsize': 8} # changes font size in the plot legend

	plt.rcParams.update(params) # reset the plot parameters

	font = {'family' : 'serif',   # Font Parameters
		'color'  : 'black',
	        'weight' : 'bold',
	        'size'   : 8,
	        }

# Types of plots. This will be changed so that the options can be automated. 
	if fig_type == 2.3: # Plotting details for two figures with a total of three sets of data
		ax1 = subplot(2,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		#ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		#plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

		################## Dont need ###################################################
		if num_of_data[0] == 2:
			ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )
			plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
		################################################################################
		# Plotting details for two figures with one set of data
		ax2 = subplot(2,1,2, sharex=ax1)
		plt.ylabel( ylabel_2 )
		ax2.plot(xaxis_3,yaxis_3_1,'b',label = legend_3 )
		ax2.plot(xaxis_3_2,yaxis_3_2,'g',label = legend_4 )

	elif fig_type == 1.2: # Plotting details for single figure with two sets of data
		ax1 = subplot(1,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		#ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		#plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

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
#	wave_min = []
#	wave_max = []
#	for m in range(len(yaxis_1)):
#		wave_min.append(max(yaxis_1[m][1]))
#		wave_max.append(min(yaxis_1[m][1]))
	
	#xmin = min(wave_min)
	#xmax = max(wave_max)
	#print xmin, xmax

	xmin = int(float(xaxis_1[0]))
	xmax = int(float(xaxis_1[-1]))
	plt.xlim((xmin,xmax))

#Label the figure and show
	plt.xlabel( xlabel, fontdict = font )
	#plt.savefig( image_title )
	plt.show()

"""	

	if num_of_data[1] == 2:
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
	if num_of_data[2] == 2:
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
	if num_of_data[3] == 2:
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
	if num_of_data[4] == 2:
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
"""	

"""	xaxis_1   = Show_Data[:][0][0] 
	yaxis_1   = Show_Data[:][0][1]
	err_p_1   = Show_Data[:][0][2] 
	err_n_1   = Show_Data[:][0][3]

	xaxis_2   = Show_Data[:][1][0]
	yaxis_2   = Show_Data[:][1][1]
	err_p_2   = Show_Data[:][1][2]
	err_n_2   = Show_Data[:][1][3]

	xaxis_3   = Show_Data[:][2][0]
	yaxis_3_1 = Show_Data[:][2][1]
	yaxis_3_2 = Show_Data[:][2][2]
"""
"""	if num_of_plots == 1:
		ax1 = subplot(1,1,1) # Create the plot

		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 ) # Read legend labels from files
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

# Error surrounds data
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')
"""
#	if num_of_plots == 2:
"""
# Scale the flux. Temporary method, will try other methods. 
#for i in range(plot_data_1):
	yaxis_1   /= np.median(plot_data_1[:][1])
	yaxis_2   /= np.median(plot_data_2[:][1])
	yaxis_3_1 /= np.median(plot_data_3[:][1])
	yaxis_3_2 /= np.median(plot_data_3[:][2])
#	new_spectra.append([wavelength, yaxis])
	print yaxis_3_2
	# There is a group interpolating and de-reddening the data
	#Interpolating
	#interpolate(data_pos, fit_flux_pos)
	#interpolate(data_neg, fit_flux_neg)
	
	#avg_flux_p = sum(fit_flux_pos)/len(data_pos)
	#avg_flux_n = sum(fit_flux_neg)/len(data_neg)


	#residual and RMS 
	res_flux_p = []
	res_flux_n = []

	for i in range(len(data_pos)):
		res_flux_p.append(fit_flux_pos[i]-avg_flux_p)

	rms_p     = np.sqrt(np.mean(np.square(res_flux_p),axis = 0))
	plus_p    = avg_flux_p + rms_p
	minus_p   = avg_flux_p - rms_p
	scatter_p = np.divide(rms_p,avg_flux_p)


	for i in range(len(data_neg)):
		res_flux_n.append(fit_flux_neg[i]-avg_flux_n)

	rms_n     = np.sqrt(np.mean(np.square(res_flux_n),axis = 0))
	plus_n    = avg_flux_n + rms_n
	minus_n   = avg_flux_n - rms_n
	scatter_n = np.divide(rms_n,avg_flux_n)

		
# Removes x-axis labels for ax1. May need later. 	
	# 	Commented code is to remove labels for both figures
	#	xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
	#	setp(xticklabels, visible=False)
		setp(ax1.get_xticklabels(), visible = False)
	

"""
"""
Light_cur_x = [] 
	Light_cur_x = []

	xaxis_2   = [] 
	yaxis_2   = []

	xaxis_3   = [] 
	yaxis_3   = []

	xaxis_4   = [] 
	yaxis_4   = []

	xaxis_5   = [] 
	yaxis_5   = []

	for i in range(3):
		xaxis_i = Show_Data[:][i][0]
		yaxis_i = Show_Data[:][i][1]

		print xaxis_2
"""
"""
# Removes the xaxis labels
if Plots[len(Plots)-1] != 0:
    RFxticklabels = Rel_flux.get_xticklabels()
    plt.setp(RFxticklabels, visible=False)
if Plots[len(Plots)-1] != 1:
    RSxticklabels = Resid.get_xticklabels()
    plt.setp(RSxticklabels, visible=False)
if Plots[len(Plots)-1] != 2:
    SBxticklabels = SpecBin.get_xticklabels()
    plt.setp(SBxticklabels, visible=False)
if Plots[len(Plots)-1] != 3:
    AGxticklabels = Age.get_xticklabels()
    plt.setp(AGxticklabels, visible=False)
if Plots[len(Plots)-1] != 4:
    DLxticklabels = Delta.get_xticklabels()
    plt.setp(DLxticklabels, visible=False)
"""
