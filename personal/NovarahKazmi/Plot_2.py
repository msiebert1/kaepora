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
[X] what figures to plot - atm just main light curve plot 
[X] flexible data label 
[X] need to account for spaces in labels
[X] pipe functions
[X] personalizing the plots from the data
[] zoom

[] Write using latex font or similar
[] stacked data (different y axis) 

"""
def main(plot_data_1, plot_data_2, title, image_title, xlabel, ylabel, legend_1,legend_2):
#def main(xaxis, yaxis,err_p, err_n, title, image_title, xlabel, ylabel, legend1):
# open .txt file. Has the following: 
# file location , title for plot, title for saved figure, x axis label, y axis label, legend label(s)
#	name = np.loadtxt( 'Labels_for_plot.txt' , dtype = str , delimiter="\n" )
#	test_data = np.loadtxt( name[0] ,dtype = str)

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
#	ax1.plot(xaxis[:],yaxis[:],'k',label = name[5] ) # (X_comp, Y_comp, 'color', curve name)
#	ax1.plot(xaxis[:],yaxis[:],'k',label = legend1 ) # still reading in test_data
	ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
	ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

# Error surrounds data
#	plt.fill_between(xaxis[:],err_p[:],err_n[:],alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')
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
	
	
### Only plotting a single figure
### Subplots can be introduced at any point
"""
# axes rect in relative 0,1 coords left, bottom, width, height.  Turn
# off xtick labels on all but the lower plot

f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(3,1,1)
ax1.plot(t,s1)
ax1.text(1.1, 0.2, r"Plot_1", fontsize=20, color="b")
yticks(arange(-0.9, 1.0, 0.4))
ylim(-1,1)

ax2 = subplot(3,1,2, sharex=ax1)
ax2.plot(t,s2)
ax2.text(1.1, 0.2, r"Plot_2", fontsize=20, color="k")
yticks(arange(0.1, 1.0,  0.2))
ylim(0,1)

ax3 = subplot(3,1,3, sharex=ax1)
ax3.plot(t,s3)
ax3.text(1.1, 0.2, r"Plot_3", fontsize=20, color="b")
yticks(arange(-0.9, 1.0, 0.4))
ylim(-1,1)


xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
setp(xticklabels, visible=False)
"""

