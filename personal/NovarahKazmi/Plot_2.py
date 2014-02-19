#import os
#import glob
import numpy as np
import matplotlib.pyplot as plt
#import scipy.interpolate as intp
import math
from pylab import *

"""
### Tasks
[X] what figures to plot - Just main light curve plot
[X] flexible data label 
[]  need to account for spaces in labels
[] zoom
[] stacked data (different y axis) 
"""

# open .txt file that has
# file location , title for plot, x axis label, y axis label, legend label
name = np.loadtxt('Labels_for_plot.txt',dtype = str)
test_data = np.loadtxt( name[0] ,dtype = str)

# Legend parameters, must come before plotting
params = {'legend.fontsize': 15, 'legend.linewidth':2}
plt.rcParams.update(params)

# Create the plot
ax1 = subplot(1,1,1)

# Don't like manual labels # ax1.plot(test_data[:,0],test_data[:,1],'k',label = leg_1 )
ax1.plot(test_data[:,0],test_data[:,1],'k',label = name[4] )

# FIX THIS : Need error to surround curve
ax1.fill_between(test_data[:,0],test_data[:,2],test_data[:,2],alpha=0.5, edgecolor='#1B2ACC', facecolor='#089FFF')

# Not sure how to standardize curve labeling when all cuves will be different and have different labels
#ax1.text(6000, 1000000, r"Plot_1", fontsize=20, color="k")

# Grid lines parameters
ax1.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)

# Remove legend box frame 
l = legend()
l.draw_frame(False)
draw()

#Label the figure
plt.title( name[1] )
plt.xlabel( name[2] )
plt.ylabel( name[3] )

# Label the saved plot

plt.savefig('title')

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
############## Below this line is garbage ###############
######### Final code for class won't have this ##########
# Not needed for what the code is doing at the moment
#f = figure()
#subplots_adjust(hspace=0.001)

# Allow customization of plots 
"""
# pauses created for ease
# Better to have users make a text file and write in the titles they want.
# This method is very time consuming 
title   = raw_input('Input title name: ' )
raw_input()
print "Symbol for Angstrom is: $\AA$ ( In case you needed ) "
x_label = raw_input('Label X Axis: ' )
raw_input()
y_label = raw_input('Label Y Axis: ' )
raw_input()
leg_1 = raw_input('Label for legend ')
"""

""" Dummy names
plt.title( "Hat" )
plt.xlabel( "Floor" )
plt.ylabel( "Lamp" )
"""
""" Don't like manual labels
plt.title( title )
plt.xlabel( x_label )
plt.ylabel( y_label )
plt.savefig('Plot_Test_2.png')
"""

