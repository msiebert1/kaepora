import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.font_manager
from matplotlib.font_manager import FontProperties

#Doesn't look professional yet. Neet better fonts. Test file is also childish looking. 
"""
### Tasks
[X] what figures to plot - atm just main light curve plot 
[X] flexible data label 
[X]  need to account for spaces in labels
[]  Write using latex font or similar
[] stacked data (different y axis) 
[] zoom
"""
"""
#Testing out different font commands
import tkFont
helv36 = tkFont.Font(family="Helvetica",size=36,weight="bold")
"""

# open .txt file. Has the following: 
# file location , title for plot, title for saved figure, x axis label, y axis label, legend label(s)
name = np.loadtxt( 'Labels_for_plot.txt' , dtype = str , delimiter="\n" )
test_data = np.loadtxt( name[0] ,dtype = str)

# This version of python doesn't like using LaTex
# 'usetex = True' causes the code to crash
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=False)

# re-name variables
xaxis = test_data[:,0] 
yaxis = test_data[:,1]
err_p = test_data[:,2] 
err_n = test_data[:,3]

# Legend parameters, must come before plotting
params = {'legend.fontsize': 15, 'legend.linewidth':2}
plt.rcParams.update(params)

# Create the plot
ax1 = subplot(1,1,1)

# Read legend labels from files
# (X_comp, Y_comp, 'color', curve name)
ax1.plot(xaxis[:],yaxis[:],'k',label = name[5] )

# Error surrounds data
plt.fill_between(xaxis[:],err_p[:],err_n[:],alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')

# Not sure how to standardize curve labeling when all cuves will be different and have coordinates
#ax1.text(6000, 1000000, r"Plot_1", fontsize=20, color="k")

# Grid lines parameters
#ax1.grid(color='b', alpha=0.5, linestyle='dashed', linewidth=0.5)

# Remove legend box frame 
l = legend()
l.draw_frame(False)
draw()

#Label the figure
plt.title( name[2] )
plt.xlabel( name[3] )
plt.ylabel( name[4] )

# Label the saved plot
plt.savefig( name[1] )

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

