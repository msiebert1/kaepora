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
[] flexible data label 
[] zoom
[] stacked data (different y axis) 
"""
# Read in test data 
test_data = np.loadtxt('../../data/cfa/sn2006E/sn2006E-20060126.53-fast.flm',dtype = str)

# Legend parameters, must come before plotting
params = {'legend.fontsize': 15, 'legend.linewidth':2}
plt.rcParams.update(params)

#f = figure()
#subplots_adjust(hspace=0.001)

ax1 = subplot(1,1,1)
ax1.plot(test_data[:,0],test_data[:,1],'k',label = "SN2006E")
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
plt.title('Make Title an input')
plt.xlabel('Wavelength (probably wont change)')
plt.ylabel('Make into an input statement')

plt.show()

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

