#import os
#import glob
import numpy as np
import matplotlib.pyplot as plt
#import scipy.interpolate as intp
import math

"""
### Tasks
- what figures to plot
- flexible data label 
- zoom
- stacked data (different y axis) 
"""
from pylab import *

params = {'legend.fontsize': 10, 'legend.linewidth':2}
plt.rcParams.update(params)

t = arange(0.0, 2.0, 0.01)

s1 = sin(2*pi*t)
s2 = exp(-t)
s3 = s1*s2

f = figure()
subplots_adjust(hspace=0.001)

ax1 = subplot(1,1,1)
ax1.plot(t,s1,label = "Sin")
ax1.text(1.1, 0.2, r"Plot_1", fontsize=20, color="b")
yticks(arange(-0.9, 1.0, 0.4))
ylim(-1,1)
l = legend()
l.draw_frame(False)
draw()
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
show()

