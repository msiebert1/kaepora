# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 14:53:37 2014

@author: QuantumMonkey
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# Initial condition setup
h = [8,2]
gs = gridspec.GridSpec(2, 1, height_ratios = h, hspace = 0.001)
font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 10,
        }
        
        
#Test Data        
xaxis_1 = np.linspace(0,10000,10000)
xaxis_2 = np.linspace(0,10000,10000)
comp_data1 = 2*xaxis_1
comp_data2 = comp_data1**0.5
rms_data1 = xaxis_1/xaxis_1*1000
rms_data2 = rms_data1*2

print xaxis_1,comp_data1


# The composite spectra plotting
Rel_flux = plt.subplot(gs[0])
plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
plt.minorticks_on()
plt.xlim(4000,7000)
#plt.ylim()
plt.fill_between(xaxis_1, comp_data1+rms_data1, comp_data1-rms_data1, color = 'cyan')        
plt.plot(xaxis_1, comp_data1, label = "Composite1",color = 'b')
plt.plot(xaxis_1, comp_data1+rms_data1, label = "+ RMS1",color = 'g',ls = '-')
plt.plot(xaxis_1, comp_data1-rms_data1, label = "- RMS1",color = 'g',ls = '-')
plt.fill_between(xaxis_2, comp_data2+rms_data2, comp_data2-rms_data2, color = 'orange')        
plt.plot(xaxis_2, comp_data2, label = "Composite2",color = 'r')
plt.plot(xaxis_2, comp_data2+rms_data2, label = "+ RMS2",color = 'y')
plt.plot(xaxis_2, comp_data2-rms_data2, label = "- RMS2",color = 'y')
plt.legend(prop = {'family' : 'serif'})
RFxticklabels = Rel_flux.get_xticklabels() 


# Residual plotting

Resid = plt.subplot(gs[1], sharex = Rel_flux)
plt.xlim(4000,7000)
plt.ylim(0,3000)
plt.ylabel('Residuals', fontdict = font)
plt.plot(xaxis_1, rms_data1, label = "RMS of residuals1", ls = '-')
plt.plot(xaxis_2, rms_data2, label = "RMS of residuals2", ls = '-')
RSxticklabels = Resid.get_xticklabels() 

plt.show()
