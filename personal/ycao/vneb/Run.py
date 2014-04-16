#!/usr/bin/env python -i  

import composite
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
#import targeted
"""
Here's the main function that anyone can use to run some analysis.

The first line runs the composite program.
This grabs the selected data from SNe.db.
Change the input to composite.main to whatever you want your database query to be.
This is the code that outputs a text file of data for later use.
It also shows/saves a basic plot of the composite data and associated errors.

The next section sets up the plots to be used, and then plots using Plotting.py
Right now it doesn't work...but it will

More can be added as our programs evolve

Here are a few parameters you should set beforehand.
plot_name is where the plot showing both composite spectra together will be saved.
wmin and wmax define the boundary of the plot.
"""

def loadPara(parafile):
    f = open(parafile)
    lines = f.readlines()
    f.close()

    names = [] # SN name
    vnebs = [] # heliocentric redshift
    verrs = [] # MJD at B-band maximum light (99999.9=no estimate)


    for line in lines:
        p = line.split()
        if p[0][0] != '#':
            names.append(p[0])
            vnebs.append(float(p[1]))
            verrs.append(float(p[2]))
    return names,vnebs,verrs

def vsamples(vmin,vmax):
    allnames,allvnebs,allverrs = loadPara("nebular.dat")
    names = []
    vnebs = []
    for vneb,name in zip(allvnebs,allnames):
        if vneb >= vmin and vneb < vmax:
           names.append(name)
           vnebs.append(vneb)
    return names, vnebs

def vcomp(vmin,vmax):

    names,vnebs = vsamples(vmin,vmax)

    snnames = str(names)
    import string
    snnames = string.replace(snnames,'[','(')
    snnames = string.replace(snnames,']',')')


    selection = " ".join(('SELECT * FROM Supernovae WHERE',
                      'SN IN', snnames,
                      'AND Phase > 150 AND snr > 2.0',
                      'GROUP BY SN')) 

    print selection
    comp = composite.main(selection)
    return comp


composite1 = vcomp(-5000,0)
composite2 = vcomp(0,5000)

wmin = 3500
wmax = 8000

# #Read whatever you sasved the table as
# Data = Table.read(composite.savedname, format='ascii')


####################
#plotting

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
        
xaxis_1 = composite1.wavelength
xaxis_2 = composite2.wavelength
comp_data1 = composite1.flux
comp_data2 = composite2.flux
rms_data1 = composite1.ivar ** 0.5
rms_data2 = composite2.ivar ** 0.5

# The composite spectra plotting
Rel_flux = plt.subplot(gs[0])
plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
plt.minorticks_on()
plt.xlim(4000,8000)
plt.ylim(0, 1.0e-15)
#plt.ylim()
plt.fill_between(xaxis_1, comp_data1+rms_data1, comp_data1-rms_data1, color = 'cyan')        
plt.plot(xaxis_1, comp_data1, label = "Composite-Blue",color = 'b')
plt.plot(xaxis_1, comp_data1+rms_data1, label = "+ RMS-Blue",color = 'g',ls = '-')
plt.plot(xaxis_1, comp_data1-rms_data1, label = "- RMS-Blue",color = 'g',ls = '-')
plt.fill_between(xaxis_2, comp_data2+rms_data2, comp_data2-rms_data2, color = 'orange')        
plt.plot(xaxis_2, comp_data2, label = "Composite-Red",color = 'r')
plt.plot(xaxis_2, comp_data2+rms_data2, label = "+ RMS-Red",color = 'y')
plt.plot(xaxis_2, comp_data2-rms_data2, label = "- RMS-Red",color = 'y')
plt.legend(prop = {'family' : 'serif'})
RFxticklabels = Rel_flux.get_xticklabels() 


# Residual plotting

Resid = plt.subplot(gs[1], sharex = Rel_flux)
plt.xlim(4000,8000)
plt.ylim(0, 1.1e-34**0.5)
plt.ylabel('Residuals', fontdict = font)
plt.plot(xaxis_1, rms_data1, label = "RMS of residuals1", color = 'b',ls = '-')
plt.plot(xaxis_2, rms_data2, label = "RMS of residuals2", color = 'r',ls = '-')
plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
RSxticklabels = Resid.get_xticklabels()

plot_name = 'vnebula'
# plt.savefig('plots/' + plot_name + '.png')
plt.savefig(plot_name + '.png')
plt.show()


# # #To be honest, I'm not sure entirely how this works.
# # #Can someone who worked on this piece of code work with it?
# Relative_Flux = [composite.wavelength, composite.flux, composite.name]  # Want to plot a composite of multiple spectrum
# Residuals     = [composite.wavelength, composite.ivar]
# # Spectra_Bin   = []
# # Age           = []
# # Delta         = []
# # Redshift      = [] 
# Show_Data     = [Relative_Flux,Residuals]
# image_title  = "plots/"+ 'Composite' + str(vmin) + 'to' + str(vmax)  # Name the image (with location)
# title        = "Composite Spectrum" 

# # # Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
# # #                   0              1          2            3    4      5         6,                 7
# Plots = [0,1] # the plots you want to create

# # # The following line will plot the data
# # #It's commented out until it works...
# Plotting.main(Show_Data , Plots, image_title , title)
