# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:43:46 2014

@author: bfry
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
#                   0              1          2            3    4      5
Height =           [8,             2,         3,           2,   2,     2]

Plots = [0,1,2,3,4,5] # Plots to generate

h = []

for m in Plots:
    h.append(Height[m])

params = {'legend.fontsize': 8, 
          'legend.linewidth': 2,
          'legend.font': 'serif',
          'xtick.labelsize': 8, 
          'ytick.labelsize': 8} # changes font size in the plot legend

plt.rcParams.update(params) # reset the plot parameters

font = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'bold',
        'size'   : 8,
        }

x_val = np.linspace(2000, 8000, 12001)

plt.figure(num = 2, dpi = 100, figsize = [8, np.sum(h)], facecolor = 'w')
gs = gridspec.GridSpec(len(Plots), 1, height_ratios = h, hspace = 0.001)

p = 0

if 0 in Plots:
    Rel_flux = plt.subplot(gs[p])
    #plt.title("".join(["$^{26}$Al / $^{27}$Al, t$_{res}$ = ", str(t_width_Al), " kyr", ", U$_{Al}$ = ", str(uptake_Al[0])]), fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('Relative f$_{\lambda}$', fontdict = font)
    #plt.axis([Start_Al-0.05, End_Al+0.05, 0, 40])
    plt.minorticks_on()
    #plt.xticks(np.arange(Start_Al, End_Al+0.05, x_tik))
    #plt.yticks(np.arange(0, 41, 5))
    plt.plot(x_val, x_val**2.0, label = "generic data", ls = '-')
    plt.legend(prop = {'family' : 'serif'})
    p = p+1

if 1 in Plots:
    Resid = plt.subplot(gs[p], sharex = Rel_flux)
    #plt.title("".join(["$^{53}$Mn / $^{55}$Mn, t$_{res}$ = ", str(t_width_Mn), " kyr", ", U$_{Mn}$ = ", str(uptake_Mn[0])]), fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('Residuals', fontdict = font)
    #plt.axis([Start_Mn-0.05, End_Mn+0.05, 0, 10])
    #plt.minorticks_on()
    #plt.xticks(np.arange(Start_Mn, End_Mn+0.05, x_tik))
    #plt.yticks(np.arange(0, 11, 1))
    plt.plot(x_val, 2.0*x_val, label = "title goes here", ls = '-')
    #plt.legend(prop = {'family' : 'serif'})
    p = p+1

if 2 in Plots:
    SpecBin = plt.subplot(gs[p], sharex = Rel_flux)
    #plt.title("".join(["$^{60}$Fe / $^{56}$Fe, t$_{res}$ = ", str(t_width_Fe), " kyr", ", U$_{Fe}$ = ", str(uptake_Fe[0])]), fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('Spectra/Bin', fontdict = font)
    #plt.axis([Start_Fe-0.05, End_Fe+0.05, 0, 40])
    #plt.minorticks_on()
    #plt.xticks(np.arange(Start_Fe, End_Fe+0.05, x_tik))
    #plt.yticks(np.arange(0, 41, 5))
    plt.plot(x_val, 0*x_val+3.0, label = "title goes here", ls = '-')
    #plt.legend(prop = {'family' : 'serif'})
    p = p+1

if 3 in Plots:
    Age = plt.subplot(gs[p], sharex = Rel_flux)
    #plt.title('Decay corrected $^{26}$Al / $^{27}$Al', fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('Age [d]', fontdict = font)
    #plt.axis([Start_Al-0.05, End_Al+0.05, 0, 30])
    #plt.minorticks_on()
    #plt.xticks(np.arange(Start_Al, End_Al+.05, x_tik))
    #plt.yticks(np.arange(0, 29, 5))
    plt.plot(x_val, x_val**3.0, label = "title goes here", ls = '-')
    #plt.legend(prop = {'family' : 'serif'})
    p = p+1

if 4 in Plots:
    Delta = plt.subplot(gs[p], sharex = Rel_flux)
    #plt.title('Decay corrected $^{53}$Mn / $^{55}$Mn', fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('$\Delta$', fontdict = font)
    #plt.axis([Start_Mn-0.05, End_Mn+0.05, 0, 13])
    #plt.minorticks_on()
    #plt.xticks(np.arange(Start_Mn, End_Mn+0.05, x_tik))
    #plt.yticks(np.arange(0, 13, 2))
    plt.plot(x_val, 1/x_val, label = "title goes here", ls = '-')
    #plt.legend(prop = {'family' : 'serif'})
    p = p+1

if 5 in Plots:
    Redshift = plt.subplot(gs[p], sharex = Rel_flux)
    #plt.title('Decay corrected $^{60}$Fe / $^{56}$Fe', fontdict = font)
    #plt.xlabel('Age [Myr]', fontdict = font)
    plt.ylabel('Redshift', fontdict = font)
    #plt.axis([Start_Fe-0.05, End_Fe+0.05, 0, 70])
    #plt.minorticks_on()
    #plt.xticks(np.arange(Start_Fe, End_Fe+0.05, x_tik))
    #plt.yticks(np.arange(0, 70, 10))
    plt.plot(x_val, x_val**2.0-x_val, label = "title goes here", ls = '-')
    #plt.legend(prop = {'family' : 'serif'})
    p = p+1

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

plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)

plt.show()


