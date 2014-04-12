import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec
import scipy.optimize as optimize
import random # new import so that colors in fill_between are random

# Eventually do scaling by choice
# naming
# comp_data
# rms_data

###########################################
#             HOW TO USE
## At the top of YOUR code write the following
# import Plotting
# 
## Place the following after you have calculated your data 
## fill in information
# 
#import Plotting
#
#Relative_Flux = [wavelengths_1, composite_1, wavelength_2, composite_2]  # Want to plot a composite of multiple spectrum
#Residuals     = []
#Spectra_Bin   = []
#Age           = []
#Delta         = []
#Redshift      = [] 
#Show_Data     = [Relative_Flux, Residuals, Spectra_Bin, Age , Delta , Redshift]
#image_title  = "WHOA.png"            # Name the image (with location)
#title        = "Title on the figure" 
#
## Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
##                   0              1          2            3    4      5         6,                 7
#Plots = [0] # the plots you want to create
#
## The following line will plot the data
#
#Plotting.main(Show_Data , Plots , image_title , title, Names)
#
#
###########################################

def main(Show_Data , Plots , image_title , title, Names):
#############################################################
# Set the height of each figure
#############################################################
    print "Begin plotting..."
    # Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
    #                   0              1          2            3    4      5         6,                 7
    Height =           [8,             2,         3,           2,   2,     2,        0,                 0]

    
    h = []

    for m in Plots:
        h.append(Height[m])

    stacked = []

    for m in arange(0,6,1):
        if m in Plots:
            stacked.append(m)
        continue
#############################################################
# Rename the variable.
#############################################################
    # Use the length of each component of Show_Data to rename and fill arrays
    len_RF = len(Show_Data[:][0])
    len_RS = len(Show_Data[:][1])
    len_SB = len(Show_Data[:][2])
    len_AG = len(Show_Data[:][3])
    len_DE = len(Show_Data[:][4])
    len_RD = len(Show_Data[:][5])

    # Even values are x. Odd are y  (Slightly confusing for the time being)
    RF = []
    RS = []
    SB = []
    AG = []
    DE = []
    RD = []

    #rms_data  = []
    #Names     = []
    # Fill each array with data that will go in each plot
    for i in range(len_RF):
        rf = Show_Data[:][0][i].T 
        RF.append(rf) 
    for i in range(len_RS):
        rs = Show_Data[:][0][i]
        RS.append(rs) 
    for i in range(len_SB):
        sb = Show_Data[:][1][i] 
        SB.append(sb)
    for i in range(len_AG):
        ag = Show_Data[:][2][i] 
        AG.append(ag) 
    for i in range(len_DE):
        de = Show_Data[:][3][i] 
        DE.append(de)  
    for i in range(len_RD):
        rd = Show_Data[:][4][i] 
        RD.append(rd) 
    
    """
    len_names = len(Names)
    if len(Names) == 0:
        Names.append("Spectrum")
    print Names
    """
#############################################################
# Changing font parameters
#############################################################
    params = {'legend.fontsize': 8, 
              'legend.linewidth': 2,
              'legend.font': 'serif',
              'mathtext.default': 'regular', 
              'xtick.labelsize': 8, 
              'ytick.labelsize': 8} # changes font size in the plot legend

    plt.rcParams.update(params)                             # reset the plot parameters

    font = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 10,
            } 

#############################################################
# Scaling:
#############################################################
    def Scaling(data):
        scaled = []
        #if len_RF == 2:
        #print "Min Data ", min(data)
        #print "Max Data ", max(data)
        #scaled = data
        scaled = ((data-min(data))/(max(data)-min(data)))
        #if len_RF > 2:
        #    for m in range(len(data)):
        #        scaled.append((data[m]-min(data[m]))/(max(data[m])-min(data[m])))
        return scaled 
            
        
#############################################################
# residual: Takes the data (Y values?) and subracts it from
# the composite
#############################################################
    """
    def residual(comp, data):
        return comp - data
    """
#############################################################
# The following function take the x,y,rms, and composite data    
# to plot the composite spectrum
#############################################################    
    def Composite(RF,RS,Names):
        # Rel_flux = plt.subplot(gs[p]) is turned into a global variable outside all the functions
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.minorticks_on()

        # This for loop plots all the relative fluxes that were read in
        # the iterating is not 
        """       
        for k in range(len_RF):
            if k % 2 == 0:
                plt.fill_between(RF[k], RF[k+1] + RS[k+1], RF[k+1] - RS[k+1], facecolor = random.choice(['g', 'r', 'c', 'm', 'y', 'k']),alpha=0.5)          
        """
        """
        if len_RF == 2:
            plt.plot(RF[0], RF[1], color ='g') 
        if len_RF == 4:
            plt.plot(RF[0], RF[1], color ='g', label = Names[0]) 
            plt.plot(RF[2], RF[3], color ='k', label = Names[1]) 
        if len_RF == 6:
            plt.plot(RF[0], RF[1], color ='g', label = Names[0]) 
            plt.plot(RF[2], RF[3], color ='k', label = Names[1])
            plt.plot(RF[4], RF[5], color ='b', label = Names[2]) 
        """
        
        for k in range(len_RF):
            if k % 2 == 0:
                plt.plot(RF[k], RF[k+1], color = random.choice(['g', 'r', 'c', 'm', 'y', 'k']), label = Names)
                #plt.plot(RF[k], RF[k+1] + RS[1], label = "+ RMS")
                #plt.plot(RF[k], RF[k+1] - RS[1], label = "- RMS")
       
        plt.legend(prop = {'family' : 'serif'})
        """        
        for n in range(len_names):
            legend([Names[n]])
        """   
        RFxticklabels = Rel_flux.get_xticklabels()  
        if max(stacked) == 0:
            plt.setp(RFxticklabels, visible=True)
        else:
            plt.setp(RFxticklabels, visible=False)
#############################################################
# The following function uses the x and rms to plot the RMS 
#############################################################
    def Residual(RS):
        Resid = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Residuals', fontdict = font)
        plt.plot(RS[0], RS[1], label = "RMS of residuals", ls = '-')
        RSxticklabels = Resid.get_xticklabels()
        if max(stacked) == 1:
            plt.setp(RSxticklabels, visible=True)
        else:
            plt.setp(RSxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows 
# the Spectra per bin with respect to wavelength.
#############################################################            
    def SpecBin(SB):
        SpecBin = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Spectra/Bin', fontdict = font)
        plt.plot(SB[0], SB[1], label = "title goes here", ls = '-')
        SBxticklabels = SpecBin.get_xticklabels()        
        if max(stacked) == 2:
            plt.setp(SBxticklabels, visible=True)
        else:
            plt.setp(SBxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the age(?) with respect to wavelength(?).
#############################################################              
    def Age(AG):
        Age = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Age [d]', fontdict = font)
        plt.plot(AG[0], AG[1], label = "title goes here", ls = '-')
        AGxticklabels = Age.get_xticklabels()        
        if max(stacked) == 3:
            plt.setp(AGxticklabels, visible=True)
        else:
            plt.setp(AGxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the delta(?) with respect to wavelength(?).
#############################################################              
    def Delta(DE):
        Delta = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('$\Delta$', fontdict = font)
        plt.plot(DE[0], DE[1], label = "title goes here", ls = '-')
        DLxticklabels = Delta.get_xticklabels()
        if max(stacked) == 4:            
            plt.setp(DLxticklabels, visible=True)
        else:
            plt.setp(DLxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the redshift with respect to wavelength(?).
#############################################################              
    def Redshift(RE):
        Redshift = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Redshift', fontdict = font)
        plt.plot(RE[0], RE[1], label = "title goes here", ls = '-')
        Zxticklabels = Redshift.get_xticklabels()        
        if max(stacked) == 5:            
            plt.setp(Zxticklabels, visible=True)
        else:
            plt.setp(Zxticklabels, visible=False)
#############################################################
# The following function will take all the data sets to build  
# the composite spectrum and lay each one over the next. 
# This will appear on a separate figure.
############################################################# 
    # Temporary comment until we get data for multiple arrays
    """
    def Multi(RF,rms_data):
        plt.figure(num = 2, dpi = 100, figsize = [8, 8], facecolor = 'w')
        for m in range(len(RF[1])):
            plt.plot(xtrunc, RF[1][m], label = str(Names[m]))
        plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.minorticks_on()
        plt.legend(prop = {'family' : 'serif'})        
        plt.savefig('multi-spectrum plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
#############################################################
# The following function stacks all the data sets used to  
# build the composite spectrum. This will appear on a 
# separate figure.
#############################################################  
    def Stacked(RF,rms_data): 
       plt.figure(num = 3, dpi = 100, figsize = [8, 4*len(RF[1])], facecolor = 'w')
       plt.plot(RF[0], RF[1][0])
       plt.annotate(str(Names[0]), xy = (max(RF[0]), max(RF[1][0])), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       buffer = (max(RF[1][0])-min(RF[1][0]))/2.0
       for m in range(len(RF[1])-1):
            plt.plot(RF[0], RF[1][m+1]+(m+1)*(1+buffer))
            plt.annotate(str(Names[m+1]), xy = (max(RF[0]), max(RF[1][m+1])+(m+1)*(1+buffer)), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
       plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
       plt.minorticks_on()
       plt.savefig('stacked plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
       """
#############################################################
########################      End of function section
########################              :)   
########################  Proceed to function application
############################################################# 
# What does this section do? I was having problems with it
# so I've temporarily commented it out. 
#############################################################
    
    for j in range(len_RF):
        if j % 2 != 0:
            RF[j] = Scaling(RF[j])
    # Commented out for the time being because we're only given a single array
    """
    source = (np.array(RF[1])).T

    sbins = []
    guess = []
    print source
    
    if len(source) !=1:
        for m in range(len(source)):
            comp, cov = optimize.leastsq(residual, guess[m], args = (source[m]), full_output = False)
            comp.append(comp[0])
            sbins.append(len(source[m]))
        else:
            comp_data = RF[1]
        
    delta_data = []     # difference from the mean value
    
    for m in range(len(RF[1])):   # finds difference between composite data and interpolated data sets
        delta_data.append(comp-RF[1][m]) 
        rms_data = np.sqrt(np.mean(np.square(delta_data), axis = 0))    # finds root-mean-square of differences at each wavelength within overlap
    """
#############################################################
# The following section sets up the plotting figure information. 
# it sets the size, color, title. 
# gs allows multiple plots to be shown on the same image
# with their own x axis but in the same figure.
# Rel_Flux is the name of the shared y axis, that will be 
# stacked on
# (not sure if I understand this 100%)
#############################################################
    plt.figure(num = 1, dpi = 100, figsize = [8, np.sum(h)], facecolor = 'w')
    gs = gridspec.GridSpec(len(Plots), 1, height_ratios = h, hspace = 0.001)
    
    p = 0
    Rel_flux = plt.subplot(gs[p])
    plt.title(title, fontdict = font)
#############################################################
# The following series of if statments runs a specific portion 
# of the ploting functions. The p = p+1 allows the code to  
# iterate through all the possible plot options      
# if 0 or 1 or 2 or 3 or 4 or 5 in Plots.
############################################################# 
    
    if 0 in Plots: # will always plot a composite spectrum if any 1-5 are selected   
        Composite(RF,RS,Names) 
        p = p+1
    if 1 in Plots:
        Residual(RS) 
        p = p+1
    if 2 in Plots:  
        SpecBin(SB)  
        p = p+1               
    if 3 in Plots:
        Age(AG)
        p = p+1
    if 4 in Plots:
        Delta(DE)
        p = p+1
    if 5 in Plots:
        Redshift(RD)
        p = p+1
    # Regardless of what is plotted, we label the Xaxis and save the plot image  
    
    plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
    
    plt.savefig(image_title, dpi = 600, facecolor='w', edgecolor='w', pad_inches = 0.1)
    print "Plotting complete..."    
    #plt.savefig('composite plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
    
#############################################################
# Other figures to plot, we put these after the first figure
# is saved. So as to not confuse the code.     
#############################################################      
    """
   if 6 in Plots:
        Multi(RF[0],RF[1])
           
    if 7 in Plots:
        Stacked(RF[0],RF[1])
    """
#############################################################
# Last steps are removing the legend frame *crosses fingers*
#############################################################                  
        
# Remove legend box frame 
    #l = legend()
    #l.draw_frame(False)
    #plt.draw()

#Set the visual range. Automatic range is ugly. 
#    xmin = int(float(RF[0][0])) 
#    xmax = int(float(RF[0][-1])) 
#    plt.xlim((xmin,xmax))

#Label the figure and show
    #plt.xlabel( "Wavelength $\AA$" )
    #plt.savefig( image_title )
    plt.show()

