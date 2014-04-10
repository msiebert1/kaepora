import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec
import scipy.optimize as optimize
import random # new import so that colors in fill_between are random


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
#Relative_Flux = [wavelengths, avg_flux_p]  # Want to plot a composite of multiple spectrum
#Residuals     = []
#Spectra_Bin   = []
#Age           = []
#Delta         = []
#Redshift      = [] 
#Show_Data     = [Relative_Flux,Residuals]
#image_title  = "WHOA.png"            # Name the image (with location)
#title        = "Image is called this" 
#
## Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
##                   0              1          2            3    4      5         6,                 7
#Plots = [] # the plots you want to create
#
## The following line will plot the data
#
#Plotting.main(Show_Data , Plots, image_title , title)
#
#
###########################################

def main(Show_Data , Plots , image_title , title):
#############################################################
# Set the height of each figure
#############################################################
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
# I want a function/for-loop/if-statement
# that turns each component of Show_Data (which has each 
# composite, residual...etc) to its own variable name that 
# can be iterated through
# This is so that each of the plotting functions (Composite,
# Residual, Spec.... etc) can be run multiple time and the 
# final figure has multiple data sets that can be shown 
# overlapping on the same figure. 
#############################################################
    # Use the length of each component of Show_Data to rename and fill arrays
    len_RF = len(Show_Data[:][0])
    len_SB = len(Show_Data[:][1])
    len_AG = len(Show_Data[:][2])
    len_DE = len(Show_Data[:][3])
    len_RD = len(Show_Data[:][4])

    # Even values are x. Odd are y  (Slightly confusing for the time being)
    RF = []
    SB = []
    AG = []
    DE = []
    RD = []
    
    # Fill each array with data that will go in each plot
    for i in range(len_RF):
        data = Show_Data[:][0][i] 
        RF.append(data) 
    for i in range(len_SB):
        data = Show_Data[:][1][i] 
        SB.append(data)
    for i in range(len_AG):
        data = Show_Data[:][2][i] 
        AG.append(data) 
    for i in range(len_DE):
        data = Show_Data[:][3][i] 
        DE.append(data)  
    for i in range(len_RD):
        data = Show_Data[:][4][i] 
        RD.append(data)  
    
    names_1   = "Random String"

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
# smoothGauss Function: 
#############################################################
   # Error with SmoothGauss, it turns small data sets to none.. Not sure what the issue is
    deg = 5
    def smoothGauss(list, strippedXs = False, degree = deg):
        window = degree*2-1
        weight = np.array([1.0]*window)  
        weightGauss = []  

        for m in range(window):  
            m = m-degree+1  
            frac = m/float(window)  
            gauss = 1/(np.exp((4*(frac))**2))  
            weightGauss.append(gauss)  
        
        weight = np.array(weightGauss)*weight  
        smooth = [0.0]*(len(list)-window)  

        for m in range(len(smooth)):  
            smooth[m] = sum(np.array(list[m:m+window])*weight)/sum(weight)  

        return smooth     
#############################################################
# Scaling:
#############################################################
    def Scaling(data):
        scaled = []
        for m in range(len(data)):
            scaled.append((data[m]-min(data[m]))/(max(data[m])-min(data[m])))
        return scaled 
#############################################################
# residual: Takes the data (Y values?) and subracts it from
# the composite
#############################################################
    def residual(comp, data):
        return comp - data
#############################################################
# shift:
#############################################################
    def shift(key, array):
        a = []
        a = array
        return a[-key:]+a[:-key]
#############################################################
# The following function take the x,y,rms, and composite data    
# to plot the composite spectrum
#############################################################    
    def Composite(RF,comp_data,rms_data):
        # Rel_flux = plt.subplot(gs[p]) is turned into a global variable outside all the functions
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.minorticks_on()

        # This for loop plots all the relative fluxes that were read in
        # the iterating is not 
        for k in range(len_RF):
            if k % 2 == 0:
                plt.fill_between(RF[k], comp_data[k+1] + rms_data[k+1], comp_data[k+1] - rms_data[k+1], facecolor = random.choice(['g', 'r', 'c', 'm', 'y', 'k']),alpha=0.5)          
        plt.plot(RF[0], comp_data[1], label = "Composite")
        plt.plot(RF[0], comp_data[1] + rms_data[1], label = "+ RMS")
        plt.plot(RF[0], comp_data[1] - rms_data[1], label = "- RMS")
        plt.legend(prop = {'family' : 'serif'})
        RFxticklabels = Rel_flux.get_xticklabels()  
        if max(stacked) == 0:
            plt.setp(RFxticklabels, visible=True)
        else:
            plt.setp(RFxticklabels, visible=False)
#############################################################
# The following function uses the x and rms to plot the RMS 
#############################################################
    def Residual(RF, rms_data):
        Resid = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Residuals', fontdict = font)
        plt.plot(RF[0], rms_data[1], label = "RMS of residuals", ls = '-')
        RSxticklabels = Resid.get_xticklabels()
        if max(stacked) == 1:
            plt.setp(RSxticklabels, visible=True)
        else:
            plt.setp(RSxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows 
# the Spectra per bin with respect to wavelength.
#############################################################            
    def SpecBin(RF, rms_data):
        SpecBin = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Spectra/Bin', fontdict = font)
        plt.plot(RF[0], 0*RF[0]+3.0, label = "title goes here", ls = '-')
        SBxticklabels = SpecBin.get_xticklabels()        
        if max(stacked) == 2:
            plt.setp(SBxticklabels, visible=True)
        else:
            plt.setp(SBxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the age(?) with respect to wavelength(?).
#############################################################              
    def Age(RF, rms_data):
        Age = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Age [d]', fontdict = font)
        plt.plot(RF[0], RF[0]**3.0, label = "title goes here", ls = '-')
        AGxticklabels = Age.get_xticklabels()        
        if max(stacked) == 3:
            plt.setp(AGxticklabels, visible=True)
        else:
            plt.setp(AGxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the delta(?) with respect to wavelength(?).
#############################################################              
    def Delta(RF, rms_data):
        Delta = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('$\Delta$', fontdict = font)
        plt.plot(RF[0], 1/RF[0], label = "title goes here", ls = '-')
        DLxticklabels = Delta.get_xticklabels()
        if max(stacked) == 4:            
            plt.setp(DLxticklabels, visible=True)
        else:
            plt.setp(DLxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the redshift with respect to wavelength(?).
#############################################################              
    def Redshift(RF, rms_data):
        Redshift = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Redshift', fontdict = font)
        plt.plot(RF[0], RF[0]**2.0-RF[0], label = "title goes here", ls = '-')
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
    def Multi(xtrunc,RF):
        plt.figure(num = 2, dpi = 100, figsize = [8, 8], facecolor = 'w')
        for m in range(len(RF[1])):
            plt.plot(xtrunc, RF[1][m], label = str(names_1[m]))
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
       plt.annotate(str(names_1[0]), xy = (max(RF[0]), max(RF[1][0])), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       buffer = (max(RF[1][0])-min(RF[1][0]))/2.0
       for m in range(len(RF[1])-1):
            plt.plot(RF[0], RF[1][m+1]+(m+1)*(1+buffer))
            plt.annotate(str(names_1[m+1]), xy = (max(RF[0]), max(RF[1][m+1])+(m+1)*(1+buffer)), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
       plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
       plt.minorticks_on()
       plt.savefig('stacked plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
       
#############################################################
########################      End of function section
########################              :)   
########################  Proceed to function application
############################################################# 
# What does this section do? I was having problems with it
# so I've temporarily commented it out. 
#############################################################
    xtrunc = RF[0][deg-1:len(RF[0])-deg]
    """    
    RF[0] = []
    RF[0] = xtrunc  
    
    smoothed = []
    for m in range(len(RF[1])):
        smoothed.append(smoothGauss(RF[1][m]))
    
    RF[1] = []
    RF[1] = Scaling(smoothed)
    
    source = (np.array(RF[1])).T
    guess = RF[1][0]
    comp_data = []
    sbins = []
    print source
    for m in range(len(source)):
        comp, cov = optimize.leastsq(residual, guess[m], args = (source[m]), full_output = False)
        print "This is cov = ", cov
        print comp
        comp_data.append(comp[0])
        sbins.append(len(source[m]))
        
    delta_data = []     # difference from the mean value
    
    for m in range(len(RF[1])):   # finds difference between composite data and interpolated data sets
        delta_data.append(comp_data-RF[1][m]) 
        rms_data = np.sqrt(np.mean(np.square(delta_data), axis = 0))    # finds root-mean-square of differences at each wavelength within overlap
    """
    
############################################################# 
# Delete this later. I was having trouble with smoothing and
# needed data to play with. 
#############################################################
    comp_data = []
    rms_data  = []
    
    for j in range(len_RF):
        c = RF[j]*3
        r = RF[j]*1
        #if j % 2 != 0:
        comp_data.append(c) 
        rms_data.append(r)
    
    """    
    print comp_data[2]   
    
    comp_data[1] = RF[1]*3
    rms_data_1 = RF[1]*1
    # RMS and comp should easier to iterate so that multiple runs can be done
    
    comp_data_2 = []
    comp_data_2 = RF[3]*2
    rms_data_2 = []
    rms_data_2 = RF[3]*3
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
    plt.title(title)
    gs = gridspec.GridSpec(len(Plots), 1, height_ratios = h, hspace = 0.001)
    p = 0
    
    Rel_flux = plt.subplot(gs[p])

#############################################################
# The following series of if statments runs a specific portion 
# of the ploting functions. The p = p+1 allows the code to  
# iterate through all the possible plot options      
# if 0 or 1 or 2 or 3 or 4 or 5 in Plots.
############################################################# 
    
    if 0 in Plots: # will always plot a composite spectrum if any 1-5 are selected   
        Composite(RF, comp_data,rms_data) 
        p = p+1
    if 1 in Plots:
        Residual(RF, rms_data) 
        p = p+1
    if 2 in Plots:  
        SpecBin(RF, rms_data) # Currently set to RF, needs to be changed 
        p = p+1               # to it's associated array. 
    if 3 in Plots:
        Age(RF, rms_data)
        p = p+1
    if 4 in Plots:
        Delta(RF, rms_data)
        p = p+1
    if 5 in Plots:
        Redshift(RF, rms_data)
        p = p+1
    # Regardless of what is plotted, we label the Xaxis and save the plot image       
    plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
    plt.savefig('image_title', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
    #plt.savefig('composite plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
    
#############################################################
# Other figures to plot, we put these after the first figure
# is saved. So as to not confuse the code.     
#############################################################      
    if 6 in Plots:
        Multi(xtrunc,RF[1])
           
    if 7 in Plots:
        Stacked(RF[0],RF[1])
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
    #plt.show()


