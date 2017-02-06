import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec
#from matplotlib.figure import Figure
import scipy.optimize as optimize
import random # new import so that colors in fill_between are random
import operator
# Eventually do scaling by choice
# naming
# rms data

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
#image_title   = "WHOA.png"            # Name the image (with location)
#title         = "Title on the figure" 
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

def main(Show_Data , Plots , image_title , title , Names , xmin , xmax):
#############################################################
# Set the height of each figure
#############################################################
    print "Begin plotting..."
    # Available Plots:  Relative Flux,Variance, Residuals, Spectra/Bin, Age, Delta, Redshift, Stacked  
    #                   0              1          2            3         4      5     6         7
    Height =           [6,             2,         2,           2,        2,     2,    2,        0,      0,      0]
    
    h = []
    
    for m in Plots:
        h.append(Height[m])

    stacked = []

    for m in arange(0,7,1):
        if m in Plots:
            stacked.append(m)
        continue
#############################################################
# Rename the variable.
#############################################################
    # Use the length of each component of Show_Data to rename and fill arrays
    len_RF = len(Show_Data[:][0])    
    len_VA = len(Show_Data[:][1])
    len_RS = len(Show_Data[:][2])
    len_SB = len(Show_Data[:][3])
    len_AG = len(Show_Data[:][4])
    len_DE = len(Show_Data[:][5])
    len_RD = len(Show_Data[:][6])
    len_LC = len(Show_Data[:][7])
    len_UC = len(Show_Data[:][8])

    # Even values are x. Odd are y  (Slightly confusing for the time being)
    RF = []
    RS = []
    SB = []
    AG = []
    DE = []
    RD = []
    VA = []
    LC = []
    UC = []

    # Fill each array with data that will go in each plot
    for i in range(len_RF):
        rf = Show_Data[:][0][i].T 
        RF.append(rf)
    for i in range(len_VA):
        va = Show_Data[:][1][i].T
        VA.append(va) 
    for i in range(len_RS): 
        rs = Show_Data[:][2][i].T
        RS.append(rs)    
    for i in range(len_SB):
        sb = Show_Data[:][3][i].T 
        SB.append(sb)
    for i in range(len_AG):
        ag = Show_Data[:][4][i].T
        AG.append(ag) 
    for i in range(len_DE):
        de = Show_Data[:][5][i].T
        DE.append(de)  
    for i in range(len_RD):
        rd = Show_Data[:][6][i].T
        RD.append(rd) 
    for i in range(len_LC):
        lc = Show_Data[:][7][i].T
        LC.append(lc)
    for i in range(len_UC):
        uc = Show_Data[:][8][i].T
        UC.append(uc)
    
    print len(LC)
    """    
    def remove_zero(data):
        delete = []
        for i in range(data[0]):
            if (data[i] == 0):
                delete.append(i)
                print delete
                for j in range(len(delete)):
                    del data[0][delete[len(delete)-1-j]]
                    del data[1][delete[len(delete)-1-j]]
            
    remove_zero(RS)
    print RS
    
    print "Length of RS :", len_RS
    print "Length of VA :", len_VA
    print "Length of x1 :", len(RF[0])
    print "Length of y1 :", len(RF[1])
    print "Length of x2 :", len(VA[0])
    print "Length of y2 :", len(VA[1])
    """
    """
    len_names = len(Names)
    for n in range(len_names*2):
        Names.append("Spectrum")
    print Names
    """
#############################################################
# Changing font parameters
#############################################################
#    params = {'legend.fontsize': 10, 
#              'legend.linewidth': 2,
#              'legend.font': 'serif',
#              'mathtext.default': 'regular', 
#              'xtick.labelsize': 10, 
#              'ytick.labelsize': 10} # changes font size in the plot legend
#
#    plt.rcParams.update(params)                             # reset the plot parameters

    font = {'family' : 'serif',
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 10,
            } 

#############################################################
# Scaling: Normalize the y values. For reference the error 
# in scaling (I think) was caused by calculating rows 
# verses colums. So, I changed the file renaming and added .T
#############################################################
    def Scaling(data, medvalue):
        scaled = []
        #still getting an error here, which might be causing the 'nan'
        #both data.max() and data.min() are 0.0, which is a problem.
        #print max(data), min(data)
        #scaled = (data-min(data))/(max(data)-min(data))
        #scaled = (data)/(max(data))
        scaled = (data)/(medvalue)
        return scaled

#############################################################
# residual: Takes the data (Y values?) and subracts it from
# the composite. 
#############################################################
    """
    def residual(comp, data):
        return comp - data
    """
#############################################################
# The following function take the x,y,rms, and composite data    
# to plot the composite spectrum
#############################################################    
    def Composite(RF,RS,Names, LC, UC):

        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        #plt.axis([xmin, xmax, 0, max(RFmax)])
        #plt.minorticks_on()
        #plt.yticks(np.arange(0, 1.1, 0.2))
       
        for k in range(len_RF):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                plt.plot(RF[k][good], RF[k+1][good],label = Names[k] )
                if k == 0:
                    plt.fill_between(LC[k][good], LC[k+1][good], UC[k+1][good], color = 'c', alpha = 0.5)
                else:
                    plt.fill_between(LC[k][good], LC[k+1][good], UC[k+1][good], color = 'g', alpha = 0.5)
                #plt.plot(RF[k], RF[k+1], color = random.choice(['g', 'r', 'c', 'm', 'y', 'k']), label = Names[k] )
                #plt.fill_between(RF[k], RF[k+1] + RS[k+1], RF[k+1] - RS[k+1], facecolor = random.choice(['g', 'r', 'c', 'm', 'y', 'k']),alpha=0.5)                
                #plt.plot(RF[k], RF[k+1] + RS[1], label = "+ RMS")
                #plt.plot(RF[k], RF[k+1] - RS[1], label = "- RMS")

        """  
        for n in range(len_names):
            legend([Names[n]])
        """
        
        # Remove legend box frame        
        l = plt.legend(prop = {'family' : 'serif'})
        l.draw_frame(False)
        plt.draw()
          
        RFxticklabels = Rel_flux.get_xticklabels() 
        
        if max(stacked) == 0:
            plt.setp(RFxticklabels, visible=True)
        else:
            plt.setp(RFxticklabels, visible=False)
#############################################################
# The following function plots the variance 
#############################################################
    def Variance(VA):
        Variance = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Variance', fontdict = font)
        plt.yscale('log')
        #plt.axis([xmin, xmax, 1e-35, 1e-30])
        #plt.yticks(np.arange(1e-33, 1e-29))
        """
        for j in range(len_VA):
            if j % 2 == 0:
                VA[j], VA[j+1] = remove_extremes(VA[j],VA[j+1])
            
        VA = np.array(VA)
        """
        for k in range(len_VA):
            if k % 2 == 0:
                good = np.where((1/VA[k+1]) > 0)
                plt.plot(VA[k][good], VA[k+1][good], label = "Variance", ls = '-')
        VAxticklabels = Variance.get_xticklabels()
        
        # This isn't working the way I'd like it to         
        
        gca().yaxis.get_major_locator()._nbins=6
        
        if max(stacked) == 1:            
            plt.setp(VAxticklabels, visible=True)
        else:
            plt.setp(VAxticklabels, visible=False)
#############################################################
# The following function uses the x and rms to plot the RMS 
#############################################################
    def Residual(RS,LC, UC):
        Resid = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Residuals', fontdict = font)
        #plt.yticks(np.arange(0, 0.9, 0.2))
        for k in range(len_RS):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                comp = np.concatenate(RF[1][good])
		#There's something wrong with the dimensionality of x and y here
		#I added the extra [0]s because the arrays were 3 dimensional somehow, so now they're both 1-D
		#But they're still full of 'nan' so the plot gets messed up. But it runs through.
                plt.plot(RF[k][good], RF[k+1][good]-RF[1][good], label = "RMS of residuals", ls = '-')
                if k ==0:
                    plt.fill_between(LC[k][good], LC[k+1][good] - comp, UC[k+1][good] - comp, color = 'c', alpha = 0.5)
                else:
                    plt.fill_between(LC[k][good], LC[k+1][good] - comp, UC[k+1][good] - comp, color = 'g', alpha = 0.5)
        #plt.plot(RS[0], RS[1], label = "RMS of residuals", ls = '-')
        RSxticklabels = Resid.get_xticklabels()
        
        if max(stacked) == 2:
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
        #plt.yticks(np.arange(0, 0.9, 0.2))
        
        for k in range(len_SB):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                plt.plot(SB[k][good], SB[k+1][good], label = "Spectra per Bin", ls = '-')
        SBxticklabels = SpecBin.get_xticklabels()     
        
        gca().yaxis.get_major_locator()._nbins=6
        
        if max(stacked) == 3:
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
        #plt.yticks(np.arange(0, 0.9, 0.2))
        """
        for j in range(len_AG):
            if j % 2 == 0:
                AG[j], AG[j+1] = remove_extremes(AG[j],AG[j+1])
                
        AG = np.array(AG)
        """
        for k in range(len_AG):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                plt.plot(AG[k][good], AG[k+1][good], label = "Age", ls = '-')
        AGxticklabels = Age.get_xticklabels()        

        gca().yaxis.get_major_locator()._nbins=8
        
        if max(stacked) == 4:
            plt.setp(AGxticklabels, visible=True)
        else:
            plt.setp(AGxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the delta(?) with respect to wavelength(?).
#############################################################              
    def Delta(DE):
        Delta = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('$\Delta$m$_{15}$', fontdict = font)
        #plt.yticks(np.arange(0, 0.9, 0.2))
        """
        for j in range(len_DE):
            if j % 2 == 0:
                DE[j], DE[j+1] = remove_extremes(DE[j],DE[j+1])
        """      
        DE = np.array(DE)
        
        for k in range(len_DE):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                plt.plot(DE[k][good], DE[k+1][good], label = "Delta", ls = '-')
        DLxticklabels = Delta.get_xticklabels()
        
        gca().yaxis.get_major_locator()._nbins=8
        
        if max(stacked) == 5:            
            plt.setp(DLxticklabels, visible=True)
        else:
            plt.setp(DLxticklabels, visible=False)
#############################################################
# The following function will fill a smaller plot that shows  
# the redshift with respect to wavelength(?).
#############################################################              
    def Redshift(RD):
        Redshift = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Redshift', fontdict = font)
        #plt.yticks(np.arange(0, 0.9, 0.2))
        """
        for j in range(len_RD):
            if j % 2 == 0:
                RD[j], RD[j+1] = remove_extremes(RD[j],RD[j+1])
               
        RD = np.array(RD)
        """
        for k in range(len_RD):
            if k % 2 == 0:
                good = np.where(VA[k+1] != 0)
                plt.plot(RD[k][good], RD[k+1][good], label = "Redshift", ls = '-')
        Zxticklabels = Redshift.get_xticklabels()        

        gca().yaxis.get_major_locator()._nbins=6
        
        if max(stacked) >= 6:            
            plt.setp(Zxticklabels, visible=True)
        else:
            plt.setp(Zxticklabels, visible=False)

#############################################################
# The following function stacks all the data sets used to  
# build the composite spectrum. This will appear on a 
# separate figure.
#############################################################  
    def Stacked(RF):
#        plt.figure(num = 2, dpi = 100, figsize = [6, 4*len(RF[1])], facecolor = 'w')
        plt.figure(num = 2, dpi = 100, figsize = [10, 10], facecolor = 'w')

        
        buff = (max(RF[1])-min(RF[1]))/2.0
        
        plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.title('Stacked Figures', fontdict = font)  
        labels = ['+10 - +12 days', '', '+4 - +6 days', '', '-1 - +1 days', '', '-6 - -4 days', '', '-12 - 10 days', '']
        for m in range(len_RF):
            if m % 2 == 0:
                good = np.where(VA[m+1] != 0)
                plt.plot(RF[m][good], RF[m+1][good] + (m+1)*buff, label = labels[m])
                if m == 0:
                    plt.fill_between(LC[m][good], LC[m+1][good] + (m+1)*buff, UC[m+1][good] + (m+1)*buff, color = 'c', alpha = 0.5)
                else:
                    plt.fill_between(LC[m][good], LC[m+1][good] + (m+1)*buff, UC[m+1][good] + (m+1)*buff, color = 'g', alpha = 0.5)
#                plt.annotate(str(Names[m]), xy = (max(RF[m]), max(RF[m+1]+(m+1)*(1+buff))), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
                #plt.annotate(str(Names[m+1]), xy = (max(RF[0]), max(RF[1][m+1])+(m+1)*(1+buff)), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
        #plt.xlim(xmin, xmax) 
        plt.minorticks_on()
        plt.legend()
        plt.show(block = False)
#        plt.savefig(image_title[:-4] +'_stack.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
    
#############################################################
#remove_extremes will remove the peaks and dips from plotting
#############################################################
    def remove_extremes(dataX,dataY):
        # sort the data array by the flux
        sortedArray = sorted(zip(dataX,dataY), key = operator.itemgetter(1))

        # find the length of the array data. locate 95% and 05% indice
        length = len(dataX)
        newMax = np.int_(floor(length*.95))
        newMin = np.int_(ceil(length*.05))
        
        # cut the sorted array down to the remaining 5-95% range
        chopArray = []
        for i in range(newMin,newMax):
            chopArray.append(sortedArray[i])
            # sort by wavelength of the chopped array
        orderArray = sorted(chopArray, key = operator.itemgetter(0))
        # unzip the array. make it a single 2D array
        final = zip(*orderArray)
        return final
    
#############################################################
# The following function will take all the data sets to build  
# the composite spectrum and lay each one over the next. 
# This will appear on a separate figure.
############################################################# 
    # Temporary comment until we get data for multiple arrays
    """
    def Multi(RF,rms_data):
        plt.figure(num = 3, dpi = 100, figsize = [8, 8], facecolor = 'w')
        for m in range(len(RF[1])):
            plt.plot(xtrunc, RF[1][m], label = str(Names[m]))
        plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.minorticks_on()
        plt.legend(prop = {'family' : 'serif'})        
        plt.savefig('multi-spectrum plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
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
            #print median(RF[j])
            #print median(RF[j-1])
            RFtrunc = []
            for i in range(len(RF[j-1])):
#                if (RF[j-1][i] >= 4500) & (RF[j-1][i] <= 7000):
                if RF[j][i] != 0.0:
                    #print RF[j-1][i]
                    #print RF[j][i]
                    RFtrunc.append(RF[j][i])
            RF[j] = Scaling(RF[j], median(RFtrunc))
            LC[j] = Scaling(LC[j], median(RFtrunc))
            UC[j] = Scaling(UC[j], median(RFtrunc))
    
    
    #Not implemented until it can be fully tested
    
    """
    for j in range(len_RF):
        if j % 2 != 0:
            RF[j] = Scaling(RF[j])
    """
    """
    for j in range(len_VA):
        if j % 2 != 0:
            VA[j] = Scaling(VA[j])
    for j in range(len_RS):
        if j % 2 != 0:
            RS[j] = Scaling(RS[j])
    for j in range(len_SB):
        if j % 2 != 0:
            SB[j] = Scaling(SB[j])
    for j in range(len_AG):
        if j % 2 != 0:
            AG[j] = Scaling(AG[j])
    for j in range(len_DE):
        if j % 2 != 0:
            DE[j] = Scaling(DE[j])
    for j in range(len_SB):
        if j % 2 != 0:
            RD[j] = Scaling(RD[j])
    """
    
    
    
    """
    RFmax = []
    VAmax = []
    RSmax = []
    SBmax = []
    AGmax = []
    DEmax = []
    RDmax = []
    """
    
   
    """
    for j in range(len_VA):
        if j % 2 != 0:
            VAmax.append(max(VA[j])[0])
    for j in range(len_RS):
        if j % 2 != 0:
            RSmax.append(max(RS[j])[0])
    for j in range(len_SB):
        if j % 2 != 0:
            SBmax.append(max(SB[j])[0])
    for j in range(len_AG):
        if j % 2 != 0:
            AGmax.append(max(AG[j])[0])
    for j in range(len_DE):
        if j % 2 != 0:
            DEmax.append(max(DE[j])[0])
    for j in range(len_SB):
        if j % 2 != 0:
            RDmax.append(max(RD[j])[0])
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
#    plt.figure(num = 1, dpi = 100, figsize = [6, np.sum(h)], facecolor = 'w')
    plt.figure(num = 1, dpi = 100, figsize = [10,10], facecolor = 'w')
    gs = gridspec.GridSpec(len(Plots), 1, height_ratios = h, hspace = 0.001)
    
    mpl.rcParams['ytick.major.pad'] = 8
        
    p = 0
    Rel_flux = plt.subplot(gs[p])
    plt.title(title, fontdict = font)
    
#############################################################
# The following series of if statements runs a specific portion 
# of the ploting functions. The p = p+1 allows the code to  
# iterate through all the possible plot options      
# if 0 or 1 or 2 or 3 or 4 or 5 in Plots.
############################################################# 
    
    if 0 in Plots: # will always plot a composite spectrum if any 1-5 are selected   
        Composite(RF,RS,Names, LC, UC) 
        p = p+1
    if 1 in Plots:
        Variance(VA)
        p = p+1
    if 2 in Plots:
        Residual(RS, LC, UC) 
        p = p+1
    if 3 in Plots:  
        SpecBin(SB)  
        p = p+1               
    if 4 in Plots:
        Age(AG)
        p = p+1
    if 5 in Plots:
        Delta(DE)
        p = p+1
    if 6 in Plots:
        Redshift(RD)
        p = p+1
    plt.xlim(xmin, xmax)     
    plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
#    plt.savefig(image_title, dpi = 600, facecolor='w', edgecolor='w', pad_inches = 0.1)
    plt.show()
    # Initiates stacked plot
    if 7 in Plots:        
        Stacked(RF)
        p = p+1
    plt.show()
    #plt.axis([xmin, xmax, 0, max(RFmax)])    
    
    #fig = Figure()
    #fig.tight_layout()
    print "Plotting complete..."    