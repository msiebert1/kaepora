import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec
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
    Height =           [6,             2,         2,           2,        2,     2,    2,        0]
    
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

    # Even values are x. Odd are y  (Slightly confusing for the time being)
    RF = []
    RS = []
    SB = []
    AG = []
    DE = []
    RD = []
    VA = []

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
    params = {'legend.fontsize': 10, 
              'legend.linewidth': 2,
              'legend.font': 'serif',
              'mathtext.default': 'regular', 
              'xtick.labelsize': 10, 
              'ytick.labelsize': 10} # changes font size in the plot legend

    plt.rcParams.update(params)                             # reset the plot parameters

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
    def Composite(RF,RS,Names):
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.axis([xmin, xmax, 0, max(RFmax)])
        plt.minorticks_on()
        #plt.yticks(np.arange(0, 1.1, 0.2))
       
        for k in range(len_RF):
            if k % 2 == 0:
                plt.plot(RF[k], RF[k+1], label = Names[k] )
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
        #plt.yticks(np.arange(0, 0.9, 0.2))

        for k in range(len_VA):
            if k % 2 == 0:
                plt.plot(VA[k], VA[k+1], label = "Variance", ls = '-')
        VAxticklabels = Variance.get_xticklabels()
        if max(stacked) == 1:            
            plt.setp(VAxticklabels, visible=True)
        else:
            plt.setp(VAxticklabels, visible=False)
#############################################################
# The following function uses the x and rms to plot the RMS 
#############################################################
    def Residual(RS):
        Resid = plt.subplot(gs[p], sharex = Rel_flux)
        plt.ylabel('Residuals', fontdict = font)
        #plt.yticks(np.arange(0, 0.9, 0.2))
        
        for j in range(len_RS):
            if j % 2 == 0:		
		#There's something wrong with the dimensionality of x and y here
		#I added the extra [0]s because the arrays were 3 dimensional somehow, so now they're both 1-D
		#But they're still full of 'nan' so the plot gets messed up. But it runs through.
                plt.plot(RF[j], RF[j+1]-RF[1], label = "RMS of residuals", ls = '-')
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
                plt.plot(SB[k], SB[k+1], label = "Spectra per Bin", ls = '-')
        SBxticklabels = SpecBin.get_xticklabels()        
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
        
        for k in range(len_AG):
            if k % 2 == 0:
                plt.plot(AG[k], AG[k+1], label = "Age", ls = '-')
        AGxticklabels = Age.get_xticklabels()        
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
        plt.ylabel('$\Delta$', fontdict = font)
        #plt.yticks(np.arange(0, 0.9, 0.2))
        
        for k in range(len_DE):
            if k % 2 == 0:
                plt.plot(DE[k], DE[k+1], label = "Delta", ls = '-')
        DLxticklabels = Delta.get_xticklabels()
        if max(stacked) == 5:            
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
        #plt.yticks(np.arange(0, 0.9, 0.2))
        
        for k in range(len_RE):
            if k % 2 == 0:
                plt.plot(RE[k], RE[k+1], label = "Redshift", ls = '-')
        Zxticklabels = Redshift.get_xticklabels()        
        if max(stacked) == 6:            
            plt.setp(Zxticklabels, visible=True)
        else:
            plt.setp(Zxticklabels, visible=False)
#############################################################
#remove_extremes will remove the peaks and dips from plotting
#############################################################
    def remove_extremes(dataX,dataY):
        # sort the data array by the flux
        sortedArray = sorted(zip(dataX,dataY), key = operator.itemgetter(1))

        # find the length of the array data. locate 95% and 05% indice
        length = len(data[0])
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
        plt.figure(num = 2, dpi = 100, figsize = [8, 8], facecolor = 'w')
        for m in range(len(RF[1])):
            plt.plot(xtrunc, RF[1][m], label = str(Names[m]))
        plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
        plt.ylabel('Relative, f$_{\lambda}$', fontdict = font)
        plt.minorticks_on()
        plt.legend(prop = {'family' : 'serif'})        
        plt.savefig('multi-spectrum plot.png', dpi = 100, facecolor='w', edgecolor='w', pad_inches = 0.1)
        """
#############################################################
# The following function stacks all the data sets used to  
# build the composite spectrum. This will appear on a 
# separate figure.
#############################################################  
    """  
    # Commented out till testing is complete.
    def Stacked(RF): 
       len_RF = 2
       plt.figure(num = 2, dpi = 100, figsize = [8, 4*len(RF[1])], facecolor = 'w')
       plt.plot(RF[0], RF[1])
       #plt.annotate(str(Names[0]), xy = (max(RF[0]), max(RF[1][0])), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       buffer = (max(RF[1])-min(RF[1]))/2.0
       for m in range(len_RF-1):
           if m % 2 == 0:
               plt.plot(RF[m], RF[m+1]+(m+1)*(1+buffer))
               #plt.annotate(str(Names[m+1]), xy = (max(RF[0]), max(RF[1][m+1])+(m+1)*(1+buffer)), xytext = (-10, 0), textcoords = 'offset points', fontsize = 8, family  = 'serif', weight = 'bold', ha = 'right')
       
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
            #print median(RF[j])
            #print median(RF[j-1])
            RFtrunc = []
            for i in range(len(RF[j-1])):
                if (RF[j-1][i] >= 4500) & (RF[j-1][i] <= 7000):
                    #print RF[j-1][i]
                    #print RF[j][i]
                    RFtrunc.append(RF[j][i])
            RF[j] = Scaling(RF[j], median(RFtrunc))
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
    
    RFmax = []
    
    """
    VAmax = []
    RSmax = []
    SBmax = []
    AGmax = []
    DEmax = []
    RDmax = []
    """
    
    for j in range(len_RF):
        if j % 2 != 0:
            RFmax.append(max(RF[j])[0])
    # The following code is to remove extreme data points
    # the random peaks and dips. This will not be 
    # implemented until code can be run :(
    """
    for k in range(len_RF):
	if j % 2 != 0:
	    combinArray_RF = []
	    combinArray_RF = remove_extremes(RF[k],RF[k+1])
	    RF[k].append(combinArray_RF[0])
	    RF[k+1].append(combinArray_RF[1])

    for k in range(len_RS):
	if j % 2 != 0:
	    combinArray_RS = []
 	    combinArray_RS = remove_extremes(RS[k],RS[k+1])
	    RS[k].append(combinArray_RS[0])
	    RS[k+1].append(combinArray_RS[1])

    for k in range(len_SB):
	if j % 2 != 0:
	    combinArray_SB = []
	    combinArray_SB = remove_extremes(SB[k],SB[k+1])
	    SB[k].append(combinArray_SB[0])
	    SB[k+1].append(combinArray_SB[1])

    for k in range(len_AG):
	if j % 2 != 0:
	    combinArray_AG = []
 	    combinArray_AG = remove_extremes(AG[k],AG[k+1])
	    AG[k].append(combinArray_AG[0])
	    AG[k+1].append(combinArray_AG[1])
    
    for k in range(len_DE):
	if j % 2 != 0:
	    combinArray_DE = []
	    combinArray_DE = remove_extremes(DE[k],DE[k+1])
	    DE[k].append(combinArray_DE[0])
	    DE[k+1].append(combinArray_DE[1])

    for k in range(len_RD):
	if j % 2 != 0:
	    combinArray_RD = []
 	    combinArray_RD = remove_extremes(RD[k],RD[k+1])
	    RD[k].append(combinArray_RD[0])
	    RD[k+1].append(combinArray_RD[1])
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
    plt.figure(num = 1, dpi = 100, figsize = [6, np.sum(h)], facecolor = 'w')
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
        Variance(VA)
        p = p+1
    if 2 in Plots:
        Residual(RS) 
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
    """ # Initiate stacked plot
    if 7 in Plots:
        Stacked(RF)
        p = p+1
    """   
    # Regardless of what is plotted, we label the Xaxis and save the plot image
    plt.xlim(xmin, xmax)       
    plt.xlabel('Rest Wavelength [$\AA$]', fontdict = font)
    #plt.axis([xmin, xmax, 0, max(RFmax)])    
    plt.savefig(image_title, dpi = 600, facecolor='w', edgecolor='w', pad_inches = 0.1) # CHANGE dpi = 600
    print "Plotting complete..."    
    
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
# Show the plots!
#############################################################                  

    plt.show()


    """
    # Commented out for the time being because we're only given a single array
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

    #Will implement once testing is complete
    """
    def remove_extremes(data):
        new_data = []
        bad      = []
        new_data = data
        length = len(data) #do i need it?
        sort(data[1])
        new_max = data[1][-1]*.95
        new_min = data[1][0]*.05
        for i in range(length):
            if (new_data[1][i] > new_max) or (new_data[1][i] < new_min):
                bad.append(i)
        for j in range(len(bad)):
            del new_data[0][bad[len(bad)-1-j]]
            del new_data[1][bad[len(bad)-1-j]]
        return new_data
        
        
        def chop_data(d,p):
	delete = []	#temporary variable to track which indices to remove

	for i in range(len(d)):
		if d[i][0][0] > wave_min or d[i][-1][0] < wave_max:
			delete.append(i)
			#print "\nwaves #",i,":",data_pos[i][:,0]
			#print "from path:",path_pos[i],"does not make the cut"
	for i in range(len(delete)):
		#print "remove",d[delete[len(delete)-1-i]]
		del d[delete[len(delete)-1-i]]	#have to remove from end of array to account for moving indexes
		del p[delete[len(delete)-1-i]]
		#print d[delete[len(delete)-1-i]]
      """

