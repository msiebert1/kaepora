import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec


###########################################
#             HOW TO USE
## At the top of YOUR code write the following
# import Plotting
# 
## Place the following after you have calculated your data 
## fill in information
# 
# plot_data_1  = [wavelength_1, Spectra_1, plus, minus]     # Expects [X, Y, Pos_Y Err, Neg_Y Err]
# plot_data_2  = [wavelength_2, Spectra_2, plus, minus]     # Expects [X, Y, Pos_Y Err, Neg_Y Err]
# plot_data_3  = [wavelengths, scatter_p, scatter_n,]       # [ X , First Y_value, Second Y_value]
# Show_Data   = [plot_data_1, plot_data_2, plot_data_3]
# image_title  = "../personal/YOUR_REP/WHOA.png"            # Name the image (with location)
# plot_labels  = ["Clever Title","Wavelength ($\AA$)","AVG Spectrum","Scatter"]    #[ Figure Title, X title, Y-top Title, Y-bottom Title] 
# legend       = ["First","Second","Third","Fouth"]         # Names for the legend
#
# fig_type     = 1.1 or 1.2 or 2.3       # The types will change, this is only to give differnt plot options
#
## The following line will plot the data
#
# Plotting.main(Show_Data , image_title , plot_labels , legend)
#
#
###########################################

def main(Show_Data,Plots,image_title,title,legend_names):

    # Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
    #                   0              1          2            3    4      5
    Height =           [8,             2,         3,           2,   2,     2]

    #Plots = [0,1] # Plots to generate

    h = []

    for m in Plots:
        h.append(Height[m])
        
# re-name variables
    xaxis_1   = Show_Data[:][0] 
    yaxis_1   = Show_Data[:][1]
    err_p_1   = Show_Data[:][2] 
    err_n_1   = Show_Data[:][3]

#	xaxis_2   = Show_Data[:][1][0]
#	yaxis_2   = Show_Data[:][1][1]
	#err_p_2   = Show_Data[:][1][2]
	#err_n_2   = Show_Data[:][1][3]

	#xaxis_3   = Show_Data[:][2][0]
	#yaxis_3_1 = Show_Data[:][2][1]
	#yaxis_3_2 = Show_Data[:][2][2]

    legend_1  = legend_names[0]
    legend_2  = legend_names[1]
    legend_3  = legend_names[2]
    legend_4  = legend_names[3]


    #f = figure()
    #subplots_adjust(hspace=0.001)

# Changing font parameters. 
    params = {'legend.fontsize': 8,  
              'legend.linewidth': 2,
              'legend.font': 'serif',
              'xtick.labelsize': 12, 
              'ytick.labelsize': 12} # changes font size in the plot legend

    plt.rcParams.update(params) # reset the plot parameters

    font = {'family' : 'serif',   # Font Parameters
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 8,
            }
        
    plt.figure(num = 2, dpi = 100, figsize = [8, np.sum(h)], facecolor = 'w')
    plt.title(title)
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
        plt.plot(xaxis_1, yaxis_1, label = "generic data", ls = '-')
        plt.legend(prop = {'family' : 'serif'})
        plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
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
        plt.plot(xaxis_1, 2.0*xaxis_1, label = "title goes here", ls = '-')
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
        plt.plot(xaxis_1, 0*xaxis_1+3.0, label = "title goes here", ls = '-')
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
        plt.plot(xaxis_1, xaxis_1**3.0, label = "title goes here", ls = '-')
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
        plt.plot(xaxis_1, 1/xaxis_1, label = "title goes here", ls = '-')
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
        plt.plot(xaxis_1, xaxis_1**2.0-xaxis_1, label = "title goes here", ls = '-')
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

# Remove legend box frame 
    l = legend()
    l.draw_frame(False)
    plt.draw()

#Set the visual range. Automatic range is ugly. 
    xmin = int(float(xaxis_1[0]))
    xmax = int(float(xaxis_2[-1]))
    plt.xlim((xmin,xmax))

#Label the figure and show
    #plt.xlabel( xlabel )
    plt.savefig( image_title )
    plt.show()


"""
# Data should be read in and then plotted, make it so the user doens't have to choose the figure to plot
	if fig_type == 2.3: # Plotting details for two figures with a total of three sets of data
		ax1 = subplot(2,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

		# Plotting details for two figures with one set of data
		ax2 = subplot(2,1,2, sharex=ax1)
		plt.ylabel( ylabel_2 )
		ax2.plot(xaxis_3,yaxis_3_1,'b',label = legend_3 )
		ax2.plot(xaxis_3,yaxis_3_2,'g',label = legend_4 )
		
#Removes x-axis labels for ax1. May need later. 	
	# 	Commented code is to remove labels for both figures
	#	xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
	#	setp(xticklabels, visible=False)
		setp(ax1.get_xticklabels(), visible = False)
	
	elif fig_type == 1.2: # Plotting details for single figure with two sets of data
		ax1 = subplot(1,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'b',label = legend_1 )
		ax1.plot(xaxis_2,yaxis_2,'g',label = legend_2 )

		# Below fills the error around the plot
		plt.fill_between(xaxis_1,err_p_1,err_n_1,alpha=1.5, edgecolor='#000080', facecolor='#AFEEEE')
		plt.fill_between(xaxis_2,err_p_2,err_n_2,alpha=1.5, edgecolor='#006400', facecolor='#98FB98')

	else : # Plotting single plot with single set of data
		ax1 = subplot(1,1,1)
		plt.title( title )  # Placement of figure title in code is important
		plt.ylabel( ylabel_1 )
		ax1.plot(xaxis_1,yaxis_1,'k',label = legend_1 )
		plt.fill_between(xaxis_1,err_p_1,err_n_2,alpha=1.5, edgecolor='#000080', facecolor='#5F9EA0')

"""
"""

# Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
#                   0              1          2            3    4      5
Height =           [8,             2,         3,           2,   2,     2]

Plots = [0,1] # Plots to generate

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

xaxis_1 = np.linspace(2000, 8000, 12001)

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
    plt.plot(xaxis_1, xaxis_1**2.0, label = "generic data", ls = '-')
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
    plt.plot(xaxis_1, 2.0*xaxis_1, label = "title goes here", ls = '-')
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
    plt.plot(xaxis_1, 0*xaxis_1+3.0, label = "title goes here", ls = '-')
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
    plt.plot(xaxis_1, xaxis_1**3.0, label = "title goes here", ls = '-')
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
    plt.plot(xaxis_1, 1/xaxis_1, label = "title goes here", ls = '-')
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
    plt.plot(xaxis_1, xaxis_1**2.0-xaxis_1, label = "title goes here", ls = '-')
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



plt.show()


"""

