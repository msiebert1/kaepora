import composite
import Plotting
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
import time
#import bootstrap
"""
Here's the main function that anyone can use to run some analysis.

It iterates through a loop based on how many queries were entered
Each one gets its data saved and the composite spectrum returned
Then all of the composites that were created get plotted

The next section sets up the plots to be used, and then plots using Plotting.py
Right now it doesn't work all the time...but it will

More can be added as our programs evolve

There are a few parameters you should set beforehand.
plot_name is where the plot showing both composite spectra together will be saved.
wmin and wmax define the boundary of the plot.
"""
def main(queries,plot_name,labels):
	num = int(queries[1])
	#Now the file name of the plot is labeled by time of creation, but you can personalize it if you want.
	#Or rename it once it's been saved.
	#plot_name = str(queries) + '_composite_comparison, ' + (time.strftime("%H,%M,%S"))
	wmin = 3000
	wmax = 10000

	d={}
	for n in range(num):
		if len(queries) == num+3:
			d["composite{0}".format(n+1)] = composite.main(queries[n+2], queries[num+2])
		if len(queries) == num+4:
			d["composite{0}".format(n+1)] = composite.main(queries[n+2], queries[num+2], queries[num+3])
		if len(queries) == num+5:
			d["composite{0}".format(n+1)] = composite.main(queries[n+2], queries[num+2], queries[num+3], queries[num+4])
		else:
			d["composite{0}".format(n+1)] = composite.main(queries[n+2])
			
	xmin =0
	xmax=100000
	for n in range(num):
		SN=d["composite{0}".format(n+1)]
		if (SN.minwave > xmin):
			xmin=SN.minwave
		if (SN.maxwave < xmax):
			xmax=SN.maxwave

	#Read whatever you sasved the table as, iterates over however many composites you used.
	#This is how you have to address things if you want to iterate over queries.
	#The n+1 makes the first item composite1, not composite0.
	for n in range(num):
		d["data{0}".format(n+1)]         = Table.read(d["composite{0}".format(n+1)].savedname, format='ascii')
		d["wavelengths{0}".format(n+1)]  = np.array([d["data{0}".format(n+1)]["Wavelength"]])
		d["fluxes{0}".format(n+1)]       = np.array([d["data{0}".format(n+1)]["Flux"]])
		d["variances{0}".format(n+1)]    = np.array([d["data{0}".format(n+1)]["Variance"]])
		d["ages{0}".format(n+1)]         = np.array([d["data{0}".format(n+1)]["Age"]])
		d["dm15s{0}".format(n+1)]        = np.array([d["data{0}".format(n+1)]["Dm_15s"]])
		

	# From now on, list the data you want to plot as [ Xdata, Ydata, Xdata_2, Ydata_2]
	#This chunk will create an array that's the right length for however many queries you used.
	plot_array     = []
	name_array     = []
	residual_array = []
	variance_array = []
	age_array      = []
	dm15_array     = []
	for n in range(num):
		plot_array.append(d["wavelengths{0}".format(n+1)])
		plot_array.append(d["fluxes{0}".format(n+1)])
		residual_array.append(d["wavelengths{0}".format(n+1)][0])
		residual_list = np.array([d["fluxes{0}".format(n+1)]-d["fluxes1"]])
		residual_array.append(residual_list[0][0])
		variance_array.append(d["wavelengths{0}".format(n+1)][0])
		variance_array.append(d["variances{0}".format(n+1)][0])
		age_array.append(d["wavelengths{0}".format(n+1)][0])
		age_array.append(d["ages{0}".format(n+1)][0])
		dm15_array.append(d["wavelengths{0}".format(n+1)][0])
		dm15_array.append(d["dm15s{0}".format(n+1)][0])
		name_array.append(labels[n])
		name_array.append(" ")

	#print variance_array # there were some problems with dimensionality, fixed now.

	##################
	#If you want to use custom names for your composites,
	#fill out and uncomment this next line
	#name_array = ["composite1name", " ", "composite2name", " ", etc]
	##################

	Relative_Flux   = plot_array #plots all given composites
	Variance        = variance_array # Check it out! Variances plot now.
	Residuals       = residual_array # Check it out! Residuals plot now.
	Spectra_Bin     = [] 
	Age             = age_array 
	Delta           = dm15_array
	Redshift        = []
	## If you want custom names, uncomment and use line 83, for consistency.
	##Otherwise it'll default to just labeling composites in order.
	Names           = name_array
	Show_Data       = [Relative_Flux,Variance,Residuals,Spectra_Bin,Age,Delta,Redshift]

	## Available Plots:  Relative Flux, Variances, Residuals, Spectra/Bin, Age, Delta, Redshift
	##                   0              1          2          3            4    5      6
	# the plots you want to create

	# All of these worked for me. Test with your own queries. (Sam, 4/16)
	# Choose the plot range and plot type!
	Plots        = [0] 

	image_title  = "../plots/" + plot_name + "_composites.png"
	title        = plot_name
	Plotting.main(Show_Data, Plots, image_title, title, Names, xmin,xmax)
