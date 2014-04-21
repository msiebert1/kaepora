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
It mostly works at this point

More can be added as our programs evolve

There are a few parameters you should set beforehand.
plot_name is where the plot showing all composite spectra together will be saved.
title is the header on your plot
name_array is the list of legend labels
xmin and xmax define the boundary of the plot.
Plots actually defines which plots get shown
"""
queries = int(sys.argv[1])

d={}
for n in range(queries):
    if len(sys.argv) == queries+3:
        d["composite{0}".format(n+1)] = composite.main(sys.argv[n+2], sys.argv[queries+2])
    if len(sys.argv) == queries+4:
        d["composite{0}".format(n+1)] = composite.main(sys.argv[n+2], sys.argv[queries+2], sys.argv[queries+3])
    if len(sys.argv) == queries+5:
        d["composite{0}".format(n+1)] = composite.main(sys.argv[n+2], sys.argv[queries+2], sys.argv[queries+3], sys.argv[queries+4])
    else:
        d["composite{0}".format(n+1)] = composite.main(sys.argv[n+2])

#Read whatever you sasved the table as, iterates over however many composites you used.
#This is how you have to address things if you want to iterate over queries.
#The n+1 makes the first item composite1, not composite0.
for n in range(queries):
    d["data{0}".format(n+1)]         = Table.read(d["composite{0}".format(n+1)].savedname, format='ascii')
    d["wavelengths{0}".format(n+1)]  = np.array([d["data{0}".format(n+1)]["Wavelength"]])
    d["fluxes{0}".format(n+1)]       = np.array([d["data{0}".format(n+1)]["Flux"]])
    d["variances{0}".format(n+1)]    = np.array([d["data{0}".format(n+1)]["Variance"]])
    d["ages{0}".format(n+1)]         = np.array([d["data{0}".format(n+1)]["Age"]])
    d["dm15s{0}".format(n+1)]        = np.array([d["data{0}".format(n+1)]["Dm_15s"]])
    

#This chunk will create an array that's the right length for however many queries you used.
plot_array     = []
name_array     = []
residual_array = []
variance_array = []
age_array      = []
dm15_array     = []
for n in range(queries):
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
    name_array.append("composite{0}".format(n+1))
    name_array.append(" ")


#################
#Names of each line are, by default, composite1, composite2, etc.
#This is what shows up in the plot legend.
#If you want to use custom names for your composites,
#fill out and uncomment this next line (you need to have spaces between each name)
#name_array = ["name 1", " ", "name2", " ", etc]
#################

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
# Choose the plot range and plot type!
xmin         = 3000 
xmax         = 10100
Plots        = [0,1,2,4,5] 

#################
#Now the file name of the plot is labeled by time of creation, but you can personalize it if you want.
#Or rename it once it's been saved.

plot_name = str(queries) + '_composite_comparison, ' + (time.strftime("%H,%M,%S"))
image_title  = "../plots/" + plot_name + ".png"

# Want a custom saved image name? Uncomment this next line...but leave the one above alone.
#image_title  = "../plots/" + "Custome Title" + ".png"
print "Plot saved as: " + image_title

#Next line is the header on your plot
title        = "Insert Title Here"
################

Plotting.main(Show_Data, Plots, image_title, title, Names, xmin,xmax)
