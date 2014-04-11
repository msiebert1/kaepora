import composite
import Plotting
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
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

avgphase = (d["composite1"].phase + d["composite2"].phase)/2
avgred   = (d["composite1"].redshift + d["composite2"].redshift)/2

plot_name = '2_composite_comparison, ' + 'avgphase, ' + str(avgphase) + ', avgred, ' + str(avgred)
wmin = 3000
wmax = 10000

#This makes, shows, and saves a quick comparison plot...we can probably get rid of this when plotting.main works.
lowindex  = np.where(d["composite1"].wavelength == composite.find_nearest(d["composite1"].wavelength, wmin))
highindex = np.where(d["composite1"].wavelength == composite.find_nearest(d["composite1"].wavelength, wmax))
#factor = np.mean(d["composite1"].flux[lowindex[0]:highindex[0]]/d["composite2"].flux[lowindex[0]:highindex[0]])
for n in range(queries):
    plt.plot(d["composite{0}".format(n+1)].wavelength[lowindex[0]:highindex[0]], d["composite{0}".format(n+1)].flux[lowindex[0]:highindex[0]])
plt.savefig('../plots/' + plot_name + '.png')
plt.show()

#Read whatever you sasved the table as
Data1 = Table.read(d["composite1"].savedname, format='ascii')
Data2 = Table.read(d["composite2"].savedname, format='ascii')

#Checking to see how the table reads..right now it has a header that might be screwing things up.
#print Data
wavelengths  = np.array([Data1["Wavelength"]])
fluxes1      = np.array([Data1["Flux"]])
variances1   = np.array([Data1["Variance"]])
fluxes2      = np.array([Data2["Flux"]])
variances2   = np.array([Data2["Variance"]])

#print wavelengths

# From now on, list the data you want to plot as [ Xdata, Ydata, Xdata_2, Ydata_2]
Relative_Flux   = [wavelengths, fluxes1, wavelengths, fluxes2]  # Want to plot a composite of multiple spectrum
Residuals       = [wavelengths, variances1]
Spectra_Bin     = [] 
Age             = [] 
Delta           = [] 
Redshift        = []
Names           = [d["composite1"].name,d["composite1"].name] # the names corresponding to each composite go here
Show_Data       = [Relative_Flux, Residuals, Spectra_Bin, Age , Delta , Redshift] # Removed residual section 
## Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
##                   0              1          2            3    4      5         
# the plots you want to create
Plots        = [0,1] 
image_title  = "../plots/Composite_Spectrum_plotted.png"				 
title        = "Composite Spectrum"	
# The following line will plot the data
Plotting.main(Show_Data , Plots , image_title , title, Names)
