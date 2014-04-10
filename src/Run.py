import composite
import Plotting
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

There are a few parameters you should set beforehand.
plot_name is where the plot showing both composite spectra together will be saved.
wmin and wmax define the boundary of the plot.
"""



#This part works just fine
#composite_full = composite.main("SELECT * FROM Supernovae WHERE Redshift > 0 AND Phase > -100")
"""
Here we set the queries that get used to find the spectra for compositing.
We only want to select spectra that have data for both redshift and phase, so both of them need to be in the query.
But you can change the values to whatever you want, and add more parameters.
"""
composite1 = composite.main("SELECT * FROM Supernovae WHERE Redshift > -100 AND Phase BETWEEN -7 AND -3 AND Velocity BETWEEN -11 AND -10")
composite2 = composite.main("SELECT * FROM Supernovae WHERE Redshift > -100 AND Phase BETWEEN -3 AND 3 AND Velocity BETWEEN -11 AND -10")
#composite3 = composite.main("SELECT * FROM Supernovae WHERE Redshift > -100 AND Phase BETWEEN 3 AND 7 AND Velocity BETWEEN -11 AND -10")
#composite4 = composite.main("SELECT * FROM Supernovae WHERE Redshift > -100 AND Phase BETWEEN 7 AND 14 AND Velocity BETWEEN -11 AND -10")
#composite5 = composite.main("SELECT * FROM Supernovae WHERE Redshift > -100 AND Phase BETWEEN 14 AND 24 AND Velocity BETWEEN -11 AND -10")

avgphase = (composite1.phase + composite2.phase)/2
avgred = (composite1.redshift + composite2.redshift)/2

plot_name = '2_composite_comparison, ' + 'avgphase, ' + str(avgphase) + ', avgred, ' + str(avgred)
wmin = 3000
wmax = 10000

#This makes, shows, and saves a quick comparison plot...we can probably get rid of this when plotting.main works.
lowindex = np.where(composite1.wavelength == composite.find_nearest(composite1.wavelength, wmin))
highindex = np.where(composite1.wavelength == composite.find_nearest(composite1.wavelength, wmax))
#factor = np.mean(composite1.flux[lowindex[0]:highindex[0]]/composite2.flux[lowindex[0]:highindex[0]])
plt.plot(composite1.wavelength[lowindex[0]:highindex[0]], composite1.flux[lowindex[0]:highindex[0]])
plt.plot(composite2.wavelength[lowindex[0]:highindex[0]], composite2.flux[lowindex[0]:highindex[0]])
#plt.plot(composite3.wavelength[lowindex[0]:highindex[0]], composite3.flux[lowindex[0]:highindex[0]])
#plt.plot(composite4.wavelength[lowindex[0]:highindex[0]], composite4.flux[lowindex[0]:highindex[0]])
#plt.plot(composite5.wavelength[lowindex[0]:highindex[0]], composite5.flux[lowindex[0]:highindex[0]])
plt.savefig('../plots/' + plot_name + '.png')
plt.show()

#Read whatever you sasved the table as
Data1 = Table.read(composite1.savedname, format='ascii')
Data2 = Table.read(composite2.savedname, format='ascii')

#Checking to see how the table reads..right now it has a header that might be screwing things up.
#print Data
wavelengths  = np.array([Data1["Wavelength"]])
fluxes1      = np.array([Data1["Flux"]])
variances1   = np.array([Data1["Variance"]])
fluxes2      = np.array([Data2["Flux"]])
variances2   = np.array([Data2["Variance"]])

#print wavelengths

#To be honest, I'm not sure entirely how this works.
#Can someone who worked on this piece of code work with it?
Relative_Flux = [wavelengths, fluxes1, composite1.name]  # Want to plot a composite of multiple spectrum
Residuals     = [wavelengths, variances1, composite1.name]
Spectra_Bin   = []
Age           = []
Delta         = []
Redshift      = [] 
Show_Data     = [Relative_Flux,Residuals]
image_title  = "../plots/Composite_Spectrum_plotted.png"            # Name the image (with location)
title        = "Composite Spectrum" 

# Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift, Multiple Spectrum, Stacked Spectrum
#                   0              1          2            3    4      5         6,                 7
Plots = [0] # the plots you want to create

# The following line will plot the data
#It's commented out until it works...
Plotting.main(Show_Data , Plots, image_title , title)