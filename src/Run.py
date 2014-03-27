import composite
import Plotting
from astropy.table import Table
#import targeted
"""
Here's the main function that anyone can use to run some analysis.

The first line runs the composite program.
This grabs the selected data from SNe.db.
Change the input to composite.main to whatever you want your database query to be.
This is the code that outputs a text file of data for later use.
It also shows/saves a basic plot of the composite data and associated errors.

The next section sets up the plots to be used, and then plots using Plotting.py
Right now it doesn't work...

More can be added as our programs evolve
"""

#This part works just fine
composite = composite.main("SELECT * FROM Supernovae WHERE snr > 8")

#Read whatever you sasved the table as
Data = Table.read('../plots/TestComposite', format='ascii')

#Checking to see how the table reads..right now it has a header that might be screwing things up.
print Data

#To be honest, I'm not sure entirely how this works.
#Can someone who worked on this piece of code work with it?
Relative_Flux = [Data[0], Data[1], composite.name]  # Want to plot a composite of multiple spectrum
Residuals     = [Data[0], Data[2], composite.name]
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
#Plotting.main(Show_Data , Plots, image_title , title)