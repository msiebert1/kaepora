import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *

#import matplotlib.gridspec as gridspec # Doesn't exist on EWS Linux

#import matplotlib.font_manage
#from matplotlib.font_manager import FontProperties
#from matplotlib.ticker import FuncFormatter

# Show_Data = [Light_curve, Residual, Age, Spectra_Bin, Age, Delta, Redshift] 

#def main(fig_type,height,plot_data_1,plot_data_2,plot_data_3,image_title,plot_labels,legend_names):

# Available Plots:  Relative Flux, Residuals, Spectra/Bin, Age, Delta, Redshift
	#                    0              1          2            3    4      5
	#Height =           [8,             2,         3,           2,   2,     2]

	#Plots = [0,1,2,3,4,5] # Plots to generate

	#h = []

	#for m in Plots:
	#    h.append(Height[m])

# re-name variables

#	for i in range(len(Plots)):
#		xaxis[] = Show_Data[:][0][0] 
def main(fig_type, Show_Data, Plots, image_title , plot_labels , legend):

    def remove_zero(s):
	print "Read in Show_Data"
	for i in range(len(s[:])):
		if any(s[:][i][1] == 0) :
"""
			delete.append(i)
			print "filling in delete.append"
			print len(delete)
			for i in range(len(delete)):
				print "Now remove zero data from array"
				del s[:][i][delete[len(delete)-1-i]]	#remove zero indices from the end of the data array
				del s[:][i][delete[len(delete)-1-i]]
				print "Done"
				print len(s)
"""

    remove_zero(Show_Data)

