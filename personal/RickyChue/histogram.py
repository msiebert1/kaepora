### Assuming the input in an array that you wish to make a histogram.
### Input: Lower bound (argv[1]), upper bound (argv[2]), bin size (argv[3]), parameter that is to be binned (argv[4])
### Updated: Mar 12, 2014


import numpy as np
import matplotlib.pyplot as plt
import sys

low = float(sys.argv[1])   ### Lower bound of the bins
up = float(sys.argv[2])    ### Upper bound of the bins
bin = float(sys.argv[3])   ### Bin size
x_bin = sys.argv[4]        ### Parameter that is to be binned
bin_num = int(np.ceil((up - low) / bin)) ### Number of bins

a = np.loadtxt('redshift_hist.dat', unpack = True)  ### Input array for the histogram

hist, bins = np.histogram(a, bins = bin_num)
width = 0.7 * (bins[1] - bins[0])   ### To leave some space among stripes
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width = width)
plt.ylabel(r'# of $SNe Ia$')
plt.xlabel(x_bin)
plt.show()
