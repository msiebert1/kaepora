# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:17:39 2014

@author: Brian
"""

"""
blank spaces need to be removed from in front of the first column of the data files in order for the code to properly read the data
"""

import numpy as np
import matplotlib.pyplot as plt

d1 = np.genfromtxt("sn2011by-hst+lick.txt", skiprows = 1, unpack = True)
d2 = np.genfromtxt("sn2011fe-visit3-hst.txt", skiprows = 1, unpack = True)

freq1 = d1[0]
flux1 = d1[1]
freq2 = d2[0]
flux2 = d2[1]

xmin = max(min(freq1),min(freq2))
xmax = min(max(freq1),max(freq2))

x_val = np.linspace(xmin, xmax, 1000)

int1 = np.interp(x_val, freq1, flux1)
int2 = np.interp(x_val, freq2, flux2)

int_ave = (int1 + int2)/2.0

scale1 = np.average(int1/int_ave)
scale2 = np.average(int2/int_ave)

plt.xlabel('freq (A)')
plt.ylabel('relative flux')
Plot = [plt.plot(freq1, flux1*1e13, label = "SN2011BY"), plt.plot(freq2, flux2*1e13, label = "SN2011FE"), plt.plot(x_val, int_ave*1e13, label = "Average")]
plt.legend()

plt.savefig('Fry spectra.png')
plt.show()

file = open("average spectra.txt", "w")
for w in range(len(x_val)):
    file.write(str(x_val[w]))
    file.write('\t')
    file.write(str(int_ave[w]))
    file.write('\n')
file.close()
