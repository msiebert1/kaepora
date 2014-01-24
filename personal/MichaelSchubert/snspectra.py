from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('sn2011by-hst+lick.flm')
data2 = np.loadtxt('sn2011fe-visit3-hst.flm')

len1 = len(data1)
len2 = len(data2)



plt.plot(data1[:,0], data1[:,1], 'r')
plt.plot(data2[:,0], data2[:,1], 'b')
plt.show()
