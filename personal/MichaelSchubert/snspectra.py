from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

def find_boundaries(arr1, arr2):
    """
    Takes 2 1d arrays and returns the indicies between which
    they overlap.
    """
    start = arr1[0]
    if arr1[0] > arr2[0]:
        start = arr2[0]

    end = arr1[len(arr1) - 1]
    if arr1[len(arr1) - 1] > arr2[len(arr2) - 1]:
        end = arr2[len(arr2) - 1]

    return start, end

data1 = np.loadtxt('sn2011by-hst+lick.flm')
data2 = np.loadtxt('sn2011fe-visit3-hst.flm')

start, end = find_boundaries(data1[:,0], data2[:,0])
print start, end

plt.plot(data1[:,0], data1[:,1], 'r')
plt.plot(data2[:,0], data2[:,1], 'b')
plt.show()
