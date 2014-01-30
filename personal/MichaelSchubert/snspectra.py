from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

def find_boundaries(arr1, arr2):
    """
    Takes 2 1d arrays and returns the indicies between which
    they overlap.
    """
    start = arr1[0]
    if arr1[0] < arr2[0]:
        start = arr2[0]

    end = arr1[len(arr1) - 1]
    if arr1[len(arr1) - 1] > arr2[len(arr2) - 1]:
        end = arr2[len(arr2) - 1]

    return start, end


def trim_arrays(arr1, arr2):
    """
    Returns section of each array that overlaps
    """
    start, end = find_boundaries(arr1[:,0], arr2[:,0])

    ind1 = np.where(arr1 == start)
    ind2 = np.where(arr2 == start)
    end1 = np.where(arr1 == end)
    end2 = np.where(arr2 == end)

    g1 = data1[ind1[0]:end1[0]]
    g2 = data2[ind2[0]:end2[0]]

    return g1, g2


data1 = np.loadtxt('sn2011by-hst+lick.flm')
data2 = np.loadtxt('sn2011fe-visit3-hst.flm')

g1, g2 = trim_arrays(data1, data2)

#scales data by dividing by median
gmed1 = np.median(g1[:,1])
gmed2 = np.median(g2[:,1])
sg1 = np.divide(g1[:, 1], gmed1)
sg2 = np.divide(g2[:, 1], gmed2)

#average flux
avx = g1[:, 0]
avy = np.divide(np.add(sg1, sg2), 2)

#plot trimmed wavelength vs. trimmed/scaled flux
p1 ,= plt.plot(avx, avy, 'g')
p2 ,= plt.plot(g1[:,0], sg1, 'r')
p3 ,= plt.plot(g2[:,0], sg2, 'b')

plt.legend([p1, p2, p3], ['Average', 'sn2011by', 'sn2011fe'], 4)
plt.xlabel('Wavelength [A]')
plt.ylabel('Scaled Flux')
plt.yscale('log')
plt.title('SNe spectra (2011by, 2011fe)')
plt.savefig('spectra.png')
plt.show()