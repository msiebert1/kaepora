from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os

#consider changing to binary search for O(logn)
def find_lo_hi(lo, hi, arr):
    d = arr[:, 0]

    if d[0] > lo:
        lo = d[0]
    if d[len(d) - 1] < hi:
        hi = d[len(d) - 1]

    print lo, hi
    return lo, hi

def trim_data(lo, hi, arr):
    """
    Returns section of each array between given high and low values
    """
    start = np.where(arr == lo)
    end = np.where(arr == hi)
    trimmed = arr[start[0]:end[0]]
    return trimmed

root = '../../data/'

#initialize lo and hi to min/max possible
lo = 0.0
hi = float('inf')
bad = 0
count = 0
for path, subdirs, files in os.walk(root):
    for name in files:
        f = os.path.join(path, name)
        count += 1
        if f.endswith('.flm'):
            #ignore spectra that produce loaderrors
            try:
                data = np.loadtxt(f)
            except:
                bad += 1
                continue
        else:
            continue
        lo, hi = find_lo_hi(lo, hi, data)
print bad
print count
spectra = {}

for f in files:
    if f.endswith('.flm'):
        data = np.loadtxt(f)
    else:
        continue
    trimmed = trim_data(lo, hi, data)
    #scale data by dividing by median
    med = np.median(data[:, 1])
    scaled = np.divide(trimmed, med)
    spectra[f] = trimmed
    print spectra

#average flux
avx = g1[:, 0]
#avy = np.divide(np.add(sg1, sg2), 2)
avy = np.divide((sg1+sg2), 2)

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