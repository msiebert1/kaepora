# Load the data files and produce an average spectrum
#
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import interpolate

def halfSearch(arr,find):
    median = 0
    bottom = 0
    ceil = len(arr) - 1
    while (bottom <= ceil):
        median = (bottom + ceil) / 2
        temp  = arr[median]
        if (temp == find):
            return median 
        else:
            if (temp < find):
                bottom = median + 1
            else:
                ceil = median - 1
    return -1 #not found


#Loading data from data files
def loadData(file):
    datadir = '../../../data/'
    f = open(datadir+file)
    lines = f.readlines()
    f.close()

    x = []
    y = []

    for line in lines:
        p = line.split()
        x.append(float(p[0]))
        y.append(float(p[1]))

    x = np.array(x)
    y = np.array(y)        

    return np.array([x,y])


file1 = 'sn2011by-hst+lick.flm'
data = loadData(file1)
x1 = data[0]
y1 = data[1]

file2 = 'sn2011fe-visit3-hst.flm'
data = loadData(file2)
x2 = data[0]
y2 = data[1]

#Truncate the extra wavelength of SN2011by
low = halfSearch(x1,x2[0])
high = halfSearch(x1,x2[len(x2)-1])
x1 = x1[low:high+1] #Caution: wired +1 here,different from IDL
y1 = y1[low:high+1]

#De-redshift and resample by B-spline interpolation 
z1 = 0.003402
z2 = 0.001208
x1 /= 1+z1
x2 /= 1+z2

xmin = max(x1[0],x2[0])
xmax = min(x1[-1],x2[-1])
xs = scipy.linspace(xmin,xmax,len(x1)*2)
tck1 = interpolate.splrep(x1, y1)
tck2 = interpolate.splrep(x2, y2)
x1 = xs
x2 = xs
y1 = interpolate.splev(xs,tck1)
y2 = interpolate.splev(xs,tck2)

#Normalize y
nfac = np.median(y1) #normalization factor
y1 = y1/nfac

nfac = np.median(y2)
y2 = y2/nfac

#Average two spectrum
ax = xs
ay = np.mean(np.array([y1,y2]), axis=0)

#plot
pltdir = '../plots/'

plt.yscale('log')
p1,=plt.plot(x1,y1)
p2,=plt.plot(x2,y2)
p3,=plt.plot(ax,ay)

plt.xlabel('Wavelength [A]')
plt.ylabel('Scaled Flux')
plt.legend([p1,p2,p3],['SN2011BY','SN2011FE','Average'],
           4,)
plt.savefig(pltdir+'spectrum.eps')
plt.show()

