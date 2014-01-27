# Load the data files and produce an average spectrum
#
import matplotlib.pyplot as plt
import numpy as np

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

#Load data
datadir = '../../../data/'

f1 = open(datadir+'sn2011by-hst+lick.flm')
lines = f1.readlines()
f1.close()

x1 = []
y1 = []
for line in lines:
    p = line.split()
    x1.append(float(p[0]))
    y1.append(float(p[1]))

x1 = np.array(x1)
y1 = np.array(y1)

f2 = open(datadir+'sn2011fe-visit3-hst.flm')
lines = f2.readlines()
f2.close()

x2 = []
y2 = []
for line in lines:
    p = line.split()
    x2.append(float(p[0]))
    y2.append(float(p[1]))

x2 = np.array(x2)
y2 = np.array(y2)

#Match the frequencies
low = halfSearch(x1,x2[0])
high = halfSearch(x1,x2[len(x2)-1])
x1 = x1[low:high+1] #Caution: wired +1 here,different from IDL
y1 = y1[low:high+1]

#Normalize y
mediany1 = np.median(y1)
y1 = y1/mediany1 

mediany2 = np.median(y2)
y2 = y2/mediany2

#Average two spectrum
ax = x1
ay = (y1+y2)/2

#plot
pltdir = '../plots/'

plt.yscale('log')
p1,=plt.plot(x1,y1)
p2,=plt.plot(x2,y2)
p3,=plt.plot(ax,ay)

plt.xlabel('Wavelength [A]')
plt.ylabel('Scaled Flux')
plt.legend([p1,p2,p3],['SN2011BY','SN2011FE','Avarage'],
           4,)
plt.savefig(pltdir+'spectrum.eps')
plt.show()

