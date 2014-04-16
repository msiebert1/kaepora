
"""
How to use:

To plot the historgram, first one needs to find the collection of samples needed to be plotted. They are selected via standard SQL query. One can make multiple queries at once. Then one needs to specify over which parameter the histogram is plotted and then input the bin size and range of the plot. The format of the input is:

   ' number of '

  The selection of parameters are: 'Filename' 'SN' 'Redshift' 'Phase' 'MinWave' 'MaxWave' 'RedMin' 'RedMax'
  'Dm15' 'M_B' and 'B_mMinusV_m'
  
  An example is:
  
            'SELECT * FROM Supernovae ORDER BY Phase DESC LIMIT 100' 0.001 Redshift
  
  
  Which will give you 100 spectra with highest phases plotted in the histogram based on their redshift, with bin size of 0.001

Be careful with the combination of the last two inputs, bin size and parameters. For example, '0.001' is a good bin size for redshift, but '0.1' is a good bin size for Delta15(Dm15), and '100' is a good bin size of MaxWave.



"""

'''
import numpy as np
import matplotlib.pyplot as plt
import sys
import sqlite3 as sq3
from matplotlib.patches import Ellipse, Polygon



#connect to the database
con = sq3.connect('SNe.db')

#create a cursor object to execute commands
cur = con.cursor()
cur2 = con.cursor()


#The fields in the database are as follows
#Filename, SN name, Redshift, Phase, MinWave, MaxWave, Dm15, M_B, B_mMinusV_m

#Suppose you want to return the 10 SNe with the highest redshifts
#cur.execute('SELECT * FROM Supernovae ORDER BY ' + par + ' DESC LIMIT ' + sne_num)
cur.execute(sys.argv[1])
top = cur2.execute('Select * from Supernovae')
d=list(top)

b2 = np.asarray([x[0] for x in top.description], dtype='object')
b = []
for i in range(b2.size):
    b.append(b2[i])

a = []


# the parameter wanted to be plotted.
param = sys.argv[3]

if any(param for param in b):
    index = b.index(param)



#The output is stored in the cursor object, based on the input parameter.
for row in cur:
    if row[index] == None:
        continue
    else:
        a.append(row[index])
print a, 'there are ' + str(len(a)) + ' spectra with measured data'



bin = float(sys.argv[2])   ### Bin size
up = np.max(a)
low = np.min(a)
bin_num = int(np.ceil((up - low) / bin)) ### Number of bins


params = {'legend.fontsize': 8,
    'legend.linewidth': 2,
    'legend.font': 'serif',
    'mathtext.default': 'regular',
    'xtick.labelsize': 8,
    'ytick.labelsize': 8} # changes font size in the plot legend

plt.rcParams.update(params)                             # reset the plot parameters

font = {'family' : 'serif',
    'color'  : 'black',
    'weight' : 'bold',
    'size'   : 10,
}

hist, bins = np.histogram(a, bins = bin_num)
width = 0.7 * (bins[1] - bins[0])   ### To leave some space among stripes
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width = width, color='cyan', edgecolor='blue', hatch="//")
plt.ylabel(r'# of $SNe Ia$', fontdict = font)
plt.xlabel(param, fontdict = font)
plt.show()
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import sqlite3 as sq3

# Number of queries you want
que_num = int(sys.argv[1])

#connect to the database
con = sq3.connect('SNe.db')

#create a cursor object to execute commands
cur_header = con.cursor()

#The fields in the database are as follows
#Filename, SN name, Redshift, Phase, MinWave, MaxWave, Dm15, M_B, B_mMinusV_m

# Name of the headers of the table
top = cur_header.execute('Select * from Supernovae')
d=list(top)

b2 = np.asarray([x[0] for x in top.description], dtype='object')
b = []

for i in range(b2.size):
    b.append(b2[i])

# the parameter wanted to be plotted.
param = sys.argv[2]

if any(param for param in b):
    index = b.index(param)

cur = [0] * que_num
a = [0] * que_num

# A for loop for all the queries.
for m in range(que_num):
    cur[m] = con.cursor()
    cur[m].execute(sys.argv[m + 6])
    a[m] = []
    
    #The output is stored in the cursor object, based on the input parameter.
    for col in cur[m]:
        if col[index] == None:
            continue
        else:
            a[m].append(col[index])
    
#print 'There are ' + str(len(a[m])) + ' spectra with measured data'

bin = float(sys.argv[3])    ### Bin size
up = float(sys.argv[4])     ### Upper bound of the histogram
low = float(sys.argv[5])    ### Lower bound of the histogram
bin_num = int(np.ceil((up - low) / bin)) ### Number of bins

hist = [0] * que_num
bins = [0] * que_num
leg_name = [0] * que_num

colors = iter(cm.rainbow(np.linspace(0, 1, len(a))))

hatch = ["////","\\\\","||||","___","++"]

hatindex = hatch.index

bold = "\033[1m"

for m in range(que_num):
    hist[m], bins[m] = np.histogram(a[m], bins = bin_num, range=(low, up))
    width = 0.7 * bin   ### To leave some space among stripes
    center = (bins[m][:-1] + bins[m][1:]) / 2
    plt.bar(center, hist[m], align='center', width = width, alpha = 0.3, color = next(colors), hatch= hatch[m])
    plt.xlim((low, up))
    ### Insert the name of the legends
    if ('ORDER' in sys.argv[m + 6].split()):
        leg_name[m] = sys.argv[m + 6].split()[sys.argv[m + 6].split().index('ORDER') + 2] + "   " + str(len(a[m])) + ' spectra available'
    else:
        leg_name[m] = str(len(a[m])) + ' spectra available'
    #leg_name[m] = sys.argv[m + 6].split()[6] + ' ' + sys.argv[m + 6].split()[7] + ' ' + sys.argv[m + 6].split()[9]+ "   " + str(len(a[m])) + ' spectra available'

params = {'legend.fontsize': 8,
    'legend.linewidth': 2,
    'legend.font': 'serif',
    'mathtext.default': 'regular',
    'xtick.labelsize': 8,
    'ytick.labelsize': 8} # changes font size in the plot legend

plt.rcParams.update(params)                             # reset the plot parameters

font = {'family' : 'serif',
    'color'  : 'black',
    'weight' : 'bold',
    'size'   : 10,
}

plt.ylabel(r'# of $SNe Ia$', fontdict = font)
plt.xlabel(param, fontdict = font)
plt.legend(leg_name)
plt.show()



