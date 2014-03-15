### Assuming the input in an array that you wish to make a histogram.
### Input: Query (argv[1]), bin size (argv[2])
### Updated: Mar 12, 2014


import numpy as np
import matplotlib.pyplot as plt
import sys
import sqlite3 as sq3

#connect to the database
con = sq3.connect('/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/MichaelSchubert/SNe.db')

#create a cursor object to execute commands
cur = con.cursor()


#The fields in the database are as follows
#Filename, SN name, Redshift, Phase, MinWave, MaxWave, Dm15, M_B, B_mMinusV_m

#Suppose you want to return the 10 SNe with the highest redshifts
#cur.execute('SELECT * FROM Supernovae ORDER BY ' + par + ' DESC LIMIT ' + sne_num)
cur.execute(sys.argv[1])

a = []

b = ['SN Name', 'Redshift', 'Julian Date', 'Julian Date Error', 'Light Curve', 'Delta 15', 'Delta 15 Error']
print b[0]

#The output is stored in the cursor object.
for row in cur:
    a.append(row[2])

bin = float(sys.argv[2])   ### Bin size
up = np.max(a)
low = np.min(a)
bin_num = int(np.ceil((up - low) / bin)) ### Number of bins


hist, bins = np.histogram(a, bins = bin_num)
width = 0.7 * (bins[1] - bins[0])   ### To leave some space among stripes
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width = width)
plt.ylabel(r'# of $SNe Ia$')
#plt.xlabel(par)
plt.show()
