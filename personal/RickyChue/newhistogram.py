"""
    How to use:
    
    The input command is similar with finding desired spectra from the database, with addtional commands of bin size and parameters needed to be plotted.
    
    Input:
    
    1.  Number of queries.
    2.  Parameter that you want to plot a histogram for.
    3.  Bin size.
    4.  Upper bound of the histogram.
    5.  Lower bound of the histogram.
    
    6 and onwards.  SQL command to download the queries. Say if the number you input in argv[1] is 3, then you input 3 SQL commands to download the data.
    
    Input format for SQL queries:
    'SELECT * FROM Supernovae ORDER BY (the parameter you want to sort the spectra)  DESC LIMIT (how many spectra wanted)'
    
    The selection of param
    An example is:
    
    python newhistogram.py 2 Redshift 0.01 0.15 0 'SELECT * FROM Supernovae ORDER BY Phase DESC LIMIT 100' 'SELECT * FROM Supernovae ORDER BY Redshift ASC LIMIT 100'
    
    The code will plot a redshift histogram based on the data of 2 queries, with bin size of 0.01, upper bound = 0.15 and lower bound = 0. The two queries are respectively selecting the 100 SNe with max phase and 100 SNe with min redshift.
    
    The fields in the database are as follows
    Filename, SN name, Redshift, Phase, MinWave, MaxWave, Dm15, M_B, B_mMinusV_m
    
    Updated: Mar 14, 2014
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import sqlite3 as sq3

# Number of queries you want
que_num = int(sys.argv[1])

#connect to the database
con = sq3.connect('/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/MichaelSchubert/SNe.db')

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

    print 'There are ' + str(len(a[m])) + ' spectra with measured data'

bin = float(sys.argv[3])    ### Bin size
up = float(sys.argv[4])     ### Upper bound of the histogram
low = float(sys.argv[5])    ### Lower bound of the histogram
bin_num = int(np.ceil((up - low) / bin)) ### Number of bins

hist = [0] * que_num
bins = [0] * que_num
leg_name = [0] * que_num

colors = iter(cm.rainbow(np.linspace(0, 1, len(a))))

for m in range(que_num):
    hist[m], bins[m] = np.histogram(a[m], bins = bin_num)
    width = 0.7 * (bins[m][1] - bins[m][0])   ### To leave some space among stripes
    center = (bins[m][:-1] + bins[m][1:]) / 2
    plt.bar(center, hist[m], align='center', width = width, alpha = 0.3, color = next(colors))
    ### Insert the name of the legends
    leg_name[m] = sys.argv[m + 6].split()[6] + ' ' + sys.argv[m + 6].split()[7] + ' ' + sys.argv[m + 6].split()[9]

plt.ylabel(r'# of $SNe Ia$')
plt.xlabel(param)
plt.legend(leg_name)
plt.show()

