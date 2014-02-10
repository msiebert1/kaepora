'''
Created on Feb 9, 2014

@author: Sam Rubin
'''
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sqlite3 as sq3

#connect to the database
con = sq3.connect('../MichaelSchubert/SNe.db')
cur = con.cursor()

#checking to make sure everything is being read correctly
#for row in cur.execute('SELECT * FROM Supernovae ORDER BY redshift'):
#    print row

#Creates an array containing only the used information, which can be easily referenced 
data = []
for row in cur.execute('SELECT Filename, Redshift, MinWave, MaxWave FROM Supernovae'):
    data.append(row)
data = np.array(data)

#open the velocity file, make it useable
v_data = open('foley_master_data')
lines = v_data.readlines()
lines = lines[1:-1]
lines = np.array(lines)
v_data.close()
velocity = []
for line in lines:
    p = line.split()
    velocity.append(p)
velocity = np.array(velocity)
print velocity 

for line in velocity:
    name = line[0]
    v_si = line[2]
    dv_si = line[3]
    v_ca = line[6]
    dv_ca = line[7]

#create another database for velocity?
con2 = sq3.connect('velocity.db')
cur2 = con2.cursor()
cur2.execute('DROP TABLE IF EXISTS Velocity')
cur2.execute('CREATE TABLE IF NOT EXISTS Velocity (Filename TEXT, V_SiII REAL, DV_Si REAL, V_CaII REAL, DV_Ca REAL)')
cur2.execute('INSERT INTO Velocity(Filename, V_SiII, DV_Si, V_CaII, DV_Ca) VALUES (?, ?, ?, ?, ?)', (name, v_si, dv_si, v_ca, dv_ca))


con.close()