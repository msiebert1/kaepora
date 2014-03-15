import sqlite3 as sq3

#connect to the database
con = sq3.connect('/Users/rickyccy/Documents/Urbana-Champaign/Courses/ASTR596_Spring2014/astr596/personal/RickyChue/personal/MichaelSchubert/SNe.db')

#create a cursor object to execute commands
cur = con.cursor()

#The fields in the database are as follows
#Filename, SN name, Redshift, Phase, MinWave, MaxWave, Dm15, M_B, B_mMinusV_m

#Suppose you want to return the 10 SNe with the highest redshifts
cur.execute('SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 10')

#The output is stored in the cursor object.
for row in cur:
     #print sn filename
     print row[0], row[1], row[2]
