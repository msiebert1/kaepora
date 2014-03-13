#!/usr/bin/env python -i 
# Provide a function to return the SNs index for a SQL query.
#
import sqlite3 as sq3

con = sq3.connect('../../MichaelSchubert/SNe.db')

cur = con.cursor()

cur.execute('SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 10')

for row in cur:
     #print sn filename
     print row
