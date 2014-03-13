#!/usr/bin/env python -i 
# Provide a function to return the SNs index for a SQL query.
#
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn

mn.patch()

def selectsn(sndb,sqlstr):
     con = sq3.connect(sndb)

     cur = con.cursor()
     
     cur.execute(sqlstr)
     
     return cur

# Calling example:
# from selectsn import selectsn
# sndb = '../../MichaelSchubert/SNe.db'
# sqlstr = 'SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 10'
# cur = selectsn(sndb,sqlstr)

# for row in cur:
#      print row
#      spectra = msg.unpackb(row[5])

# print spectra
