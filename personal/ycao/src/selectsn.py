#!/usr/bin/env python -i 
# Provide a function to return the SNs index for a SQL query.
#
import numpy as np
import sqlite3 as sq3
import msgpack as msg
import msgpack_numpy as mn

mn.patch()

class supernova(object):
     """Attributes can be added"""

def selectsn(sndb,sqlstr):

     con = sq3.connect(sndb)
     cur = con.cursor()
     
     cur.execute(sqlstr)

     SN_Array = []
     names = []
     for row in cur:
          SN = supernova()
          SN.filename = row[0]
          SN.name = row[1]
          SN.redshift = row[2]
          SN.minwave = row[3]
#           spectra = msg.unpackb(row[5])
#           SN.spectrum = spectra
          SN_Array.append(SN)
          names.append(SN.name)

     names, indices = np.unique(names, return_index=True)
     print indices, names
     SN_Array = SN_Array[indices]
     return SN_Array

# Calling example:
# from selectsn import selectsn
sndb = '../../MichaelSchubert/SNe.db'
sqlstr = 'SELECT * FROM Supernovae ORDER BY Redshift DESC LIMIT 10'
SN_Array = selectsn(sndb,sqlstr)

# for row in cur:
#      print row
#      spectra = msg.unpackb(row[5])

# print spectra
