import sqlite3 as sq3
import msgpack
import msgpack_numpy as mn

mn.patch()

con = sq3.connect('SNe.db')
cur = con.cursor()

#cur.execute('SELECT * FROM Supernovae WHERE MinWave > 4000')
cur.execute('SELECT * FROM Supernovae ORDER BY M_B DESC LIMIT 1')
for row in cur:
    print msgpack.unpackb((row[len(row) - 1]))