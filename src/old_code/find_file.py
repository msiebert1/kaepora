import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table
import msgpack as msg
import msgpack_numpy as mn
from scipy.optimize import leastsq
import file_name

mn.patch()

con = sq3.connect('../data/SNe.db')
cur = con.cursor()

cur.execute("SELECT Filename FROM Supernovae WHERE Morphology=5 AND PHASE BETWEEN .82 AND .85")
sp1=[]
for row in cur:
	sp1.append(row[0])
"""
cur.execute("SELECT Filename FROM Supernovae WHERE Morphology=4 AND PHASE BETWEEN -.5 AND 0")
sp2=[]
for row in cur:
	sp2.append(row[0])
cur.execute("SELECT Filename FROM Supernovae WHERE Morphology=4 AND PHASE BETWEEN 0 AND 1")
sp3=[]
for row in cur:
	sp3.append(row[0])
	
	"""

print sp1