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

class supernova(object):
    """Attributes can be added"""

con = sq3.connect('../../data/SNe.db')
cur = con.cursor()
i = 0

SN_Array = []
cur.execute('SELECT * FROM Supernovae WHERE Phase BETWEEN -3 AND 3 AND Velocity BETWEEN -12 AND -8')
for row in cur:
    SN = supernova()
    SN.filename = row[0]
    interp = msg.unpackb(row[15])
    SN.interp = interp
    print interp

print "done"
