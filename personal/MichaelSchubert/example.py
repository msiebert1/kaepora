import sqlite3 as sq3
import msgpack
import msgpack_numpy as mn
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import time
mn.patch()
start = time.clock()
"""

f = open('sn2011by-hst+lick.flm')
data = np.loadtxt(f)
minwave = min(data[:, 0])
maxwave = max(data[:, 0])
shift = .05

"""
con = sq3.connect('SNe.db')
cur = con.cursor()

#cur.execute('SELECT * FROM Supernovae WHERE MinWave > 4000')

cur.execute('SELECT Phase FROM Supernovae ORDER BY M_B DESC LIMIT 1')
for row in cur:
    """
    shift = row[2]
    minwave = row[4]
    maxwave = row[5]
    data = msgpack.unpackb(row[len(row)-1])
    waves = data[:, 0]
    flux = data[:, 1]
    shifted = np.divide(data[:, 0], (shift))
    """
"""
x_val = np.linspace(minwave, maxwave, 10000)
f = interp1d(waves, flux)
print f(shifted)

x = waves
y = flux
f = interp1d(x, y)
#xnew = np.linspace(minwave, maxwave, 10)
xnew = np.linspace(minwave, 1000, 4000)
interpd = f(xnew)



plt.plot(x, y, 'r', xnew, f(xnew), 'b')
plt.show()
"""
end = time.clock()
print end-start