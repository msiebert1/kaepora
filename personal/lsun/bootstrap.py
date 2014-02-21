# Authors: Ricky, Lunan
# Input number of tries (argv[1])

import numpy as np
import re, sys

filename = open('name.txt', 'r').read()

list = filename.split('\n') # Name of the files.
num = len(list)             # Number of samples.

spec = [0] * num

for i in range(num):
    spec[i] = np.loadtxt('sn1995al/' + list[i], unpack = True, usecols = (0, 1,))

num_arr = np.arange(0, num, 1)

tries = (int)(sys.argv[1])                 # Number of tries.

sel_spec = [0] * tries      # An array for the bootstraped spectra.

for i in range(tries):
    sel_spec[i] = np.round(np.random.uniform(0, num, num))

for i in range(tries):
    print sel_spec[i]
#for m in range(tries):
