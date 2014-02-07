import numpy as np

data = {}
#sn_data = np.genfromtxt("cfasnIa_param.dat", dtype=None, skip_header=42, names=True)
sn_data = np.loadtxt("new_params.dat")
arr = sn_data[0]
print arr
#data = arr[1:]
print len(arr)
print len(data)
#print sn_data.dtype