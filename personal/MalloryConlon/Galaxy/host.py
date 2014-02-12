import numpy as np
import os
import glob
import matplotlib.pyplot as plt

#Pre-allocate arrays
sn_name = []
host_type = []
host_name = []
elliptical = []
S0 = []
spiral = []
irregular = []
anon = []
spectra_files = []
good_files = []
bad_files = []
sn_subdir = []
sn_path = []

#Load in host galaxy information with morphology as reported in NED.
## Retrieved from http://www.cfa.harvard.edu/supernova/CfA3/snmorph.txt
## Morphology translation:
## 1,2,3,4,5,6,7,8,9,10,11 = E,E/S0,S0,S0a,Sa,Sab,Sb,Sbc,Sc,Scd,Sd/Irr
## 0 = No morphology available

host = np.genfromtxt("/users/malloryconlon/astr596/personal/malloryconlon/galaxy/host_info.dat", dtype = None, unpack = True)

for i in range(len(host)):
    load_host = host[i]
    sn_name.append(load_host[0])
    host_type.append(load_host[1])
    host_name.append(load_host[2])

for j in range(len(host)):
    if host_type[j]==1 or host_type[j]==2:
        elliptical.append(sn_name[j])
    if host_type[j]==3 or host_type[j]==4:
        S0.append(sn_name[j])
    if host_type[j]==5 or host_type[j]==6 or host_type[j]==7 or host_type[j]==8 or host_type[j]==9 or host_type[j]==10:
        spiral.append(sn_name[j])
    if host_type[j]==11:
        irregular.append(sn_name[j])
    if host_type[j]==0:
        anon.append(sn_name[j])

print elliptical
print len(S0)
print len(spiral)
print len(irregular)
print len(anon)


