import numpy as np
import matplotlib.pyplot as plt

SN1 = np.loadtxt("../../data/sn2011by-hst+lick.flm")
SN2 = np.loadtxt("../../data/sn2011fe-visit3-hst.flm")

f1 = plt.figure()
plt.plot(SN1[:,0], SN1[:,1], SN2[:,0], SN2[:,1])
SN3 = []

for i in range(len(SN1[:,0])):
    for j in range(len(SN2[:,0])):
        if SN1[i,0] == SN2[j,0]:
            SN3.append([SN1[i, 0], (SN1[i,1]+SN2[j, 1])/2.0])

SN3 = np.array(SN3, dtype=float)
np.savetxt("averagespectra.flm", SN3, fmt=["%.2f", "%.5e"], delimiter=",")

f2 = plt.figure()
plt.plot(SN3[:, 0], SN3[:, 1])
plt.show()
    
        



    
