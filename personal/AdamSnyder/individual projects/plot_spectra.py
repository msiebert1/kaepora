import numpy as np
import matplotlib.pyplot as plt

# function two average 2 arbitrary spectra
# Plan to expand to average arbitrary numbers of spectra, need better way to average
def average_spectra(spectra1, spectra2):
    if len(spectra1[:,0]) < len(spectra2[:, 0]):
        S1 = spectra2
        S2 = spectra1

    else:
        S1 = spectra1
        S2 = spectra2
    j = 0
    S3 = []

    for i in range(len(S1[:,0])):
        if S1[i,0] == S2[j,0]:
            S3.append([S1[i, 0], (S1[i,1]+S2[j, 1])/2.0])
            if j < len(S2)-1: j+= 1

    S3 = np.array(S3, dtype=float)
    return S3
            

# Load data from files

SN1 = np.loadtxt("../../data/sn2011by-hst+lick.flm")
SN2 = np.loadtxt("../../data/sn2011fe-visit3-hst.flm")

# Plot both spectra on one graph

f1 = plt.figure()
plt.plot(SN1[:,0], SN1[:,1], SN2[:,0], SN2[:,1])
plt.ylabel("Flux")
plt.xlabel("Wavelength")

# Construct an "average" spectrum graph

SN3 = average_spectra(SN1, SN2)

np.savetxt("averagespectra.flm", SN3, fmt=["%.2f", "%.5e"])

f2 = plt.figure()
plt.plot(SN3[:, 0], SN3[:, 1])
plt.ylabel("Flux")
plt.xlabel("Wavelength")
plt.show()