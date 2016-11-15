import matplotlib.pyplot as plt
import numpy as np

# file =  '..\data\spectra\cfa\sn2002g\sn2002g-20020211-ui.flm'
file =  '..\data\spectra\/bsnip\sn2002g-20020211-ui.flm'
sn = np.genfromtxt(file, dtype=None)

wave = []
flux = []
for i in range(len(sn)):
	wave.append(float(sn[i][0]))
	flux.append(float(sn[i][1]))

plt.plot(wave,flux)
plt.show()