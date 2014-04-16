import composite
import Plotting
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys

composites=[]		
composites.append(composite.main("SELECT * FROM Supernovae WHERE Morphology=5 AND PHASE BETWEEN .78 AND .80"))
composites.append(composite.main("SELECT * FROM Supernovae WHERE Morphology=5 AND PHASE BETWEEN .80 AND .82"))
composites.append(composite.main("SELECT * FROM Supernovae WHERE Morphology=5 AND PHASE BETWEEN .82 AND .85"))


labels=['Phase: .78 to .82','Phase: .82 to .85','Phase: .85 to .89']

scales=composite.find_scales(composites,composites[0].flux,composites[0].ivar)

wmin = 0
wmax = 100000
for comp in composites:
	SN=comp
	if (SN.minwave > wmin):
		wmin=SN.minwave
	if (SN.maxwave < wmax):
		wmax=SN.maxwave

plt.figure()
i=0
for comp in composites:
	plt.plot(comp.wavelength,scales[i]*comp.flux,label=labels[i])
	i+=1

plt.title('SprialA')
plt.xlim(wmin,wmax)
legend = plt.legend(loc='upper right')
plt.show()	
	