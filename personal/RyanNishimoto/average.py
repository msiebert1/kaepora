import numpy as np
import matplotlib.pyplot as mp


#read spectra 
sn1 = np.loadtxt('sn2011by-hst+lick.flm')
sn2 = np.loadtxt('sn2011fe-visit3-hst.flm')

#create lists for wavelengths (wl) and fluxes (fx)
wl1 = sn1[:,0]
wl2 = sn2[:,0]
fx1 = sn1[:,1]
fx2 = sn2[:,1]

#trim data of spectra for comparison

#determines common starting wavelength and respective indexes
if wl1[0] > wl2[0]:
	start = wl1[0]
else:
	start = wl2[0]

first1 = np.where(wl1==start)
first2 = np.where(wl2==start)
#determines common ending wavelength and respective indexes
if wl1[len(wl1)-1] < wl2[len(wl2)-1]:
	end = wl1[len(wl1)-1]

else:
	end = wl2[len(wl2)-1]

last1 = np.where(wl1==end)
last2 = np.where(wl2==end)

#trim spectra
wl1 = wl1[first1[0]:last1[0]]
wl2 = wl2[first2[0]:last2[0]]
fx1 = fx1[first1[0]:last1[0]]
fx2 = fx2[first2[0]:last2[0]]

#average spectra
avgwl = (wl1+wl2)/2
avgfx = (fx1+fx2)/2

#plot spectras
mp.title("Compared Spectras of SN2011by and SN2011fe")
mp.xlabel("Wavelength")
mp.ylabel("Flux(scaled)")
mp.yscale('log')
mp.plot(wl1,fx1,'r')
mp.plot(wl2,fx2,'b')
mp.plot(avgwl,avgfx,'k')
mp.savefig("avgspectra.png")
mp.show()


