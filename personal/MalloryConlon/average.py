##This code finds the average spectrum for two provided spectra.  The required input is a list of files in a particular directory that you wish to be averaged.

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

#Read in files.  This will need to be automated once we are looking at a ton of spectra...

for dirs,subdirs,files in os.walk('../data/'):
    for subdirs in dirs:
        list= glob.glob("*.flm")
print list

file1_red=np.loadtxt(list[0])
file2_red=np.loadtxt(list[1])

wave=file1_red[:,0]
wave2=file2_red[:,0]
flux=file1_red[:,1]
flux2=file2_red[:,1]

#Find the first common minimum wavelength of the two spectra
for i in range(len(wave)):
    if wave[i]==wave2[0]:
        start=i

end=start+len(flux2)

#Truncate the spectra such that they have the same range of wavelengths
#Working on generalizing this...
wave1=[]
flux1=[]
wave1=wave[start:end]
flux1=flux[start:end]


#Uncomment the next 2 lines to de-redden the spectra.  Need to consider how to get redshifts automatically without looking up for each SNe

#wave1=wave1/(1+0.003402)
#wave2=wave2/(1+0.001208)

#Scale the two spectra such that they are comparable.  Again, this will need to be automated...
medflux1=np.median(flux1)
medflux2=np.median(flux2)

flux1=flux1/medflux1
flux2=flux2/medflux2

#Average the two spectra
average_wave=(wave1+wave2)/(len(list))
average_flux=(flux1+flux2)/(len(list))

#Make plots of individual and average spectra, save plot
plt.yscale('log')
plot1,=plt.plot(wave1,flux1,'k')
plot2,=plt.plot(wave2,flux2, 'g')
plot3,=plt.plot(average_wave,average_flux, 'r')

plt.xlabel('Observed Wavelength ($\AA$)')
plt.ylabel('log[Scaled Flux]')
plt.legend([plot1,plot2,plot3],['SN2011BY','SN2011FE','Average Spectrum'],4)
plt.savefig('../personal/malloryconlon/'+'spectrum.pdf')
plt.show()

#Write average spectrum to file...need to think more about how to store this properly.

#np.savetxt('../personal/malloryconlon/'+'average.flm',(average_wave,average_flux))