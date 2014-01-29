##This code finds the average spectrum for two provided spectra.  Work is being done such that it can have an input directory with n spectra. 

import numpy as np
import matplotlib.pyplot as plt

start=0
#Read in files.  This will need to be automated once we are looking at a ton of spectra...
file1_red=np.loadtxt('sn2011by-hst+lick.flm')
file2_red=np.loadtxt('sn2011fe-visit3-hst.flm')

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
medwave1=np.median(wave1)
medwave2=np.median(wave2)
medflux1=np.median(flux1)
medflux2=np.median(flux2)



#Average the two spectra
average_wave=(wave1+wave2)/2
average_flux=(flux1+flux2)/2

#Make plots of individual and average spectra, save plot
plt.yscale('log')
plot1,=plt.plot(wave1,flux1,'k')
plot2,=plt.plot(wave2,flux2, 'g')
plot3,=plt.plot(average_wave,average_flux, 'r')

plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('log[Scaled Flux]')
plt.legend([plot1,plot2,plot3],['SN2011BY','SN2011FE','Average Spectrum'],4)
plt.savefig('../personal/malloryconlon/'+'spectrum.pdf')
plt.show()

#Write average spectrum to file...need to think more about how to store this properly.

#np.savetxt('../personal/malloryconlon/'+'average.flm',(average_wave,average_flux))