import glob
import os


for path,subdirs,files in os.walk('../../../../data/cfa/'):
    


for i in range(len(list)):
    SN = np.genfromtxt('../../../data/cfa/'+list[i])
    wavelength = SN[:,0]
    flux = SN[:,1]
    var = SN[:,2]
    var1 = genvar(wavelength,flux)
    plt.plot(wavelength,var)
    plt.plot(wavelength,var1)
    plt.show()