#This code looks for telluric lines in a given spectrum.  It clips those lines and changes the inverse variance to 0 because the flux value has been corrected.

def scale(root):

    import numpy as np
    from datafidelity import *
    import matplotlib.pyplot as plt
    import os
    list = []
    scales = []

    for path, subdirs, files in os.walk(root):
        for name in files:
            if name.endswith('.flm'):
                list.append(os.path.join(path,name))


    for i in range(len(list)):
        try:
            SN = np.genfromtxt(list[i])
            wavelength = SN[:,0]
            flux = SN[:,1]
            var = SN[:,2]
            var1 = genvar(wavelength,flux)
                #if np.sum(var)!= 0:
            new = 0
            new = var/var1
            scale = np.average(new)
            #print list[i]
            #print scale
            scales.append(scale)
            print len(scales)
        except:
            print list[i]

    new = [x for x in scales if str(x) != 'nan']

    new1 = [x for x in new if float(x) < 100]

    new2 = [x for x in new1 if float(x) != 0]

    new_scales = [x for x in new2 if str(x) != 'inf']

    print new_scales
    print np.average(new_scales)
    print np.median(new_scales)


    num_bins = 900
# the histogram of the data
    n, bins, patches = plt.hist(new_scales, num_bins, normed=1, facecolor='green', alpha=0.5)
    plt.xlabel('Scaling factors')
    plt.ylabel('Number of counts')
    plt.axis([0,100,0, 10])
    plt.show()
