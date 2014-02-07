#!/usr/bin/env python -i
# Load the data files and produce an average spectrum
#
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import interpolate
from matplotlib.ticker import AutoMinorLocator

def halfSearch(arr,find):
    median = 0
    bottom = 0
    ceil = len(arr) - 1
    while (bottom <= ceil):
        median = (bottom + ceil) / 2
        temp  = arr[median]
        if (temp == find):
            return median 
        else:
            if (temp < find):
                bottom = median + 1
            else:
                ceil = median - 1
    return -1 #not found

# Deredshifting wave  with redshift z
def dered(wave, z):
    wave /= 1. + z
    return wave

#Loading data from data files and do the deredshifting 
def loadSpec(file,z,xmax=False):
    f = open(file)
    lines = f.readlines()
    f.close()

    x = []
    y = []

    for line in lines:
        p = line.split()
        x.append(float(p[0]))
        y.append(float(p[1]))

    x = np.array(x)
    y = np.array(y)        
    
    dered(x,z)

    if xmax == True :
        return min(x),max(x)
    else:
        return x,y

    
def findBoundary(files,zhels):
    wavemin = 0.
    wavemax = 1e10
    
    for file,z in zip(files,zhels):
        wavemint, wavemaxt = loadSpec(file,z,xmax=True)
        if wavemint > wavemin:
            wavemin = wavemint
            if wavemaxt < wavemax:
                wavemax = wavemaxt
    
    wavemin = (int(wavemin/10) + 1) * 10
    wavemax = (int(wavemax/10) - 1) * 10
    return  np.array([wavemin,wavemax])


# Interpalate spectrum
def interpSpec(wave,flux,bound):
    nsample = int(1e5)
    intwave = scipy.linspace(bound[0],bound[1],nsample) 
    tck = interpolate.splrep(wave, flux)
    intflux = interpolate.splev(intwave,tck)
    return intwave,intflux    

# Get average spectrum
def averSpec(files,zhels):
    fluxs = []
    bound = findBoundary(files,zhels)

    for file,z in zip(files,zhels):
        wave,flux = loadSpec(file,z)
        intwave,intflux = interpSpec(wave,flux,bound)

        #Normalize flux
        intflux /= np.median(intflux)
        
        if len(fluxs) == 0:
            fluxs = np.array([intflux])
        else:
            fluxs = np.append(fluxs,np.array([intflux]),axis=0)
#             print np.shape(fluxs)
            
        averflux = np.mean(fluxs, axis=0)
        
        resflux = fluxs - averflux

    return intwave,averflux,fluxs,resflux

def plotSpec(wave,averflux,fluxs,resflux):

    pltdir = '../plots/'

    f, axarr = plt.subplots(2, sharex=True)

    ax1 = axarr[0]
    ax2 = axarr[1]

    minorLocator  = AutoMinorLocator(10)
    ax2.xaxis.set_minor_locator(minorLocator)

    minorLocator  = AutoMinorLocator(10)
    ax2.xaxis.set_minor_locator(minorLocator)

    ax2.tick_params(which='major', length=7)
    ax2.tick_params(which='minor', length=4)


    ax2.set_title('Wavelength [A]')
#     ax2.ylabel('Residule Flux')
#     ax1.ylabel('Scaled Flux')
    

    ax1.plot(wave, averflux,marker='o')
    ax2.plot(wave,np.zeros(len(wave)), marker='o')

    nspec = np.shape(fluxs)[0]
    index = range(0,nspec,1)

    for i in index:
        ax1.plot(wave,fluxs[i,:])
        ax2.plot(wave,resflux[i,:])

    plt.savefig(pltdir+'spectrum.eps')
    plt.show()



#     print 'flux:', fluxs[1,0:3]
#     print 'residule:',resflux[1,0:3]
#     print 'average:',averflux[0:3]

# datadir = '../../../data/'

# file1 = 'sn2011by-hst+lick.flm'
# file2 = 'sn2011fe-visit3-hst.flm'
# z1 = 0.003402
# z2 = 0.001208

# files = [datadir+file1,datadir+file2]
# zhels = [z1,z2]

def spectrum(files,zhels):
    wave,averflux,fluxs,resflux = averSpec(files,zhels)
    plotSpec(wave,averflux,fluxs,resflux)
    
    return


