import numpy as np
import os, fnmatch
import time

def get_spectra(filename):
    root = '../../data/'
    for path, subdirs, files in os.walk(root):
        for filename in fnmatch.filter(files, filename):
            return os.path.join(path, filename)

def deredshift(spectra, redshift):
    """
    Divides wavelength of given spectra by 1 + redshift
    """
    wavelengths = spectra[:, 0]
    shifted = np.divide(wavelengths, (1 + z))
    return shifted

def interpolate():
    return
def find_error(spectra):
    return 1/((spectra[:, 2])**2)

start = time.clock()
loc = get_spectra('snf20080720-001-20080724.44-fast.flm')
f = np.loadtxt(loc)
print f
e = find_error(f)
print e
end = time.clock()
print end-start

