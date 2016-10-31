# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 15:14:33 2016

@author: siebe_000
"""
import numpy as np
import prep
import matplotlib.pyplot as plt
from specutils import extinction as ex
from specutils import Spectrum1D
from astropy import units as u


def find_SN(fname, source=None, csplist=None):
    """
    Returns SN name, either from cspdict if source is a csp spectra
    or from slicing the file name if source is Cfa or bsnip
    """
    if source == 'csp':
        snname = csplist[0]
        return snname[2:]
    elif source == 'other':
        snname = fname.replace('_', '-').split('-')
        if snname[0][:2] == 'sn':
            return snname[0][2:]
        else:
            return snname[0]
    else:
        snname = fname.replace('_', '-').split('-')
        if snname[0][:3] == 'snf':
            namelist = [snname[0], snname[1]]
            snname = '-'.join(namelist).upper()
        else:
            snname = snname[0][2:]

        return snname
        
def dered(sne, snname, wave, flux):
    for j in range(len(sne)):  # go through list of SN parameters
        sn = sne[j][0]
        if sn in snname:  # SN with parameter matches the path
            b = sne[j][1].astype(float)
            v = sne[j][2].astype(float)
            bv = b-v
#            print "B(%s)-V(%s)=%s"%(b,v,bv)
#            print "R(v) =",r
            #or use fm07 model
            #test1 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='ccm89')
            #test2 = spectra_data[i][:,1] * ex.reddening(spectra_data[i][:,0],ebv = bv, model='od94')
            red = ex.reddening(wave, a_v=1.33*1.8, r_v=1.8, model='f99')
            flux *= red
            # flux *= ex.reddening(wave, a_v=bv, r_v=3.1, model='f99')

#            wave /= (1+z)

            #print "de-reddened by host galaxy\n",flux*ex.reddening(wave,ebv = 0, r_v = r, model='f99')
            #host *= ex.reddening(wave,ebv = bv, r_v = r, model='f99')

    return flux
        
fname = '..\data\spectra\cfa\sn2003cg\sn2003cg-20030331.21-fast.flm'
spectrum  = np.loadtxt(fname)
sn_name = find_SN(fname)
old_wave = spectrum[:, 0]*u.Angstrom	    # wavelengths
old_flux = spectrum[:, 1]*u.Unit('W m-2 angstrom-1 sr-1')
spec1d = Spectrum1D.from_array(old_wave, old_flux)

sne = prep.ReadExtin('extinction.dat')
plt.plot(old_wave,old_flux)
plt.show()
new_flux = dered(sne, sn_name, spec1d.wavelength, spec1d.flux)
#new_flux = dered(sne, sn_name, old_wave, old_flux)
plt.plot(old_wave, new_flux)
plt.show()