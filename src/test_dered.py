# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 15:14:33 2016

@author: siebe_000
"""
import numpy as np
import matplotlib.pyplot as plt
# from specutils import extinction as ex
from specutils import Spectrum1D
from dust_extinction.parameter_averages import F99
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
        
def dered(sne, snname, wave, flux, ivar, source='not_swift_uv', ):
    corrected = False
    for j in range(len(sne)):  # go through list of SN parameters
        sn = sne[j][0]
        if (snname.lower() == sn.lower()[2:]) or (snname.lower() == sn.lower()) or (snname.lower()[2:] == sn.lower()):  # SN with parameter matches the path
            if not corrected:
                # print 'Milky Way correction...'
                if source != 'swift_uv' and source != 'foley_hst' and source != 'foundation' and source != 'bsnip2'  and source != 'kyleplot' and source != 'marion09':
                    b = sne[j][1].astype(float)
                    v = sne[j][2].astype(float)
                    bv = b-v
                    # red = ex.reddening(wave, a_v=3.1*bv, r_v=3.1, model='f99')

                    ext = F99(Rv=3.1)
                    red = ext.extinguish(wave*u.AA, Av = 3.1*bv)

                    flux *= red
                    ivar *= 1./(red**2.)
                    corrected = True
                else:
                    Av = sne[j][1].astype(float)
                    # red = ex.reddening(wave, a_v=Av, r_v=3.1, model='f99')
                    ext = F99(Rv=3.1)
                    red = ext.extinguish(wave*u.AA, Av = Av)

                    flux *= red
                    ivar *= 1./(red**2.)
                    corrected = True

    if not corrected:
        print ('No MW extinction info for', snname)
        raise TypeError   
    return flux, ivar

def host_correction(a_v, r_v, snname, wave, flux, ivar, model = 'f99'):
    # print 'Host correction...'
    # red = ex.reddening(wave, a_v = a_v, r_v = r_v, model=model)
    ext = F99(Rv=r_v)
    red = ext.extinguish(wave*u.AA, Av = a_v)
    flux *= red
    ivar *= 1./(red**2.) #correct ivar too
    return flux, ivar




