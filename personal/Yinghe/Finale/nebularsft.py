# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 14:25:16 2014

@author: QuantumMonkey
"""

# Modified: import the nebular shift data to select sample spectra

from astropy.table import Table

# snsample,vnebsample,zsample,files

def vbin(vmin,vmax):
    nfile=[]
    nz = []
    nv = []
    nname = []
    for file,z,vneb,name in zip(files, zsample,vnebsample,snsample):
        if vneb >= vmin and vneb < vmax:
           nfile.append(file)
           nz.append(z)
           nv.append(vneb)
           nname.append(name)
    return nfile,nz,nv,nname

def selectneb():

    param = Table.read('nebular.dat')
    sn_name = param['col1']
    vneb = param['col2']
    verror = param['col3']
    

            


