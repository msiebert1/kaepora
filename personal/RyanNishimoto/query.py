import os
import numpy as np
import astropy
from astropy.io import ascii
from astroquery.ned import Ned
from specutils.extinction import extinction
from specutils.extinction import extinction
import scipy

sn = "1994m"
photometry = Ned.get_table("SN%s"%sn,table='photometry')
for p in photometry:
	print p

ascii.write(photometry,"table.dat")
