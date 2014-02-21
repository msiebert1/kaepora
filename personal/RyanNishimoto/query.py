import os
import numpy as np
import astropy
from astropy.io import ascii
import astroquery
from astroquery.ned import Ned
from astroquery.irsa_dust import IrsaDust
from specutils.extinction import extinction
from specutils.extinction import extinction
import scipy

sn = "1994m"
ext = IrsaDust.get_extinction_table("SN%s"%sn)
for e in ext:
	print e

ascii.write(ext,"extinction.dat")

