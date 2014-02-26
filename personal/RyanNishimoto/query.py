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

sn = "2007kg"
ext = IrsaDust.get_extinction_table("SN%s"%sn)
for i in range(len(ext)):
	if "B" in ext[i][0]:
		print ext[i][3]
	if "V" in ext[i][0]:
		print ext[i][3]

print ext
