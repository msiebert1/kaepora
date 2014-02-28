import os
import glob
from specutils import extinction as ex
from astropy.table import Table
from astropy.io import ascii
import astroquery
from astroquery.irsa_dust import IrsaDust
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter


#list of files
spectra_files = glob.glob ('../../data/cfa/*/*.flm')

#holds spectra data (wavelength,flux,weight)
spectra_data = []
#holds file pathname
file_path = []

junk_data = []

#number of spectra to modify
num = 1

#get data, pathnames
for i in range(num):
	try:
        	spectra_data.append(np.loadtxt(spectra_files[i]))
        	file_path.append(spectra_files[i][14:-4])
		#print file_path
             
	except ValueError:
		junk_data.append(spectra_files)

#update num to number of good spectra files
num = len(spectra_data)

#table containing sn names, redshifts, etc.
sn_parameters = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)

#holds sn name
sn = []
#holds redshift value
z = []

#get relevent parameters needed for calculations
for i in range(len(sn_parameters)):
	sn.append(sn_parameters[i][0])
	z.append(sn_parameters[i][1])

"""
NOTE:
Use NED
"""

"""
Note:
create separate code that pulls all the data and creates a "parameter" table much like the one used in past,
but it pulls all relevant data via astroquery and writes it in a table.  
must do that following:
-read in all files we want to prep
-truncate file name to format "SNyear"(this is how astroquery searches for SN)
-get relevent data into table, following a format like:
SN name		Host Galaxy		Redshift	B	V	De-redden factor	CarbonPos/Neg
"""
def fm_unred(wave, flux, ebv, *args, **kwargs):
    
    '''
    NAME:
     FM_UNRED
    PURPOSE:
     Deredden a flux vector using the Fitzpatrick (1999) parameterization
    EXPLANATION:
     The R-dependent Galactic extinction curve is that of Fitzpatrick & Massa 
     (Fitzpatrick, 1999, PASP, 111, 63; astro-ph/9809387 ).    
     Parameterization is valid from the IR to the far-UV (3.5 microns to 0.1 
     microns).  UV extinction curve is extrapolated down to 912 Angstroms.

    CALLING SEQUENCE:
     fm_unred( wave, flux, ebv [, 'LMC2', 'AVGLMC', 'ExtCurve', R_V = ,  
                                   gamma =, x0=, c1=, c2=, c3=, c4= ])
    INPUT:
      wave - wavelength vector (Angstroms)
      flux - calibrated flux vector, same number of elements as "wave"
      ebv  - color excess E(B-V), scalar.  If a negative "ebv" is supplied,
              then fluxes will be reddened rather than dereddened.

    OUTPUT:
      Unreddened flux vector, same units and number of elements as "flux"

    OPTIONAL INPUT KEYWORDS
      R_V - scalar specifying the ratio of total to selective extinction
               R(V) = A(V) / E(B - V).  If not specified, then R = 3.1
               Extreme values of R(V) range from 2.3 to 5.3

      'AVGLMC' - if set, then the default fit parameters c1,c2,c3,c4,gamma,x0 
             are set to the average values determined for reddening in the 
             general Large Magellanic Cloud (LMC) field by Misselt et al. 
            (1999, ApJ, 515, 128)
      'LMC2' - if set, then the fit parameters are set to the values determined
             for the LMC2 field (including 30 Dor) by Misselt et al.
             Note that neither /AVGLMC or /LMC2 will alter the default value 
             of R_V which is poorly known for the LMC. 
             
      The following five input keyword parameters allow the user to customize
      the adopted extinction curve.  For example, see Clayton et al. (2003,
      ApJ, 588, 871) for examples of these parameters in different interstellar
      environments.

      x0 - Centroid of 2200 A bump in microns (default = 4.596)
      gamma - Width of 2200 A bump in microns (default = 0.99)
      c3 - Strength of the 2200 A bump (default = 3.23)
      c4 - FUV curvature (default = 0.41)
      c2 - Slope of the linear UV extinction component 
           (default = -0.824 + 4.717 / R)
      c1 - Intercept of the linear UV extinction component 
           (default = 2.030 - 3.007 * c2)
            
    OPTIONAL OUTPUT KEYWORD:
      'ExtCurve' - If this keyword is set, fm_unred will return two arrays.
                  First array is the unreddend flux vector.  Second array is
                  the E(wave-V)/E(B-V) extinction curve, interpolated onto the
                  input wavelength vector.

    EXAMPLE:
       Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
       is altered by a reddening of E(B-V) = 0.1.  Assume an "average"
       reddening for the diffuse interstellar medium (R(V) = 3.1)

       >>> w = 1200 + arange(40)*50       #Create a wavelength vector
       >>> f = w*0 + 1                    #Create a "flat" flux vector
       >>> fnew = fm_unred(w, f, -0.1)    #Redden (negative E(B-V)) flux vector
       >>> plot(w, fnew)                   

    NOTES:
       (1) The following comparisons between the FM curve and that of Cardelli, 
           Clayton, & Mathis (1989), (see ccm_unred.pro):
           (a) - In the UV, the FM and CCM curves are similar for R < 4.0, but
                 diverge for larger R
           (b) - In the optical region, the FM more closely matches the
                 monochromatic extinction, especially near the R band.
       (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
               can be represented with the FM curve, if the proper value of 
               R(V) is supplied.
    REQUIRED MODULES:
       scipy, numpy
    REVISION HISTORY:
       Written   W. Landsman        Raytheon  STX   October, 1998
       Based on FMRCurve by E. Fitzpatrick (Villanova)
       Added /LMC2 and /AVGLMC keywords,  W. Landsman   August 2000
       Added ExtCurve keyword, J. Wm. Parker   August 2000
       Assume since V5.4 use COMPLEMENT to WHERE  W. Landsman April 2006
       Ported to Python, C. Theissen August 2012
    '''
    
    # Import needed modules
    from scipy.interpolate import InterpolatedUnivariateSpline as spline
    import numpy as n

    # Set defaults
    lmc2_set, avglmc_set, extcurve_set = None, None, None
    R_V, gamma, x0, c1, c2, c3, c4 = None, None, None, None, None, None, None
    
    x = 10000. / n.array([wave])                # Convert to inverse microns
    curve = x * 0.

    # Read in keywords
    for arg in args:
        if arg.lower() == 'lmc2': lmc2_set = 1
        if arg.lower() == 'avglmc': avglmc_set = 1
        if arg.lower() == 'extcurve': extcurve_set = 1
        
    for key in kwargs:
        if key.lower() == 'r_v':
            R_V = kwargs[key]
        if key.lower() == 'x0':
            x0 = kwargs[key]
        if key.lower() == 'gamma':
            gamma = kwargs[key]
        if key.lower() == 'c4':
            c4 = kwargs[key]
        if key.lower() == 'c3':
            c3 = kwargs[key]
        if key.lower() == 'c2':
            c2 = kwargs[key]
        if key.lower() == 'c1':
            c1 = kwargs[key]

    if R_V == None: R_V = 3.1

    if lmc2_set == 1:
        if x0 == None: x0 = 4.626
        if gamma == None: gamma =  1.05	
        if c4 == None: c4 = 0.42   
        if c3 == None: c3 = 1.92	
        if c2 == None: c2 = 1.31
        if c1 == None: c1 = -2.16
    elif avglmc_set == 1:
        if x0 == None: x0 = 4.596  
        if gamma == None: gamma = 0.91
        if c4 == None: c4 = 0.64  
        if c3 == None: c3 =  2.73	
        if c2 == None: c2 = 1.11
        if c1 == None: c1 = -1.28
    else:
        if x0 == None: x0 = 4.596  
        if gamma == None: gamma = 0.99
        if c4 == None: c4 = 0.41
        if c3 == None: c3 =  3.23	
        if c2 == None: c2 = -0.824 + 4.717 / R_V
        if c1 == None: c1 = 2.030 - 3.007 * c2
    
    # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
    # R-dependent coefficients
 
    xcutuv = 10000.0 / 2700.0
    xspluv = 10000.0 / n.array([2700.0, 2600.0])
   
    iuv = n.where(x >= xcutuv)
    iuv_comp = n.where(x < xcutuv)

    if len(x[iuv]) > 0: xuv = n.concatenate( (xspluv, x[iuv]) )
    else: xuv = xspluv.copy()

    yuv = c1  + c2 * xuv
    yuv = yuv + c3 * xuv**2 / ( ( xuv**2 - x0**2 )**2 + ( xuv * gamma )**2 )

    filter1 = xuv.copy()
    filter1[n.where(xuv <= 5.9)] = 5.9
    
    yuv = yuv + c4 * ( 0.5392 * ( filter1 - 5.9 )**2 + 0.05644 * ( filter1 - 5.9 )**3 )
    yuv = yuv + R_V
    yspluv = yuv[0:2].copy()                  # save spline points
    
    if len(x[iuv]) > 0: curve[iuv] = yuv[2:len(yuv)]      # remove spline points

    # Compute optical portion of A(lambda)/E(B-V) curve
    # using cubic spline anchored in UV, optical, and IR

    xsplopir = n.concatenate(([0], 10000.0 / n.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])))
    ysplir   = n.array([0.0, 0.26469, 0.82925]) * R_V / 3.1
    ysplop   = [n.polyval(n.array([2.13572e-04, 1.00270, -4.22809e-01]), R_V ), 
                n.polyval(n.array([-7.35778e-05, 1.00216, -5.13540e-02]), R_V ),
                n.polyval(n.array([-3.32598e-05, 1.00184, 7.00127e-01]), R_V ),
                n.polyval(n.array([-4.45636e-05, 7.97809e-04, -5.46959e-03, 1.01707, 1.19456] ), R_V ) ]
    
    ysplopir = n.concatenate( (ysplir, ysplop) )
    
    if len(iuv_comp) > 0:
        cubic = spline(n.concatenate( (xsplopir,xspluv) ), n.concatenate( (ysplopir,yspluv) ), k=3)
        curve[iuv_comp] = cubic( x[iuv_comp] )

    # Now apply extinction correction to input flux vector
    curve = ebv * curve[0]
    flux = flux * 10.**(0.4 * curve)
    if extcurve_set == None: return flux
    else:
        ExtCurve = curve - R_V
        return flux, ExtCurve


#deredden and deredshift the spectra
for i in range(num):#go through selected spectra data
	for j in range(len(sn)):#go through list of SN parameters
		if sn[j] in file_path[i]:#SN with parameter matches the path
			print "looking at SN",sn[j]
			ext = IrsaDust.get_extinction_table("SN%s"%sn[j])
			print ext[1][0],ext[1][3]
			print ext[2][0],ext[2][3]
			b = ext[1][3]
			v = ext[2][3]
							
			print "starting flux:\n",spectra_data[i][:,1]
			print "de-reddened with specutils",spectra_data[i][:,1]*ex.reddening(spectra_data[i][:,0],ebv=b-v,r_v=3.1,model='fm07')
			spectra_data[i][:,1] = fm_unred(spectra_data[i][:,0],spectra_data[i][:,1],b-v,R_V=3.1)
			print "de-reddened flux with fm_unred:\n",spectra_data[i][:,1]
			#print "starting wavelength:\n",spectra_data[i][:,0]
			spectra_data[i][:,0] /= (1+z[j])
			#print "z:",z[j]
			#print "de-red-shifted wavelength:\n",spectra_data[i][:,0]	

