The files ending in ".dat" in this distribution are the optical spectra
of Type Ia supernovae published by the Carnegie Supernova Project (CSP) 
in Folatelli et al. (2013). 

The file names contain the SN name, observing UT date (yymmdd), wavelength 
range indicator, telescope tag, and instrument tag. Wavelength range 
indicators are "b01", "g01" and "r01" and are used to distinguish spectra 
obtained on the same night and with the same instrument but with different 
setups.

The file contents include a seven-line header with information about the 
supernova and spectrum. The header includes: the name of the supernova, the 
original FITS file spectrum in the CSP database, the Heliocentric redshift, 
the Julian date of maximum light in the B band, the Julian date of the 
observation, and the resulting epoch of the spectrum in rest-frame days since 
maximum light.

The data themselves are given in two columns for rest-frame wavelength in 
angstrom, and flux density in erg/s/cm^2/AA.

WAVELENGTH IS EXPRESSED IN THE REST FRAME OF THE SUPERNOVA IN ALL CASES. 
Observed wavelength has been divided by (1+z), where z is the (Heliocentric) 
redshift provided in the spectrum header. 

Spectra have been corrected for telluric absorption, except in the cases 
indicated in Table 2 of Folatelli et al. (2013). Some of the spectra, as 
indicated in the same Table, have been corrected to match the observed 
photometry of the supernova.
