import datafidelity as df

def find_extrema(SN):
	wavelength = SN.wavelength[SN.x1:SN.x2]
	flux = SN.flux[SN.x1:SN.x2]
	varflux = 1./SN.ivar[SN.x1:SN.x2]
	df = gsmooth(wavelength, flux, varflux, .002, nsig)