"""
clean spectra

then put into database

selection
"""
#read in all spectra
#pixels

#teaks for composite
#de-reddening
	#dust in milky way
		#using fitzpatrick law
		#find milky way reddening for particular SN (look on NED)

	#de-redshifting

	#by host galaxy (assume zero for now)
		#color of SN matters

#interpolate
#NO TRUNCATING
#1000-20000
#set inverse variance (square of error)(weight) to zero outside of wavelength ranges 

#List of Files containing Spectra Data
spectra_files = glob.glob ('../../data/cfa/*/*.flm')	#array of files

spectra_data = []
file_path = []

junk_data = []

for file in spectra_files:
	try:
		spectra_data.append(np.loadtxt(file))
		file_path.append(file[14:-4])
	except ValueError:
		junk_data.append(file)
	
#table containing sn names, redshifts, etc.
sn_features = np.genfromtxt('../../data/cfa/cfasnIa_param.dat',dtype = None)	#stats of SNe



