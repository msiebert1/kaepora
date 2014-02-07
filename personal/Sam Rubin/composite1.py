import numpy as np
import matplotlib.pyplot as plt
import glob
import os

array = []
ignore = []
invalid = []
#cut out the portion of the wavelength that we want

def domain_cut(low, high, lines):
	if min(lines[:,0]) >= low and max(lines[:,0]) <= high:
		return "Out of Range"
	else:
		lines = np.delete(lines, np.where(lines[:,0] < low), axis=0)
		lines = lines[0:2000]
		return lines

#create a list of spectra files

files = glob.glob ('../../data/cfa/*/*.flm')


#reads valid files into composite array
for file in files:
	if file.endswith('.flm'):
		array.append(file)
	else: 
		ignore.append(file)
total = len(array)

#open files and read spectra data
composite_dict = {}
for i in range(500):
	try:
		atext = np.loadtxt(array[i-1])
		a = open(array[i-1])
		lines = domain_cut(3500,7500,atext)
		wavelength = []
		flux = []
		for line in lines:
			if type(lines) == str:
				break
			else:
				wavelength.append(line[0]) #create lists for each axis
				flux.append(line[1])
		wavelength = np.array(wavelength) #make it an array for plotting
		flux = np.array(flux)
		b = np.linspace(3500,7500,2001)
		flux = np.interp(b, wavelength, flux)
		combined = [b, flux]
		composite_dict[a.name] = combined #create dictionary for later use
		
	except ValueError:
		print "not a number"
#print composite_dict

#averaging spectra
avgflux = []
for i in range(len(b)):
	for key, value in composite_dict.items():
		try:
			x = []
			x.append(composite_dict[key[i]]) #this is still a problem
		except ValueError:
			print "also not a number"
		except KeyError:
			print "key does not exist"
	print x
	if len(x)==0:
		pass
	else:
		avgflux.append((sum(x)/len(x)))
	

#plot the averaged spectrum on a lin-log scale
p1,=plt.plot(b, avgflux)
plt.yscale('log')
plt.xlabel('Wavelength [A]')
plt.ylabel('Flux')
plt.savefig('composite_spectrum1.png')
plt.show()
