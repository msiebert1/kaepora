import composite
import Plotting
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
import galrun
#import targeted
"""
Here's the main function that anyone can use to run some analysis.

The first line runs the composite program.
This grabs the selected data from SNe.db.
Change the input to composite.main to whatever you want your database query to be.
This is the code that outputs a text file of data for later use.
It also shows/saves a basic plot of the composite data and associated errors.

The next section sets up the plots to be used, and then plots using Plotting.py
Right now it doesn't work...but it will

More can be added as our programs evolve

Here are a few parameters you should set beforehand.
plot_name is where the plot showing both composite spectra together will be saved.
wmin and wmax define the boundary of the plot.
"""



#This part works just fine
#composite_full = composite.main("SELECT * FROM Supernovae WHERE Redshift > 0 AND Phase > -100")
"""
Here we set the queries that get used to find the spectra for compositing.
We only want to select spectra that have data for both redshift and phase, so both of them need to be in the query.
But you can change the values to whatever you want, and add more parameters.
"""
def run():
	names_array=[]
	name='Composite_comparison'
	morph=input('Morphologies: (E=1,E/S0=2,S0=3,S0a=4,Sa=5,Sab=6,Sb=7,Sbc=8,Sc=9,Scd=10,Sd/Irr=11)')
	labels=[]
	if 1 in morph:
		labels.append('Elliptical')
		name+=',E'
	if 2 in morph:
		labels.append('Elliptical/S0')
		name+=',E&S0'
	if 3 in morph:
		labels.append('S0')
		name+=',S0'
	if 4 in morph:
		labels.append('S0a')
		name+=',S0a'
	if 5 in morph:
		labels.append('SpiralA')
		name+=',SA'
	if 6 in morph:
		labels.append('SpiralAB')
		name+=',SAB'
	if 7 in morph:
		labels.append('SpiralB')
		name+=',SB'
	if 8 in morph:
		labels.append('SpiralBC')
		name+=',SBC'
	if 9 in morph:
		labels.append('SpiralC')
		name+=',SC'
	if 10 in morph:
		labels.append('SprialCD')
		name+=',SCD'
	if 11 in morph:
		labels.append('SprialD/Irr')
		name+=',SD&Irr'
		
	params=input('Select parameters: (Redshift=1,Phase=2,Dm15=3,M_B=4,B_mMinusV_m=5)')
	query='SELECT * FROM Supernovae Where '
	if 1 in params:
		range=input('Select redshift range: [xmin,xmax]')
		query += 'redshift BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',redshift;['+str(range[i])+','+str(range[+1])+']'		
	if 2 in params:
		range=input('Select phase range: [xmin,xmax]')
		query += 'Phase BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',phase;['+str(range[0])+','+str(range[1])+']'
	if 3 in params:
		range=input('Select Dm15 range: [xmin,xmax]')
		query += 'Dm15 BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',dm15;['+str(range[0])+','+str(range[1])+']'
	if 4 in params:
		range=input('Select M_B range: [xmin,xmax]')
		query += 'M_B BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',M_B;['+str(range[0])+','+str(range[1])+']'
	if 5 in params:
		range=input('Select B_mMinusV_m range: [xmin,xmax]')
		query += 'B_mMinusV_m BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',B_m-V_m;['+str(range[0])+','+str(range[1])+']'

	query += 'Morphology='

	composites=[]
	queries=[]
	for host in morph:
		composites.append(composite.main(query+str(host)))
		
	queries.append(sys.argv)
	queries.append(str(len(composites)))
	for host in morph:
		queries.append(query+str(host))
		
	print queries

	#labels=['SpiralA,','SpiralAB','SpiralB','SpiralBC','SpiralC','SpiralCD']
	scales=composite.find_scales(composites,composites[0].flux,composites[0].ivar)

	sumphase=0
	sumred=0
	tot=0
	for comp in composites:
		sumphase += comp.phase
		sumred += comp.redshift
		tot += 1
		
	avgphase=sumphase/tot
	avgred=sumred/tot
		
	plot_name = name + ',avgphase-' + str("%.2f" % avgphase) + ',avgred-' + str("%.2f" % avgred)
	names_array.append(plot_name)
	wmin = 0
	wmax = 100000
	for comp in composites:
		SN=comp
		if (SN.minwave > wmin):
			wmin=SN.minwave
		if (SN.maxwave < wmax):
			wmax=SN.maxwave
	
	galrun.main(queries,plot_name,wmin,wmax)
	
cont=0
while(cont==0):
	run()
	if (raw_input('Make another query? (y/n)') == 'y'):
		cont=0
	else:
		cont +=1
	

	"""
	#This makes, shows, and saves a quick comparison plot...we can probably get rid of this when plotting.main works.
	#lowindex = np.where(composites[0].wavelength == composite.find_nearest(composites[0].wavelength, wmin))
	#highindex = np.where(composite[0].wavelength == composite.find_nearest(composites[0].wavelength, wmax))

	plt.figure(1)
	i=0
	for comp in composites:
		plt.plot(comp.wavelength, scales[i]*comp.flux,label=labels[i])
		i+=1
	plt.title(plot_name)
	plt.xlabel('Wavelength')
	plt.ylabel('Flux')
	plt.xlim(wmin,wmax)
	plt.ylim(bottom=0)
	#[lowindex[0]:highindex[0]]
	legend=plt.legend(loc='upper right')
	plt.savefig('../plots/' + plot_name + '.png')
	plt.show()
	"""