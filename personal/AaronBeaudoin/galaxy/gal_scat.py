from astropy.table import Table
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import math
from pylab import *
import matplotlib.gridspec as gridspec
import scipy.optimize as optimize
import numpy as np
import operator
import sys
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import msgpack as msg
import msgpack_numpy as mn

class supernova(object):
	""" """
	
def grab(query):
	con = sq3.connect('../../../data/SNe.db')
	cur = con.cursor()
	
	SN_Array=[]
	cur.execute(query)
	for row in cur:
		SN = supernova()
		SN.name=row[1]
		SN.phase=row[4]
		SN.dm15=row[7]
		SN.Mb=row[8]
		SN.Bminusv=row[9]
		SN.morph=row[11]
		SN_Array.append(SN)

	return SN_Array

def get_query():
	name=''
	morph_range=[]
	morph_range.append([1,4])
	morph_range.append([5,10])
	morph_range.append([11,11])
	
	params=input('Select parameters: (Redshift=1,Phase=2,Dm15=3,M_B=4,B_mMinusV_m=5)')
	query='SELECT * FROM Supernovae Where '
	if 1 in params:
		range=input('Select redshift range: [xmin,xmax]')
		query += 'redshift BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',redshift['+str(range[0])+','+str(range[1])+']'		
	if 2 in params:
		range=input('Select phase range: [xmin,xmax]')
		query += 'Phase BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',phase['+str(range[0])+','+str(range[1])+']'
	if 3 in params:
		range=input('Select Dm15 range: [xmin,xmax]')
		query += 'Dm15 BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',dm15['+str(range[0])+','+str(range[1])+']'
	if 4 in params:
		range=input('Select M_B range: [xmin,xmax]')
		query += 'M_B BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',M_B['+str(range[0])+','+str(range[1])+']'
	if 5 in params:
		range=input('Select B_mMinusV_m range: [xmin,xmax]')
		query += 'B_mMinusV_m BETWEEN ' + str(range[0]) + ' AND ' + str(range[1]) + ' AND '
		name+=',B_m-V_m['+str(range[0])+','+str(range[1])+']'

	query += 'Morphology '

	queries=[]
	for morph in morph_range:
		queries.append(query + ' BETWEEN ' + str(morph[0]) + ' AND ' + str(morph[1]))
		
	return queries, name

def main():
	queries,plot_name=get_query()
	print queries
	
	SN_E=grab(queries[0])
	SN_S=grab(queries[1])
	SN_I=grab(queries[2])
	plot_name+=',E'+str(len(SN_E))+',S'+str(len(SN_S))+',I'+str(len(SN_I))+'.png'	
	title='Ellipticals: ' + str(len(SN_E)) + ', Spirals: ' + str(len(SN_S)) + ', Irregulars: ' + str(len(SN_I))
	
	font = {'family' : 'serif','color'  : 'black','weight' : 'bold','size' : 10,} 
	
	fig=plt.figure(num=1,dpi=100,figsize=[8,10],facecolor='w')
	gs=gridspec.GridSpec(2,1)
	ax1=fig.add_subplot(gs[0])
	for SN in SN_E:
		plt.scatter(SN.dm15,SN.Bminusv,c='b',marker='^')
	for SN in SN_S:
		plt.scatter(SN.dm15,SN.Bminusv,c='g',marker='s')
	for SN in SN_I:
		plt.scatter(SN.dm15,SN.Bminusv,c='r',marker='o')
	plt.scatter(SN_E[0].dm15,SN_E[0].Bminusv,c='b',marker='^',label='Elliptical')
	plt.scatter(SN_S[0].dm15,SN_S[0].Bminusv,c='g',marker='s',label='Spiral')
	plt.scatter(SN_I[0].dm15,SN_I[0].Bminusv,c='r',marker='o',label='Irregular')
	plt.gca().invert_yaxis()
	l=plt.legend(prop={'family':'serif'})
	l.draw_frame(False)
	plt.draw()
	plt.xlabel('$\Delta$m$_{15}$(B)', fontdict = font)
	plt.ylabel('(B-v)$_{GAL}$',fontdict=font)
	ax2=fig.add_subplot(gs[1])
	for SN in SN_E:
		plt.scatter(SN.phase,SN.morph,c='b',marker='^')
	for SN in SN_S:
		plt.scatter(SN.phase,SN.morph,c='g',marker='s')
	for SN in SN_I:
		plt.scatter(SN.phase,SN.morph,c='r',marker='o')
	plt.xlabel('Age',fontdict=font)
	plt.setp(ax2.get_yticklabels(), visible=False)
	plt.savefig('dm15_vs_bminusv' + plot_name,dpi=600, facecolor='w',edgecolor='w',pad_inches=.1)
	plt.show()
	
	fig=plt.figure(num=2,dpi=100,figsize=[8,10],facecolor='w')
	gs=gridspec.GridSpec(2,1,hspace=.001)
	ax1=fig.add_subplot(gs[0])
	for SN in SN_E:
		plt.scatter(SN.Mb,SN.dm15,c='b',marker='^')
	for SN in SN_S:
		plt.scatter(SN.Mb,SN.dm15,c='g',marker='s')
	for SN in SN_I:
		plt.scatter(SN.Mb,SN.dm15,c='r',marker='o')
	plt.scatter(SN_E[0].Mb,SN_E[0].dm15,c='b',marker='^',label='Elliptical')
	plt.scatter(SN_S[0].Mb,SN_S[0].dm15,c='g',marker='s',label='Spiral')
	plt.scatter(SN_I[0].Mb,SN_I[0].dm15,c='r',marker='o',label='Irregular')
	l=plt.legend(prop={'family':'serif'})
	l.draw_frame(False)
	plt.draw()
	plt.ylabel('$\Delta$m$_{15}$(B)', fontdict = font)
	plt.gca().invert_xaxis()
	ax2=fig.add_subplot(gs[1])
	for SN in SN_E:
		plt.scatter(SN.Mb,SN.Bminusv,c='b',marker='^')
	for SN in SN_S:
		plt.scatter(SN.Mb,SN.Bminusv,c='g',marker='s')
	for SN in SN_I:
		plt.scatter(SN.Mb,SN.Bminusv,c='r',marker='o')
	plt.ylabel('(B-v)$_{GAL}$',fontdict=font)
	plt.xlabel('M$_{B}$',fontdict=font)
	plt.gca().invert_xaxis()
	plt.gca().invert_yaxis()
	plt.setp(ax1.get_xticklabels(), visible=False)
	plt.savefig('Mb_vs_dm15&bminusv'+plot_name,dpi=600,facecolor='w',edgecolor='w',pad_inches=.1)
	plt.show()
	
	
cont=0
while(cont==0):
	main()
	if (raw_input('Make another query? (y/n)') == 'y'):
		cont=0
	else:
		cont +=1	