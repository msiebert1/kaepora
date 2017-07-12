import matplotlib.pyplot as plt
import numpy as np
import query_db
import composite
from tabulate import tabulate

SN_Array = composite.grab("SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between -.2 and .2", multi_epoch = True)
tab_arr = []
for SN in SN_Array:
	wav_range = str(SN.minwave) + ' - ' + str(SN.maxwave)
	ref = '...'
	tab_arr.append([SN.name, SN.source, SN.mjd, SN.phase, wav_range, ref])
table = tabulate(tab_arr, headers=['SN Name', 'Source', 'mjd', 'Phase', 'Wavelength Range', 'Reference'], tablefmt = 'latex')

with open('../../Paper_Drafts/table_1.tex', 'w') as file:
	file.write('\\begin{table*}[t]\n')
	file.write('\centering\n')
	file.write(table)
	file.write('\n')
	file.write('\caption{CAPTION}\n')
	file.write('\label{tab:1}\n')
	file.write('\end{table*}')
	file.close()
	# print '\\begin{table*}[t]'
	# print '\centering'
	# print table
	# print '\caption{Blabla}'
	# print '\label{tab:1}'
	# print '\end{table*}'