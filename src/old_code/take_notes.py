import find_event_data as fed 
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

def data_check(SN_Array, source, start = None):
	i = start

	if source == 'all':
		file = "data_notes.txt"
		dn = np.genfromtxt(file, dtype = 'string')
		file_length = len(dn)
		if file_length > 0:
			end = int(dn[-1][0])
		else:
			end = None
	elif source == 'csp/other':
		file = "data_notes_csp_other.txt"
		dn = np.genfromtxt(file, dtype = 'string')
		file_length = len(dn)
		if file_length > 0:
			end = int(dn[-1][0])
		else:
			end = None

	if start is None:
		with open(file, "a+") as dn_append:
			if file_length == 0:
				i = 0
			else:
				i = end + 1
			print "Starting at index: " + str(i)
			print "Notes Key: \n r = reddened \n g = gap \n s = bad sky subtraction \n gc = galaxy contamination \n cr = bad cosmic ray \n d = discontinuity \n t = telluric absorption \n q = questionable"
			print "Type notes separated by commas. Include any comments at the end."
			table = []
			while i < len(SN_Array):
				print SN_Array[i].filename
				plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
				plt.show()
				notes = raw_input("Notes: ")
				if notes == '':
					notes = 'None'
				table.append([str(i), SN_Array[i].filename, notes])
				# dn_append.write(tabulate(table))
				# dn_append.write(str(i) + "\t\t" + SN_Array[i].filename + "\t\t" + "Notes: " + notes + "\n")
				cont = raw_input("Continue?(q to quit): ")
				if cont == 'q':
					break
				i += 1
			dn_append.write('\n')
			dn_append.write(tabulate(table, tablefmt = 'plain'))

	else:	
		dn_new = []
		with open(file, "w+") as dn_redo:
			table = []
			print "Notes Key: \n r = reddened \n g = gap \n s = bad sky subtraction \n gc = galaxy contamination \n cr = bad cosmic ray \n d = discontinuity \n t = telluric absorption \n q = questionable"
			print "Type notes separated by commas. Include any comments at the end."
			for line in dn:
				if float(line[0]) == i:
					print line[0], line[1], line[2]
					print SN_Array[i].filename
					plt.plot(SN_Array[i].wavelength, SN_Array[i].flux)
					plt.show()
					notes = raw_input("Notes: ")
					if notes == '':
						notes = 'None'
					table.append([str(i), SN_Array[i].filename, notes])
					# dn_redo.write(str(i) + "\t" + SN_Array[i].filename + "\t" + notes + "\n")
					cont = raw_input("Continue?(q to quit): ")
					if cont != 'q':
						i += 1
					else:
						i = None
				else:
					table.append([line[0], line [1], line[2]])
					# dn_redo.write('\t'.join(line) + "\n")
			dn_redo.write(tabulate(table, tablefmt = 'plain'))

	return 


source = raw_input('Choose data source: \n 1 = all data \n 2 = csp/other \n:')

if source == '1':
	source == 'all'
	SN_Array = fed.find_all_data()
elif source == '2':
	source = 'csp/other'
	SN_Array = fed.find_csp_other_data()

starting_point = raw_input('Where would you like to start? (last line in file): ') or None
if starting_point is not None:
	starting_point = int(starting_point)
data_check(SN_Array, source, start = starting_point)
