import find_event_data as fed 
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

def data_check(SN_Array, start = None):
	i = start
	dn = np.genfromtxt("data_notes.txt", dtype = 'string')
	file_length = len(dn)
	if file_length > 0:
		end = int(dn[-1][0])
	else:
		end = None

	if start is None:
		with open("data_notes.txt", "a+") as dn_append:
			if file_length == 0:
				i = 0
			else:
				i = end + 1
			print "Starting at index: " + str(i)
			print "Notes Key: \n r = reddened \n g = gap \n s = bad sky subtraction \n gc = galaxy contamination \n cr = bad cosmic ray \n d = discontinuity \n t = telluric absorption"
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
		with open("data_notes.txt", "w+") as dn_redo:
			table = []
			print "Notes Key: \n r = reddened \n g = gap \n s = bad sky subtraction \n gc = galaxy contamination \n cr = bad cosmic ray \n d = discontinuity \n t = telluric absorption"
			print "Type notes separated by commas. Include any comments at the end."
			for line in dn:
				print line[0]
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


SN_Array = fed.find_all_data()
starting_point = raw_input('Where would you like to start? (last line in file): ') or None
if starting_point is not None:
	starting_point = int(starting_point)
data_check(SN_Array, start = starting_point)
