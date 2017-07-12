import find_event_data as fed
import matplotlib.pyplot as plt

def plot_light_curves(PH):
	p = PH.light_curves
	times = []
	mags = []
	bands = []
	if p != None:
		for band in p:
			t = []
			m = []
			for i in range(len(p[band][0])):
				t.append(float(p[band][0][i]))
				m.append(float(p[band][1][i][0]))
				# print p[band][1][i][1] #this is the dictionary with additional fields
				# print p[band][1][i][1].get('telescope') #will print telescope field if it exists


			times.append(t)
			mags.append(m)
			bands.append(band)

		for i in range(len(bands)):
			plt.plot(times[i], mags[i], 'o', label = bands[i])
			# if bands[i] == 'I':
			# 	plt.plot(times[i], mags[i], 'o', label = bands[i])
		plt.gca().invert_yaxis()
		plt.legend()
		plt.show()
	else:
		print 'Insufficient Data'

PH_Array = fed.grab_event_phot_data("SELECT * FROM Photometry where dm15_cfa < .7 or dm15_from_fits < .7")
print len(PH_Array)
print PH_Array[1].name
plot_light_curves(PH_Array[1])

# print PH_Array[0].light_curves

