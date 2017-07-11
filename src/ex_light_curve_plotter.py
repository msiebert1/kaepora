import find_event_data as fed
import matplotlib.pyplot as plt

def plot_light_curves(PH):
	p = PH.light_curves
	times = []
	mags = []
	bands = []
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

PH_Array = fed.grab_event_phot_data("SELECT * FROM Photometry where SN = '2001bf'")
plot_light_curves(PH_Array[0])

# print PH_Array[0].light_curves

