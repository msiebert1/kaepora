import find_event_data as fed
import matplotlib.pyplot as plt

def plot_light_curves(PH):
	p = PH.light_curves
	times = []
	ms = []
	for band in p:
		t = []
		m = []
		for i in range(len(p[band][0])):
			t.append(float(p[band][0][i]))
			m.append(float(p[band][1][i][0]))
			print p[band][1][i][1] #this is the dictionary with additional fields
			print p[band][1][i][1].get('telescope') #will print telescope field if it exists


		times.append(t)
		ms.append(m)
	for i in range(len(times)):
		plt.plot(times[i], ms[i], 'o')
	plt.gca().invert_yaxis()
	plt.show()

# SN_Array, PH_Array = fed.find_data('2001bf')
# plot_light_curves(PH_Array[0])
# print PH_Array[0].light_curves

PH_Array = fed.grab_event_phot_data("SELECT * FROM Photometry where SN = '2001bf'")
print len(PH_Array)
# SN_Array = fed.grab_event_data("SELECT * FROM Supernovae where dm15")

# names = []
# dm15s = []
# for SN in SN_Array:
# 	names.append(SN.name)
# 	dm15s.append(SN.dm15)

param1 = []
param2 = []
# same_event_dm15s
for PH in PH_Array:
	# if PH.name is in names:
	# 	same_event_dm15s.append(dm15s[names.index('')])
	param1.append(PH.x1_salt2)
	param2.append(PH.delta_mlcs17)

plt.scatter(param1, param2)
plt.show()
