import numpy as np
import matplotlib.pyplot as plt

plt.yscale('log')

f1 = open('../../data/sn2011by-hst+lick.flm')
lines = f1.readlines()
f1.close()

x1 = []
y1 = []
for line in lines:
	p = line.split()
	x1.append(float(p[0]))
	y1.append(float(p[1]))
	
x1=np.array(x1)
y1=np.array(y1)
x1=x1[0:2689]
y1=y1[0:2689]
p1,=plt.plot(x1, y1)


f2=open('../../data/sn2011fe-visit3-hst.flm')
lines=f2.readlines()
f2.close()

x2=[]
y2=[]
for line in lines:
	p=line.split()
	x2.append(float(p[0]))
	y2.append(float(p[1]))
	
x2=np.array(x2)
y2=np.array(y2)
p2,=plt.plot(x2, y2)

avx = x1
avy = (y1+y2)/2
p3,=plt.plot(avx,avy)

plt.xlabel('Wavelength [A]')
plt.ylabel('Flux')
plt.legend([p1,p2,p3],['SN2011BY', 'SN2011FE','Average'],2)
plt.savefig('avgspectrum.png')
plt.show()
