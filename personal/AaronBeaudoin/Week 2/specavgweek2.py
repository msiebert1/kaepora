from astropy.table import Table
import numpy as np

data1=Table.read('sn2011by-hst+lick.flm',format='ascii')
data2=Table.read('sn2011fe-visit3-hst.flm',format='ascii')
x1=data1["col1"]
y1=data1["col2"]
x2=data2["col1"]
y2=data2["col2"]

yy2=np.linspace(0,len(x1)-1,len(x1))

lines1=np.where((x1<min(x2)))
lines2=np.where((x1>min(x2)-1) & (x1<max(x2)))
lines3=np.where(x1>max(x2)-1)
yy2[lines1]=0
yy2[lines2]=y2
yy2[lines3]=0

ya=(y1+yy2)/2


t=Table([x1,ya],names=('col1','col2'))
t.write('snavg.flm',format='ascii')
