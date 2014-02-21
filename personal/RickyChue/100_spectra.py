import numpy as np
import matplotlib.pyplot as plt


x = np.arange(0, 100, 1)

for i in range(1000):
    outfile = open('random/' + str(i) + '.dat', 'w')
    y = np.random.rand(100)
    for m in range(len(x)):
        outfile.write(str(x[m]) + '\t' + str(y[m]) + '\n')
    outfile.close()
    