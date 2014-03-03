import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy 
from math import *

def smoothing(X_array, Y_array, var_Y, vexp = 0.01, nsig = 5.0):

    for i in range(len(var_Y)):
        if var_Y[i] == 0:
            var_Y[i] = 1E-20
    new_Y = np.zeros(len(X_array), float)

    for i in range(len(X_array)):
        gaussian = np.zeros(len(X_array), float)
        sigma = vexp*X_array[i]

        sigrange = np.nonzero(abs(X_array-X_array[i]) <= nsig*sigma)

        gaussian[sigrange] = (1/(sigma*sqrt(2*pi)))*np.exp(-0.5*((X_array[sigrange]-X_array[i])/sigma)**2)

        W_lambda = gaussian / var_Y
        W0 = np.sum(W_lambda)
        W1 = np.sum(W_lambda*Y_array)

        new_Y[i] = W1/W0

    return new_Y        
