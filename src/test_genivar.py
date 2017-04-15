import datafidelity as df 
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import os
import prep as prep

def scale_composites_in_range(data, comp):
    scales = []
    guess = 2.
    s = opt.minimize(sq_residuals_in_range, guess, args = (data, comp), 
                 method = 'Nelder-Mead').x
    return s

def sq_residuals_in_range(s, data, comp):
    data = s*data
    res = data - comp
    sq_res = res*res
    return np.sum(sq_res)

#test more spectra
scales = []
i = 0

fname = '../data/spectra/cfa/sn2007f/sn2007F-20070114.49-fast.flm'
spectrum = np.loadtxt(fname)

old_wave = spectrum[:, 0]
old_flux = spectrum[:, 1]

# real_error = spectrum[:, 2]
# if real_error[0] != 0.:
old_error = np.zeros(len(old_wave), float)
new_ivar = df.genivar(old_wave, old_flux, old_error)
output = prep.Interpo(old_wave, old_flux, new_ivar)

new_var = 1./new_ivar
plt.subplot(2,1,1)
plt.plot(output[0], output[1])
plt.subplot(2,1,2)
plt.plot(output[0], output[2])
plt.plot()
plt.show()


# for path, subdirs, files in os.walk('../data/spectra/cfa'):
# 	if i < 1000:
# 		for fname in files:

# 			fname = os.path.join(path, fname)
# 			if fname.endswith('.flm'):
# 				spectrum = np.loadtxt(fname)

# 				old_wave = spectrum[:, 0]
# 				old_flux = spectrum[:, 1]

# 				real_error = spectrum[:, 2]
# 				if real_error[0] != 0.:
# 					old_error = np.zeros(len(old_wave), float)
# 					new_ivar = df.genivar(old_wave, old_flux, old_error)

# 					new_var = 1./new_ivar
# 					real_var = real_error*real_error
# 					scale = scale_composites_in_range(new_var, real_var)
# 					new_var = new_var*scale
# 					print scale
# 					scales.append(scale)
# 					i += 1
# 	else:
# 		break

			# plt.subplot(2,1,1)
			# plt.plot(old_wave, old_flux)
			# plt.subplot(2,1,2)
			# plt.plot(old_wave, real_var)
			# plt.plot(old_wave, new_var)
			# plt.show()


# print np.average(scales)
# print np.std(scales)

#avg: 2.02 std: 1.055