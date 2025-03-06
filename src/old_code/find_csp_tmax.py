import glob
import numpy as np

csp_spectra = glob.glob('../data/spectra/csp/*.dat')
file_path = 'csp_mjd_max.txt'
print file_path

sn_dict = {}
for spec in csp_spectra:
    with open(spec) as s:
        info = [s.next().rstrip().split()[1] for x in range(6)]
        sn_dict[info[0].lower()] = float(info[3]) - 2400000.5

sns = []
tmaxs = []
for sn in sn_dict:
    sns.append(sn)
    tmaxs.append(str(sn_dict[sn]))
with open('../data/info_files/csp_mjd_max.txt', 'w') as file:
  # data = data.T
  np.savetxt(file, np.asarray([sns,tmaxs]).T, fmt = '%s')