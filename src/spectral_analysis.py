import datafidelity as df
import matplotlib.pyplot as plt
import numpy as np 
import copy
from scipy.integrate import simps

def autosmooth(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = df.gsmooth(x_array, y_array, var_y, .001)
        SNR = np.mean(y_array / error)
        # SNR = np.mean(new_y_init / error)
    else:
        var_y = np.ones(len(x_array))
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002)
        error = np.absolute(y_array - new_y_init)
        sm_error = df.gsmooth(x_array, error, var_y, .008)
        SNR = np.median(new_y_init / sm_error)

    if SNR < 5:
        vexp_auto = .0045 #temp value, need to fine tune using SNR in excel spreadsheet
    elif 5 <= SNR < 20:
        vexp_auto = .004
    elif 20 <= SNR < 40:
        vexp_auto = .003
    elif 40 <= SNR < 60:
        vexp_auto = .002
    elif 60 <= SNR < 100:
        vexp_auto = .0015
    else:
        vexp_auto = .001

    return vexp_auto, SNR

def find_vexp(x_array, y_array, var_y=None):
    if var_y is not None:
        error = np.sqrt(var_y)
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002)
        SNR = np.median(new_y_init / error)
    else:
        new_y_init = df.gsmooth(x_array, y_array, var_y, .002) #this smoothing should get in right ballpark
        error = np.absolute(y_array - new_y_init)
        sm_error = df.gsmooth(x_array, error, var_y, .008)
        SNR = np.median(new_y_init / sm_error)

    #TODO: interpolate a function of SNR
    # vexp_line = np.polyfit([2.5, 80], [.0045, .001], 1)
    # coeff_0 = vexp_line[0]
    # coeff_1 = vexp_line[1]
    # results from above:
    coeff_0 = -4.51612903e-05
    coeff_1 = 4.61290323e-03
    vexp_auto = coeff_0*SNR + coeff_1

    if SNR < 2.5:
        vexp_auto = .0045
    if SNR > 80:
        vexp_auto = .001

    return vexp_auto, SNR


def find_extrema(wavelength,sm_flux):
    #input a smoothed spectrum
    dw = 2.
    deriv_flux = np.diff(sm_flux)/dw
    abs_deriv = np.absolute(deriv_flux)
    extrema_locs = np.where(abs_deriv < .0001)
    maxima = []
    for loc in extrema_locs[0]:
        if loc != 0 and loc != len(deriv_flux)-1:
            if deriv_flux[loc-1] > 0 and deriv_flux[loc+1] < 0:
                maxima.append(loc)

    return maxima


def measure_si_ratio(wavelength,flux, varflux = None, vexp=.002, smooth=True, dm15 = None, plot=True):
    # if varflux == None:
    #   varflux = np.zeros(len(wavelength), float)
    if smooth:
        sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
    else:
        sm_flux = flux

    maxima = find_extrema(wavelength,sm_flux)
    max_waves = wavelength[maxima]
    max_fluxes = sm_flux[maxima]


    si_range_1 = np.where((wavelength > 5500.) & (wavelength < 5700.))
    si_range_2 = np.where((wavelength > 5700.) & (wavelength < 6000.))
    si_range_3 = np.where((wavelength > 6000.) & (wavelength < 6500.))

    # si_wave_1 = wavelength[si_range_1]
    # si_wave_2 = wavelength[si_range_2]
    # si_wave_3 = wavelength[si_range_3]

    # si_flux_1 = sm_flux[si_range_1]
    # si_flux_2 = sm_flux[si_range_2]
    # si_flux_3 = sm_flux[si_range_3]

    # si_max_1 = np.amax(si_flux_1)
    # si_max_2 = np.amax(si_flux_2)
    # si_max_3 = np.amax(si_flux_3)

    # si_max_index_1 = np.where(si_flux_1 == si_max_1)
    # si_max_index_2 = np.where(si_flux_2 == si_max_2)  
    # si_max_index_3 = np.where(si_flux_3 == si_max_3)

    # si_max_wave_1 = si_wave_1[si_max_index_1][0]
    # si_max_wave_2 = si_wave_2[si_max_index_2][0]
    # si_max_wave_3 = si_wave_3[si_max_index_3][0]

    m1s = []
    m2s = []
    m3s = []
    for m in maxima:
        if m in si_range_1[0]:
            m1s.append(m)
        if m in si_range_2[0]:
            m2s.append(m)
        if m in si_range_3[0]:
            m3s.append(m)

    if len(m1s) == 0 or len(m2s) == 0 or len(m3s) == 0:
        m1s = []
        m2s = []
        m3s = []
        sm_flux = df.gsmooth(wavelength, flux, varflux, vexp= vexp+.001)
        maxima = find_extrema(wavelength,sm_flux)
        for m in maxima:
            if m in si_range_1[0]:
                m1s.append(m)
            if m in si_range_2[0]:
                m2s.append(m)
            if m in si_range_3[0]:
                m3s.append(m)
    if len(m1s) == 0 or len(m2s) == 0 or len(m3s) == 0:
        print "Could not find maximum in a specified range!"
        return np.nan

    f1s = sm_flux[m1s]
    f2s = sm_flux[m2s]
    f3s = sm_flux[m3s]
    m1 = m1s[-1]
    m2 = m2s[-1]
    m3 = m3s[np.argmax(f3s)]

    # m1 = m1s[0]
    # m2 = m2s[-1]
    # m3 = m3s[-1]
    if dm15 != None and dm15 > 1.7:
        m1 = m1s[0]
        m2 = m2s[np.argmax(f2s)]

    si_max_wave_1 = wavelength[m1]
    si_max_wave_2 = wavelength[m2]
    si_max_wave_3 = wavelength[m3]

    si_max_1 = sm_flux[m1]
    si_max_2 = sm_flux[m2]
    si_max_3 = sm_flux[m3]

    weak_si_trough = np.where((wavelength >= si_max_wave_1) & (wavelength < si_max_wave_2))[0]
    strong_si_trough = np.where((wavelength >= si_max_wave_2) & (wavelength < si_max_wave_3))[0]

    interp_flux = copy.copy(flux)
    interp_flux[weak_si_trough] = np.interp(wavelength[weak_si_trough], [si_max_wave_1, si_max_wave_2], [si_max_1, si_max_2])
    interp_flux[strong_si_trough] = np.interp(wavelength[strong_si_trough], [si_max_wave_2, si_max_wave_3], [si_max_2, si_max_3])


    # v_strong, si_min_wave = measure_velocity(wavelength,flux, 5900., 6300.)
    # v_weak, si_weak_min_wave = measure_weak_si_velocity(wavelength,flux)

    # si_min_index = np.where(wavelength == si_min_wave)
    # si_weak_min_index = np.where(wavelength == si_weak_min_wave)

    #find max of diffs instead of finding minimum
    strong_line_diffs = interp_flux[strong_si_trough] - sm_flux[strong_si_trough]
    weak_line_diffs = interp_flux[weak_si_trough] - sm_flux[weak_si_trough]


    #Line ratio with fractional depths
    # strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])/interp_flux[si_min_index]
    # weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])/interp_flux[si_weak_min_index]

    # strong_line = (interp_flux[si_min_index] - sm_flux[si_min_index])
    # weak_line = (interp_flux[si_weak_min_index] - sm_flux[si_weak_min_index])

    strong_line = np.amax(strong_line_diffs)
    weak_line = np.amax(weak_line_diffs)
    strong_ind = np.where(strong_line_diffs == strong_line)
    weak_ind = np.where(weak_line_diffs == weak_line)

    ratio = weak_line/strong_line
    # ratio = ratio[0]


    # if ratio > .5:
    if plot:
        plt.plot(wavelength,flux)
        plt.plot(wavelength,sm_flux)
        plt.plot(wavelength, interp_flux)
        plt.plot(max_waves, max_fluxes, 'o', color='cyan')
        plt.plot(wavelength[strong_si_trough][strong_ind], interp_flux[strong_si_trough][strong_ind], 'o', color='orange')
        plt.plot(wavelength[strong_si_trough][strong_ind], sm_flux[strong_si_trough][strong_ind], 'o', color='orange')
        plt.plot(wavelength[weak_si_trough][weak_ind], interp_flux[weak_si_trough][weak_ind], 'o', color='orange')
        plt.plot(wavelength[weak_si_trough][weak_ind], sm_flux[weak_si_trough][weak_ind], 'o', color='orange')
        plt.xlim([5000.,7000.])
        plt.ylim([0.,.6])
        plt.show()

    # plt.plot(wavelength,flux)
    # plt.plot(wavelength,sm_flux)
    # plt.plot(wavelength, interp_flux)
    # plt.plot(wavelength[si_min_index], interp_flux[si_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_min_index], sm_flux[si_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_weak_min_index], interp_flux[si_weak_min_index], 'o', color='orange')
    # plt.plot(wavelength[si_weak_min_index], sm_flux[si_weak_min_index], 'o', color='orange')
    # plt.xlim([5000.,7000.])
    # plt.show()

    return ratio

def measure_ca_ratio(wavelength,flux, varflux = None, wave1 = 3550., wave2=3680.):
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)

    # ca_range_1 = np.where((wavelength > 3550.) & (wavelength < 3680.))
    ca_range_1 = np.where((wavelength > wave1) & (wavelength < wave2))
    ca_range_2 = np.where((wavelength > 3890.) & (wavelength < 3970.))

    ca_wave_1 = wavelength[ca_range_1]
    ca_wave_2 = wavelength[ca_range_2]

    ca_flux_1 = sm_flux[ca_range_1]
    ca_flux_2 = sm_flux[ca_range_2]

    ca_max_1 = np.amax(ca_flux_1)
    ca_max_2 = np.amax(ca_flux_2)

    ca_max_index_1 = np.where(sm_flux == ca_max_1)
    ca_max_index_2 = np.where(sm_flux == ca_max_2)

    ca_max_wave_1 = wavelength[ca_max_index_1][0]
    ca_max_wave_2 = wavelength[ca_max_index_2][0]

    print ca_max_2/ca_max_1

    plt.plot(wavelength,flux)
    plt.plot(wavelength,sm_flux)
    plt.plot(wavelength[ca_max_index_1], sm_flux[ca_max_index_1], 'o', color='orange')
    plt.plot(wavelength[ca_max_index_2], sm_flux[ca_max_index_2], 'o', color='orange')
    plt.show()

    return ca_max_2/ca_max_1

def measure_velocity(wavelength, flux, wave1, wave2, vexp=.001, rest_wave=6355., varflux=None, plot=False):

    sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
    si_range = np.where((wavelength > wave1) & (wavelength < wave2))
    si_wave = wavelength[si_range]
    if len(si_wave) == 0:
        return np.nan, np.nan
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)
    if len(si_min_index[0]) > 0.:
        si_min_wave = si_wave[si_min_index][0]

        c = 299792. # km/s
        # rest_wave = 6355. #Angstroms

        v = c*((rest_wave/si_min_wave)**2. - 1)/(1+((rest_wave/si_min_wave)**2.))
        if plot:
            plt.plot(wavelength, flux)
            plt.plot(wavelength, sm_flux)
            plt.plot(si_min_wave, si_min, 'o', color='orange')
            plt.xlim([5000.,7000.])
            plt.ylim([0.,.6])
            plt.show()
    else:
        v = np.nan
        si_min_wave = np.nan
    return v, si_min_wave


#adapted from Rodrigo 
def max_wave(sn, w1, w2, w3):
    
    sm_flux= df.gsmooth(sn.wavelength[sn.x1:sn.x2], sn.flux[sn.x1:sn.x2], None, .001)

    wave_domain_1 = (sn.wavelength[sn.x1:sn.x2] > w1) & (sn.wavelength[sn.x1:sn.x2] < w2)
    elem_flux_1 = np.argmax(sm_flux[wave_domain_1]) #find minimum value within these flux vales to locate "dip
    max_wave_1 = sn.wavelength[sn.x1:sn.x2][wave_domain_1][elem_flux_1] #find the corresponding wavelength
    
    wave_domain_2 = (sn.wavelength[sn.x1:sn.x2] > w2) & (sn.wavelength[sn.x1:sn.x2] < w3)
    elem_flux_2 = np.argmax(sm_flux[wave_domain_2]) #find minimum value within these flux vales to locate "dip
    max_wave_2 = sn.wavelength[sn.x1:sn.x2][wave_domain_2][elem_flux_2] #find the corresponding wavelength

    return max_wave_1, max_wave_2, sm_flux

def measure_boot_EWs(boot_sn_array, w1=7600., w2=8200., w3=9000.):
    EWs = []
    for i, bSN in enumerate(boot_sn_array):
        ew = measure_EW(bSN, w1, w2, w3, plot=True)
        print i, ew
        if not np.isnan(ew):
            EWs.append(ew)
    mean_EW = np.nanmean(EWs)
    var_EW = np.nanstd(EWs)**2

    return mean_EW, var_EW, EWs

def measure_comp_diff_EW(boot_sn_arrays, w1=7600., w2=8200., w3=9000.):
    means = []
    varis = []
    EW_arrs = []
    for arr in boot_sn_arrays:
        mean_ew, var_EW, EWs = measure_boot_EWs(arr, w1=w1, w2=w2, w3=w3)
        means.append(mean_ew)
        varis.append(var_EW)
        EW_arrs.append(EWs)

    diff = (means[0] - means[1])
    err = np.sqrt(varis[0]+varis[1])
    print 'Diff = ', diff, '+/-', err
    plt.hist(EW_arrs[0], color = 'blue')
    plt.hist(EW_arrs[1], color = 'red', alpha=.5)
    plt.show()

    return diff, err, means, varis, EWs

#adapted from Rodrigo 
def measure_EW(sn, w1, w2, w3, plot=False):
    
    max_1, max_2, sm_flux = max_wave(sn, w1, w2, w3)
    
    domain = (sn.wavelength[sn.x1:sn.x2] >= max_1) & (sn.wavelength[sn.x1:sn.x2] <= max_2)
    wave_range = sn.wavelength[sn.x1:sn.x2][domain]
    flux_range = sn.flux[sn.x1:sn.x2][domain]

    line_elem = np.polyfit([max_1, max_2], [sm_flux[np.where(sn.wavelength[sn.x1:sn.x2] == max_1)],
                                            sm_flux[np.where(sn.wavelength[sn.x1:sn.x2] == max_2)]], 1)
    line = line_elem[0] * wave_range + line_elem[1]
    norm = flux_range / line

    a_curve = np.trapz((norm), x = wave_range)
    a_line = max(wave_range) - min(wave_range)
    eq_width = a_line - a_curve
    
    if plot==True:
        plt.figure()
        plt.plot(sn.wavelength[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>w1) & (sn.wavelength[sn.x1:sn.x2]<w3))],
            sn.flux[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>w1) & (sn.wavelength[sn.x1:sn.x2]<w3))])
        plt.plot(sn.wavelength[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>w1) & (sn.wavelength[sn.x1:sn.x2]<w3))],
            sm_flux[np.where((sn.wavelength[sn.x1:sn.x2]>w1) & (sn.wavelength[sn.x1:sn.x2]<w3))], color='g')
#         plt.xlim([6200.,6400.])
        # plt.ylim([0.,3.])
        plt.axvline (x=max_1, color = 'red')
        plt.axvline (x=max_2, color = 'red')
        plt.plot(sn.wavelength[sn.x1:sn.x2][domain],line, color='orange')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.show()
        
#         plt.figure()
#         plt.plot(sn.wavelength[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>6200.) & (sn.wavelength[sn.x1:sn.x2]<6400.))],
#           sn.flux[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>6200.) & (sn.wavelength[sn.x1:sn.x2]<6400.))])
#         plt.plot(sn.wavelength[sn.x1:sn.x2][np.where((sn.wavelength[sn.x1:sn.x2]>6200.) & (sn.wavelength[sn.x1:sn.x2]<6400.))],
#           sm_flux[np.where((sn.wavelength[sn.x1:sn.x2]>6200.) & (sn.wavelength[sn.x1:sn.x2]<6400.))], color='blue')
# #         plt.xlim([6200.,6400.])
# #         plt.ylim([20.,30.])
#         plt.axvline (x=max_1, color = 'red')
#         plt.axvline (x=max_2, color = 'red')
#         plt.plot(sn.wavelength[sn.x1:sn.x2][domain],line, color='orange')
#         plt.xlabel('Wavelength')
#         plt.ylabel('Flux')
#         plt.show()

#         plt.plot(norm, color = 'black')
#         plt.plot(wave_range, norm, color = "orange")
#         plt.plot(wave_range , 1, color = "yellow")
#         plt.axvline (x=max_1, color = 'red')
#         plt.axvline (x=max_2, color = 'blue')
    return eq_width

def measure_weak_si_velocity(wavelength, flux, varflux=None, plot=True):
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    varflux = np.zeros(len(wavelength), float)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    si_range = np.where((wavelength > 5560.) & (wavelength < 5820.))
    si_wave = wavelength[si_range]
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)
    si_min_wave = si_wave[si_min_index][0]

    c = 299792. # km/s
    si_rest_wave = 5979. #Angstroms

    v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
    # if plot:
    #   plt.plot(wavelength, flux)
    #   plt.plot(wavelength, sm_flux)
    #   plt.plot(si_min_wave, si_min, 'o')
    #   plt.xlim([5000.,7000.])
    #   plt.ylim([0.,.6])
    #   plt.show()
    return v, si_min_wave

    #carbon II is 6580A (2006bt good example)

def measure_si_velocity_from_raw(wavelength,flux,z, varflux=None):
    # varflux = 1./SN.ivar[SN.x1:SN.x2]
    # plt.plot(wavelength,flux)
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)

    wavelength = wavelength/(1.+z)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    si_range = np.where((wavelength > 5830.) & (wavelength < 6230.))
    si_wave = wavelength[si_range]
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)
    si_min_wave = si_wave[si_min_index][0]

    c = 299792. # km/s
    si_rest_wave = 6355. #Angstroms
    v = c*((si_rest_wave/si_min_wave)**2. - 1)/(1+((si_rest_wave/si_min_wave)**2.))
    # v = c*(si_min_wave - si_rest_wave)/si_rest_wave
    # plt.plot(wavelength, flux)
    # plt.plot(wavelength, sm_flux)
    # plt.plot(si_min_wave, si_min, 'o', color='orange')
    # plt.show()
    return v, si_min_wave

def measure_C_velocity_from_raw(wavelength,flux,z, varflux=None):
    # varflux = 1./SN.ivar[SN.x1:SN.x2]
    # plt.plot(wavelength,flux)
    if varflux == None:
        varflux = np.zeros(len(wavelength), float)
    wavelength = wavelength/(1.+z)
    sm_flux = df.gsmooth(wavelength, flux, varflux, .001)
    C_range = np.where((wavelength > 6200.) & (wavelength < 6320.))
    C_wave = wavelength[C_range]
    C_flux = sm_flux[C_range]
    C_min = np.amin(C_flux)
    C_min_index = np.where(C_flux == C_min)
    C_min_wave = C_wave[C_min_index][0]

    c = 299792. # km/s
    C_rest_wave = 6580. #Angstroms
    v = c*((C_rest_wave/C_min_wave)**2. - 1)/(1+((C_rest_wave/C_min_wave)**2.))
    # v = c*(si_min_wave - si_rest_wave)/si_rest_wave
    # plt.plot(wavelength, flux)
    # plt.plot(wavelength, sm_flux)
    # plt.plot(C_min_wave, C_min, 'o', color='orange')
    # plt.show()
    return v, C_min_wave

def resample_spectrum(SN):
    wavelength = SN.wavelength[SN.x1:SN.x2]
    flux = SN.flux[SN.x1:SN.x2]
    varflux = 1./SN.ivar[SN.x1:SN.x2]


