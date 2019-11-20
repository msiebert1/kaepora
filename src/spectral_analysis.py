import datafidelity as df
import matplotlib.pyplot as plt
import numpy as np 
import copy
from scipy.integrate import simps
import random 
import pyphot
import kaepora as kpora

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

def measure_verror(wavelength, flux, var_flux, n=100):
    vdist = []
    for i in range(n):
        sample_vexp = np.random.uniform(.001, .0045)
        sample_v, sample_si_min_wave = measure_velocity(wavelength, flux, 5800, 6300, vexp=sample_vexp, plot=False, error=False)
        vdist.append(sample_v)

    sigma = np.std(vdist)
    return sigma

def measure_velocity(wavelength, flux, wave1, wave2, vexp=.001, rest_wave=6355., varflux=None, plot=False, error=False):

    sm_flux = df.gsmooth(wavelength, flux, varflux, vexp)
    si_range = np.where((wavelength > wave1) & (wavelength < wave2))
    si_wave = wavelength[si_range]
    if len(si_wave) == 0:
        return np.nan, np.nan
    si_flux = sm_flux[si_range]
    si_min = np.amin(si_flux)
    si_min_index = np.where(si_flux == si_min)

    if len(si_min_index[0]) > 0. and (wavelength[-1] > wave2):
        si_min_wave = si_wave[si_min_index][0]

        c = 299792.458 # km/s
        # rest_wave = 6355. #Angstroms

        v = c*((rest_wave/si_min_wave)**2. - 1)/(1+((rest_wave/si_min_wave)**2.))
        if error:
            sigma = measure_verror(wavelength, flux, varflux)

        if plot:
            norm = 1./np.amax(flux)
            plt.plot(wavelength, norm*flux)
            plt.plot(wavelength, norm*sm_flux)
            plt.plot(si_min_wave, norm*si_min, 'o', color='orange')
            plt.xlim([5000.,7000.])
            plt.ylim([0.,.6])
            plt.show()
    else:
        v = np.nan
        si_min_wave = np.nan

    if error:
        return (-1.*v)/1000., si_min_wave, sigma
    else:
        return (-1.*v)/1000., si_min_wave

def measure_vels(comps, sn_arrs, attr_name, boot_arrs = None, plot=False):
    avg_vs = []
    for comp in comps:
        v, si_min_wave = measure_velocity(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], 5800, 6400, plot=plot)
        avg_vs.append(v)
        print v

    if boot_arrs:
        errors = []
        boot_vel_arrs = []
        
        for boots in boot_arrs:
            b_vs = []
            for b in boots:
                v, si_min_wave = measure_velocity(b.wavelength[b.x1:b.x2], b.flux[b.x1:b.x2], 5800, 6400, plot=False)
                # print v
                if ~np.isnan(v) and v != None:
                    b_vs.append(v)
            boot_vel_arrs.append(b_vs)
            # errors.append(np.std(b_vs))
            plt.hist(b_vs)
            plt.show()
            
        low_errors = []
        up_errors = []
        for vlist in boot_vel_arrs:
            p = np.percentile(vlist, [18, 50, 82])
            vlower = p[1] - p[0]
            vupper = p[2] - p[1]
            low_errors.append(vupper)
            up_errors.append(vlower)
            print 'ERR: ', (vupper+vlower)/2

        errors = [low_errors, up_errors]
        
    avg_attrs = []
    for arr in sn_arrs:
        attr_list = []
        for i, spec in enumerate(arr):
            attr_list.append(spec.other_meta_data[attr_name])
        avg_attrs.append(np.average(attr_list))
    
    #Caveat: this includes combined spectra
    vels = []
    attrs = []
    attr_errs = []
    for arr in sn_arrs:
        for i, spec in enumerate(arr):
            v, si_min_wave = measure_velocity(spec.wavelength[spec.x1:spec.x2], spec.flux[spec.x1:spec.x2], 5800, 6300, plot=plot)
            attrs.append(spec.other_meta_data[attr_name])
            attr_errs.append(0)
            vels.append(v)
            print i,spec.name, v
    
    # if boot_arrs:
    #     plt.errorbar(avg_attrs, avg_vs, yerr=errors, fmt='o')
    # else:
    #     plt.scatter(avg_attrs, avg_vs)
    # plt.scatter(attrs, vels, color='orange')
    # plt.show()

    return boot_vel_arrs, avg_vs, avg_attrs, attrs, attr_errs, vels, errors

def plot_vels(vel_data, vwidth = .1, savename=None):
    boot_vel_arrs, avg_vs, avg_attrs, attrs, attr_errs, vels, errors = vel_data[0], vel_data[1], vel_data[2], vel_data[3], vel_data[4], vel_data[5], vel_data[6]
    plt.rc('font', family='serif')
    plt.figure(num = 1, dpi = 100, figsize = [7,7])
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    #     ax.get_yaxis().set_ticks([])
    #     plt.ylim([-311,12])
    plt.ylabel('Velocity ($10^3$ km s$^{-1}$)', fontsize = 20)
    # plt.xlabel('Hubble Residual (mag)', fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=20)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=10)

    plt.errorbar(attrs, vels, fmt='o', markersize=5, color = 'black', capsize=5, zorder=-10, label='Individual SNe')
    violin_parts = plt.violinplot(boot_vel_arrs, avg_attrs, points=40, widths=vwidth,
                      showmeans=False, showextrema=False, showmedians=False,
                      bw_method='silverman')
    
    print avg_vs
    diff = avg_vs[1] - avg_vs[0]
    # err = np.sqrt(errors[0][0]**2. + errors[1][1]**2.)
    # print 'Vdiff: ', diff, 'Err:', err, 'Sig: ', diff/err
    plt.errorbar(avg_attrs, avg_vs, fmt='*', markersize=20, color = 'orange', alpha=1., zorder=1, label= 'Composite \n Spectra')
#         vp = violin_parts['cmeans']
#         vp.set_edgecolor('darkblue')
#         vp.set_linewidth(4)
#         vp.set_alpha(1.)
    for pc in violin_parts['bodies']:
        pc.set_facecolor('deepskyblue')
        pc.set_edgecolor('steelblue')
#             pc.set_color('red')
        pc.set_alpha(.5)
    # plt.xlim(-20,-7)
    # plt.axvline(x=0., linestyle = '--', color = 'deepskyblue', linewidth=3, alpha=.6)
    # plt.text(-.45, -14.5, '+4 Days', fontsize=20)
    # plt.legend(bbox_to_anchor=(0.48, 0.45, 0.48, 0.5), fontsize=15)
    plt.gca().invert_yaxis()
    if savename is not None:
        plt.savefig(savename+'.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()
#adapted from Rodrigo 
def max_wave(wavelength, flux, w1, w2, w3, vexp=.001, sys_error=False):
    
    sm_flux= df.gsmooth(wavelength, flux, None, vexp)

    wave_domain_1 = (wavelength > w1) & (wavelength < w2)
    elem_flux_1 = np.argmax(sm_flux[wave_domain_1]) #find minimum value within these flux vales to locate "dip
    max_wave_1 = wavelength[wave_domain_1][elem_flux_1] #find the corresponding wavelength
    
    wave_domain_2 = (wavelength > w2) & (wavelength < w3)
    elem_flux_2 = np.argmax(sm_flux[wave_domain_2]) #find minimum value within these flux vales to locate "dip
    max_wave_2 = wavelength[wave_domain_2][elem_flux_2] #find the corresponding wavelength

    if sys_error:
        #change continuum locations by a random even integer (at most 100 Angstroms)
        max_wave_1 = max_wave_1 + random.randint(-50, 50)*2.
        max_wave_2 = max_wave_2 + random.randint(-50, 50)*2.

    return max_wave_1, max_wave_2, wave_domain_1, wave_domain_2, elem_flux_1, elem_flux_2, sm_flux

def measure_EWs(sn_array, w1=7600., w2=8200., w3=9000., error=False):
    EWs = []
    phases = []
    errs = []
    morphs = []
    SNs = []
    for i, SN in enumerate(sn_array):
        if SN.wavelength[SN.x2] >= w3:
            vexp_auto, SNR = autosmooth(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], var_y=None)
            ew, domain = measure_EW(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], w1, w2, w3, vexp=vexp_auto, plot=True)
            if error:
                stat_err = ew_stat_error(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2], w1, w2, w3, domain, vexp=vexp_auto, num=50)
                sys_err = ew_sys_error(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], SN.ivar[SN.x1:SN.x2], w1, w2, w3, domain, vexp=vexp_auto, num=50)
            else:
                stat_err = None
                sys_err = None
            print i, ew, sys_err, stat_err, SN.phase
            err_tot = np.sqrt(sys_err**2. + stat_err**2.)
            if not np.isnan(ew) and ew < 500. and stat_err < 100.:
                EWs.append(ew)
                phases.append(SN.phase)
                errs.append(err_tot)
                morphs.append(SN.ned_host)
                SN.nir_EW = ew
                SN.nir_EW_err = err_tot
                SNs.append(SN)
    mean_EW = np.nanmean(EWs)
    var_EW = np.nanstd(EWs)**2

    return mean_EW, var_EW, EWs, errs, phases, morphs, SNs

def measure_comp_diff_EW(boot_sn_arrays, w1=7600., w2=8200., w3=9000.):
    means = []
    varis = []
    EW_arrs = []
    for arr in boot_sn_arrays:
        mean_ew, var_EW, EWs = measure_EWs(arr, w1=w1, w2=w2, w3=w3)
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
def measure_EW(wavelength, flux, w1, w2, w3, vexp=.001, plot=False, sys_error=False):
    
    max_wave_1, max_wave_2, wave_domain_1, wave_domain_2, elem_flux_1, elem_flux_2, sm_flux = max_wave(wavelength, flux, w1, w2, w3, vexp=vexp, sys_error=sys_error)
    

    domain = (wavelength >= max_wave_1) & (wavelength <= max_wave_2)
    roi = (wavelength > w1) & (wavelength < w3)
    wave_range = wavelength[domain]
    flux_range = flux[domain]


    # line_elem = np.polyfit([max_wave_1, max_wave_2], [sm_flux[np.where(wavelength == max_wave_1)],
    #                                         sm_flux[np.where(wavelength == max_wave_2)]], 1)
    line_elem = np.polyfit([max_wave_1, max_wave_2], [sm_flux[wave_domain_1][elem_flux_1],
                                            sm_flux[wave_domain_2][elem_flux_2]], 1)

    line = line_elem[0] * wave_range + line_elem[1]
    norm = flux_range / line

    a_curve = np.trapz((norm), x = wave_range)
    a_line = max(wave_range) - min(wave_range)
    eq_width = a_line - a_curve
    
    if plot==True:
        plt.figure()
        plt.plot(wavelength[np.where((wavelength>w1) & (wavelength<w3))],
            flux[np.where((wavelength>w1) & (wavelength<w3))])
        plt.plot(wavelength[np.where((wavelength>w1) & (wavelength<w3))],
            sm_flux[np.where((wavelength>w1) & (wavelength<w3))], color='g')
#         plt.xlim([6200.,6400.])
        # plt.ylim([0.,3.])
        plt.axvline (x=max_wave_1, color = 'red')
        plt.axvline (x=max_wave_2, color = 'red')
        plt.plot(wavelength[domain],line, color='orange')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.show()
    
    return eq_width, roi

def ew_stat_error(wavelength, flux, ivar, w1, w2, w3, roi, vexp=.001, num=100):
    sm_flux = df.gsmooth(wavelength[roi], flux[roi], 1/ivar[roi], vexp)
    err = np.sqrt((1./ivar[roi]))
    # err = np.absolute(flux[roi] - sm_flux)
    sig = np.median(err)
    ews = []
    for i in range(0, num):
        new_flux = copy.deepcopy(sm_flux)
        delta_flux = np.random.normal(loc=0, scale=sig, size=len(new_flux))
        new_flux = new_flux + delta_flux
        ew, roi_ignore = measure_EW(wavelength[roi], new_flux, w1, w2, w3, vexp=vexp, plot=False)
        ews.append(ew)
    stat_err = np.std(ews)

    return stat_err

def ew_sys_error(wavelength, flux, ivar, w1, w2, w3, roi, vexp=.001, num=100):
    ews = []
    for i in range(0, num):
        # vexp = random.uniform(0.001, 0.0045)
        ew, roi_ignore = measure_EW(wavelength[roi], flux[roi], w1, w2, w3, vexp=vexp, plot=False, sys_error=True)
        ews.append(ew)
    sys_err = np.std(ews)

    return sys_err


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

    c = 299792.458 # km/s
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


def measure_mags(wave, flux, filts = ['GROUND_JOHNSON_B','GROUND_JOHNSON_V']):
    lib = pyphot.get_library()
    mags = []
    for f in filts:
        func = lib[f]
        fflux = func.get_flux(wave, flux, axis = -1)
        mag = -2.5 * np.log10(fflux) - func.Vega_zero_mag
        mags.append(mag)
    return mags

def measure_comp_1m2(comps, filts = ['GROUND_JOHNSON_B','GROUND_JOHNSON_V'], boot_arrs = None, error=False):
    lib = pyphot.get_library()
    func1 = lib[filts[0]]
    func2 = lib[filts[1]]
    
    comp_1 = []
    comp_2 = []
    phases = []
    errors = []
    kpora.set_min_num_spec(comps, 5)
    for comp in comps:
        mag1, mag2 = measure_mags(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], filts = filts)
        comp_1.append(mag1)
        comp_2.append(mag2)
        phases.append(np.average(comp.phase_array[comp.x1:comp.x2]))
        
    if error:
        boot_1m2s = []
        errors = []
        for boots in boot_arrs:
            boot_1 = []
            boot_2 = []
            for b in boots:
                mag1, mag2 = measure_mags(b.wavelength[b.x1:b.x2], b.flux[b.x1:b.x2], filts=filts)
                if ~np.isnan(mag1) and ~np.isnan(mag2):
                    boot_1.append(mag1)
                    boot_2.append(mag2)
            diff = np.asarray(boot_1) - np.asarray(boot_2)
            boot_1m2s.append(diff)
        low_errors = []
        up_errors = []
        for clist in boot_1m2s:
            p = np.percentile(clist, [18, 50, 82])
            clower = p[1] - p[0]
            cupper = p[2] - p[1]
            low_errors.append(cupper)
            up_errors.append(clower)
        errors = [low_errors, up_errors]
        
    return phases, comp_1, comp_2, errors

