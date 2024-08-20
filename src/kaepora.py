import composite
import numpy as np
import sys
import argparse
import spectral_analysis as sa
import kaepora_plot as kplot
import warnings
from tabulate import tabulate
import matplotlib.pyplot as plt
import glob
warnings.filterwarnings("ignore")

"""This file contains various functions to interact with kaepora and facilitate
constructing composite spectra. Many of these are wrapper functions so that the
majority of the code can be accessed from one place. Most of the 
processing is provided by functions within composite.py.
"""

def set_min_num_spec(composites, num):
    """Modifies the x1 and x2 attributes of objects in composites to represent
        a wavelength range with >num contributing spectra.

        Args:
            composites: A list of spectrum objects output from
                kpora.make_composite()
            num: The minimum number of spectra required to define the wavelength
                range.
    """
    for comp in composites:
        # comp.spec_bin = np.array(comp.spec_bin)
        valid_range = np.where(np.asarray(comp.spec_bin) >= num)[0]
        comp.x1, comp.x2 = valid_range[0], valid_range[-1]
        # print comp.wavelength[valid_range]


def normalize_comps(composites, scale=1., w1=3500, w2=9000.):
    """Normalizes composite spectra in composites at their maximum values within
        wavelengths w1 and w2 to scale.

        Args:
            composites: An array of spectrum objects output from
                kpora.make_composite()
        Keyword Args:
            scale: Value to scale the maximum of each composite spectrum
            w1: Minimum wavelength
            w2: Maximum wavelength
        Returns:
            composites: The normalized array of composite spectra
    """
    for comp in composites:
        scale_range = np.where((comp.wavelength > w1) & (comp.wavelength < w2))
        # norm = scale/np.amax(comp.flux[comp.x1:comp.x2])
        norm = scale/np.nanmax(comp.flux[scale_range])
        comp.flux = norm*comp.flux

        if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
            comp.low_conf = norm*comp.low_conf
            comp.up_conf = norm*comp.up_conf
        comp.ivar /= (norm)**2
    return composites


def normalize_comp(comp):
    """Similar to normalize_comps() but takes a composite spectrum instead of 
    an array and scales the maximum to 1.
    """
    norm = 1./np.nanmax(comp.flux[comp.x1:comp.x2])
    comp.flux = norm*comp.flux

    if len(comp.RMSE) > 0:
        comp.RMSE = comp.RMSE*(norm)

    if len(comp.low_conf) > 0 and len(comp.up_conf) > 0:
        comp.low_conf = norm*comp.low_conf  
        comp.up_conf = norm*comp.up_conf
    comp.ivar /= (norm)**2
    return comp, norm

def grab(query, multi_epoch = True, make_corr = False, selection = 'max_coverage', grab_all=False, shape_param = 'dm15', verbose=False, db_file = None):
    """This function takes a SQL query and provides a list spectrum objects 
        (defined in composite.py) satisfy this query. 

        Args:
            query: The SQL query string

        Keyword Args:
            multi_epoch: If True, include all spectra for a given SN that satisify 
                the query. If False, choose one 1 spectrum per SN based on the 
                selection keyword.
            make_corr: If True, remove spectra that have been marked as 
                'questionable', peculiar events (Iax), and spectra that do not have 
                host extinction estimates.
            selection: If multi_epoch is False, this string defines the selection
                criteria for choosing a single spectrum from a SN. Options are:
                'maximum_coverage'(default): largest wavelength range
                'maximum_coverage_choose_uv': largest wavelength range but prioritize
                    hst and swift spectra
                'choose_bluest': smallest minimum wavelength
                'max_snr': highest signal to noise
                'accurate_phase': closest to middle of the phase bin. TODO: implement
                    this without parsing the query (currently requires uncommenting
                    code)
                'max_coverage_splice': allows multiple spectra from the same SN as 
                    long as overlap is < 500 A
            grab_all: If True, ignore other arguments and return all data that 
                satisfy the SQL query. This also ignores metadata and sets a very basic 
                list of spectrum attributes.

        Returns:
            An array of spectrum objects populated with metadata retrieved from 
            the SQL query.
    """

    if db_file is None:
        db_file = glob.glob('../data/*.db')[0]
        print ('Using: ' + db_file)

    # spec_array = composite.grab(query, multi_epoch = multi_epoch, 
    #                             make_corr = make_corr, grab_all=grab_all, db_file = db_file)
    spec_array = composite.new_grab(query, shape_param, make_corr = make_corr, db_file = db_file)
    spec_array = composite.prelim_norm(spec_array)
    if verbose:
        print ("Name", "Filename", "Source", "SNR", "Phase", "MJD", "MJD_max", "z", "Host Morphology", "Minwave", "Maxwave")
        for spec in spec_array:
            print (spec.name, spec.filename, spec.source, spec.SNR, spec.phase, spec.mjd, spec.event_data['Redshift'], spec.event_data['NED_host'], spec.wavelength[spec.x1], spec.wavelength[spec.x2])
            # print spec.name, spec.filename, spec.source, spec.SNR, spec.phase, spec.mjd, spec.redshift, spec.ned_host, spec.wavelength[spec.x1], spec.wavelength[spec.x2]

    return spec_array

def host_dereddening(SN_Array, r_v = 2.5, verbose=False, low_av_test=None, cutoff=2., norm=True):
    """Correct a list of spectrum objects for host galaxy extinction. This function assumes
    that the av_25, av_mlcs31, or av_mlcs17 spectrum attributes have values.

        Args:
            SN_Array: An array of spectrum objects.

        Keyword Args:
            r_v: The value of R_V used for the light curve fit
            verbose: (Deprecated)
            low_av_test: If True, does not correct SNe with A_V < low_av_test
            cutoff: SN_Array will only contain spectra from objects with A_V < cutoff

        Returns:
            An array of spectrum objects whose flux and ivar have been corrected
            for host extinction.
    """
    spec_array = composite.apply_host_corrections(SN_Array, r_v = r_v, verbose=verbose, low_av_test=low_av_test, cutoff=cutoff)
    if norm:
        spec_array = composite.prelim_norm(spec_array)
    return spec_array


def make_composite(query_strings, boot=False, nboots=100, medmean = 1, selection = 'max_coverage', gini_balance=False, verbose=True,
         make_corr=True, av_corr=True, scale_region=None,
         multi_epoch=True, combine=True, low_av_test=None, measure_vs = False, get_og_arr=False, shape_param='dm15', db_file=None):
    """ This is the main fuunction for constructing composite spectra from spectra stored in kaepora.
        Args:
            query_strings: A list of SQL query strings

        Keyword Args:
            boot: If 'b', estimate error via bootstrap resampling. If 'nb', create
                composite spectrum without bootstrapping (faster).
            medmean: If 1, do an inverse variance weighted average as a function 
                of wavelength. If 2, do a median.
            selection: If multi_epoch is False, this string defines the selection
                criteria for choosing a single spectrum from a SN. Options are:
                'maximum_coverage'(default): largest wavelength range
                'maximum_coverage_choose_uv': largest wavelength range but prioritize
                    hst and swift spectra
                'choose_bluest': smallest minimum wavelength
                'max_snr': highest signal to noise
                'accurate_phase': closest to middle of the phase bin. TODO: implement
                    this without parsing the query (currently requires uncommenting
                    code)
                'max_coverage_splice': allows multiple spectra from the same SN as 
                    long as overlap is < 500 A
            gini_balance: If True, finds region of the composite spectrum dominated 
                by high SNR spectra, then deweights these spectra to have the the 
                median weight in that wavelength range. TODO: improve this 
                algorithm (currently produced representative composites, but
                deweighting straight to median might not be ideal). If False, ivar
                spectra are left as is. 
            verbose: If True, print metadata of each spectrum object in SN_Array.
            multi_epoch: If True, include all spectra for a given SN that satisify 
                the query. If False, choose one 1 spectrum per SN based on the 
                selection keyword.
            combine: If True, combines spectra from the same SN with an inverse weighted
                average prior to construction fo the composite spectrum. This should be 
                True if doing an inverse variance weighted average.
            low_av_test: If True, does not correct SNe with A_V < low_av_test
            measure_vs: If True, will measure and print the Si II 6355 line velocity and 
                error of each composite spectrum
            og_arr: If True, this function will also return the original spectra in addition 
                to the new combined spectra of individual SNe.

        Returns:
            composites: Each element is a list of spectrum objects containing the composite spectrum.
            sn_arrays: Each element is a list of the spectrum objects used to create the composite spectrum.
                This includes pre-combined spectra but not necessarily every original spectrum.
            (optional) og_sn_arrays: Each element is a list of the original spectrum objects used to make 
                the composite spectra (only if og_arr=True). 
            boot_sn_arrays: Each element is a list spectrum objects containing composite spectra generated 
                from the bootstrap resampling process.
    """
    composites = []
    sn_arrays = []
    og_sn_arrays = []
    boot_sn_arrays = []
    store_boots = True
    num_queries = len(query_strings)
    for n in range(num_queries):
        if get_og_arr:
            comp, arr, og_arr, boots = composite.main(query_strings[n],boot=boot, nboots=nboots, medmean = medmean, shape_param=shape_param,
                                            selection = selection, gini_balance=gini_balance, combine=combine, db_file=db_file,
                                            make_corr=make_corr, av_corr=av_corr, scale_region=scale_region,
                                            verbose=verbose, multi_epoch=multi_epoch, low_av_test=low_av_test, get_og_arr=get_og_arr)
            og_sn_arrays.append(og_arr)
        else:
            comp, arr, boots = composite.main(query_strings[n],boot=boot, nboots=nboots, medmean = medmean, shape_param=shape_param,
                                            selection = selection, gini_balance=gini_balance, combine=combine, db_file=db_file,
                                            make_corr=make_corr, av_corr=av_corr, scale_region=scale_region,
                                            verbose=verbose, multi_epoch=multi_epoch, low_av_test=low_av_test, get_og_arr=get_og_arr)
        if store_boots:
            boot_sn_arrays.append(boots)
        composites.append(comp)
        sn_arrays.append(arr)

    if None not in composites:
        composite.optimize_scales(composites, composites[0], True)
        composites = normalize_comps(composites)
    if measure_vs:
        for comp in composites:
            dm15 = np.round(np.nanmean(comp.dm15_array[comp.x1:comp.x2]),2)
            # r = sa.measure_si_ratio(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], vexp = .001, dm15=dm15)
            v_strong, si_min_wave = sa.measure_velocity(comp.wavelength[comp.x1:comp.x2],comp.flux[comp.x1:comp.x2], 5900., 6300.)
            print ('v = ', v_strong)

        for boot in boot_sn_arrays:
            vs = []
            for b in boot:
                v_strong, si_min_wave = sa.measure_velocity(b.wavelength[b.x1:b.x2],b.flux[b.x1:b.x2], 5900., 6300.)
                vs.append(v_strong)
            v_err = np.nanstd(vs)
            print ('v_err = ', v_err)

    if get_og_arr:
        return composites, sn_arrays, og_sn_arrays, boot_sn_arrays
    else:
        return composites, sn_arrays, boot_sn_arrays

def save_spectra_as_flm(spectra, path = '', undo_z=False):

    for i, spec in enumerate(spectra):
        if undo_z and spec.event_data['Redshift'] != None:
            wave = spec.wavelength[spec.x1:spec.x2]*(1.+spec.event_data['Redshift'])
        else:
            wave = spec.wavelength[spec.x1:spec.x2]
        flux = spec.flux[spec.x1:spec.x2]
        err  = np.sqrt(1./spec.ivar[spec.x1:spec.x2])

        txt_header = 'wavelength flux fluxerr\n'
        txt_header = txt_header + 'Phase: {}\n'.format(spec.phase)
        txt_header = txt_header + 'MJD: {}\n'.format(spec.mjd)
        txt_header = txt_header + 'Redshift: {}\n'.format(spec.event_data['Redshift'])
        txt_header = txt_header + 'Velocity: {}\n'.format(spec.velocity)
        txt_header = txt_header + 'Velocity Error: {}\n'.format(spec.velocity_err)
        # txt_header = txt_header + 'AV_MW: {}\n'.format(spec.event_data['Av_MW'])
        print (i, path+spec.filename[:-4] + '-kaepora.flm')
        np.savetxt(path+spec.filename[:-4] + '-kaepora.flm',np.transpose([wave,flux,err]), header=txt_header)

def save_comp_and_boots(SN, boots, prefix, num_avg = 2, scale_region=None, folder = '../../'):
    if scale_region:
        normalize_comps([SN], scale=1., w1=scale_region[0], w2=scale_region[1])

    set_min_num_spec([SN], num_avg)
    phase = np.round(np.average(SN.phase_array[SN.x1:SN.x2]), 2)
    dm15 = np.round(np.average(SN.dm15_array[SN.x1:SN.x2]), 2)
    z = np.round(np.average(SN.red_array[SN.x1:SN.x2]), 3)
    num = np.amax(np.array(SN.spec_bin[SN.x1:SN.x2]))
    vel = np.round(np.average(SN.vel[SN.x1:SN.x2]), 2)
    print (phase, dm15, z)
    set_min_num_spec([SN], 1)

    if phase >= 0.:
        sign = 'p'
    else:
        sign = 'm'
    abs_phase = np.absolute(phase)
    phase_str = str(abs_phase)
    dm15_str = str(dm15)
    z_str = str(z)
    vel_str = str(vel)
    num_str = str(SN.num_sne)
    num_spec_str = str(SN.num_spec)

    file_path = folder + prefix + '_N=' + num_str + '_Nspec=' + num_spec_str + '_phase='+ sign + phase_str + '.txt'
    print (file_path)

    data = []

    compwave = np.array(SN.wavelength[SN.x1:SN.x2])
    compflux = np.array(SN.flux[SN.x1:SN.x2])
    
    data.append(compwave)
    data.append(compflux)

    plt.figure(figsize=[15,8])
    for boot in boots:
        wave = np.array(boot.wavelength[SN.x1:SN.x2])
        flux = np.array(boot.flux[SN.x1:SN.x2])
        plt.plot(wave,flux, 'g')
        data.append(boot.flux[SN.x1:SN.x2])
    plt.plot(compwave, compflux, 'k')
    plt.show()

    np.savetxt(file_path, np.transpose(np.asarray(data)))

def create_composites_from_files():
    s19_composites = glob.glob('../data/S19_Composite_spectra/siebert_phase*.txt')
    composites = {}
    for comp_file in s19_composites:
        data = np.genfromtxt(comp_file, unpack=True, skip_header=9)
        low_conf = data[2]
        up_conf = data[3]
        var_est = (low_conf + up_conf)/2.
        spec = composite.spectrum(wavelength = data[0], flux=data[1], ivar=(1./var_est)**2.)
        with open(comp_file) as comp:
            phase = comp.readlines()[3].split()[-1]
            spec.phase = float(phase)
            spec.mjd = float(phase)
        spec.spec_bin = data[7]
        if float(phase) < 90:
            composites[float(phase)] = spec
    return composites

def create_query_strings_for_sequence(phase_range, dm15_range, velocity_range=None, phase_step = 1, phase_bin_size=3, phase_step_bin_scale = True):

    if velocity_range:
        v_beg = velocity_range[0]
        v_end = velocity_range[1]

    p_beg = phase_range[0]
    p_end = phase_range[1]
    dm_beg = dm15_range[0]
    dm_end = dm15_range[1]
    query_strings = ["SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN and " \
                     "phase < {p_beg} and ((dm15_source >= {dm_beg} and dm15_source < {dm_end}) or " \
                     "(dm15_from_fits >= {dm_beg} and dm15_from_fits < {dm_end})) and Spectra.SN != '2011fe' and Spectra.SN != '1994d'".format(p_beg=str(p_beg), dm_beg=str(dm_beg), dm_end=str(dm_end))]
    p1 = p_beg
    k = phase_bin_size
    p2 = p1+k

    p_scale1 = 30
    p_scale2 = 60
    while p1 < p_end:
        query_strings.append("SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN and " \
                              "phase >= {p1} and phase < {p2} and ((dm15_source >= {dm_beg} and dm15_source < {dm_end}) or " \
                              "(dm15_from_fits >= {dm_beg} and dm15_from_fits < {dm_end})) and Spectra.SN != '2011fe' and Spectra.SN != '1994d'".format(p1=str(p1), p2=str(p2), dm_beg=str(dm_beg), dm_end=str(dm_end)))
        if phase_step_bin_scale:
            if p1 > p_scale1 and p1 < p_scale2:
                p1+=phase_step
                p2 = p1+7
            elif p1 > p_scale2:
                p1 += phase_step
                p2 = p1+14
            else:
                p1 += phase_step
                p2 += phase_step

        else:
            p1 += phase_step
            p2 += phase_step

    if velocity_range:
        query_strings_v = []
        for q in query_strings:
            v_query = q = q + "and v_at_max between {v_beg} and {v_end}".format(v_beg=str(v_beg), v_end=str(v_end))
            query_strings_v.append(v_query)
        query_strings = query_strings_v
    for q in query_strings:
        print (q)
    return query_strings

def save_comps_to_files(composites, prefix, num_avg = 5, boot=True, scale_region=None, folder = '../../'):
    """Saves the data contained in a composite spectrum object to a text file.
        Args:
            composites: A list of composite spectrum objects
            prefix: prefix to be added to the filename

        Keyword Args:
            num_avg: the minimum number of spectra required to determine the 
                average properties that are included in the filename

    """
    if scale_region:
        normalize_comps(composites, scale=1., w1=scale_region[0], w2=scale_region[1])

    for SN in composites:
        set_min_num_spec(composites, num_avg)
        phase = np.round(np.average(SN.phase_array[SN.x1:SN.x2]), 2)
        dm15 = np.round(np.average(SN.dm15_array[SN.x1:SN.x2]), 2)
        z = np.round(np.average(SN.red_array[SN.x1:SN.x2]), 3)
        num = np.amax(np.array(SN.spec_bin[SN.x1:SN.x2]))
        vel = np.round(np.average(SN.vel[SN.x1:SN.x2]), 2)
        print (phase, dm15, z)
        set_min_num_spec(composites, 1)

        if phase >= 0.:
            sign = 'p'
        else:
            sign = 'm'
        abs_phase = np.absolute(phase)
        phase_str = str(abs_phase)
        dm15_str = str(dm15)
        z_str = str(z)
        vel_str = str(vel)
        num_str = str(SN.num_sne)
        num_spec_str = str(SN.num_spec)

        file_path = folder + prefix + '_N=' + num_str + '_Nspec=' + num_spec_str + '_phase='+ sign + phase_str + '_x1=' + dm15_str + '.txt'
        print (file_path)
        with open(file_path, 'w') as file:
            file.write('# SQL Query: ' + SN.query + '\n')
            file.write('# N = ' + num_str + '\n')
            file.write('# N_spec = ' + num_spec_str + '\n')
            file.write('# phase = ' + str(phase) + '\n')
            file.write('# x1 = ' + dm15_str + '\n')
            # file.write('# dm15 = ' + dm15_str + '\n')
            file.write('# vel = ' + vel_str + '\n')
            file.write('# z = ' + z_str + '\n')
            file.write('\n')

            wave = np.array(SN.wavelength[SN.x1:SN.x2])
            flux = np.array(SN.flux[SN.x1:SN.x2])
            if boot:
                up_conf = np.array(SN.up_conf[SN.x1:SN.x2]) - flux
                low_conf = flux - np.array(SN.low_conf[SN.x1:SN.x2])

            phase_arr = np.array(SN.phase_array[SN.x1:SN.x2])
            dm15_arr = np.array(SN.dm15_array[SN.x1:SN.x2])
            red_arr = np.array(SN.red_array[SN.x1:SN.x2])
            vel_arr = np.array(SN.vel[SN.x1:SN.x2])
            spec_per_sn_arr = np.array(SN.spec_bin[SN.x1:SN.x2])
            if boot:
                data = np.c_[wave,flux,low_conf,up_conf,phase_arr,dm15_arr, red_arr, vel_arr, spec_per_sn_arr]
                # table = tabulate(data, headers=['Wavelength', 'Flux', '1-Sigma Lower', '1-Sigma Upper', 'Phase', 'Dm15', 'Redshift', 'SNe per Bin'], 
                #                                     tablefmt = 'ascii')
                table = tabulate(data, headers=['Wavelength', 'Flux', '1-Sigma Lower', '1-Sigma Upper', 'Phase', 'x1', 'Redshift', 'Velocity', 'SNe per Bin'], 
                                                    tablefmt = 'ascii')
            else:
                data = np.c_[wave,flux,phase_arr,dm15_arr, red_arr, vel_arr, spec_per_sn_arr]
                # table = tabulate(data, headers=['Wavelength', 'Flux', 'Phase', 'Dm15', 'Redshift', 'SNe per Bin'], 
                #                                     tablefmt = 'ascii')
                table = tabulate(data, headers=['Wavelength', 'Flux', 'Phase', 'x1', 'Redshift', 'Velocity','SNe per Bin'], 
                                                    tablefmt = 'ascii')
            file.write(table)

if __name__ == "__main__":
    composites = []
    SN_Arrays = []
    boot_sn_arrays = []
    store_boots = True

    boot = sys.argv[1]
    if boot == 'nb':
        boot = False
    else:
        boot = True

    query_strings = sys.argv[2:]

    num_queries = len(query_strings)

    for n in range(num_queries):
        c, sn_arr, boots = composite.main(query_strings[n], boot, nboots=100, medmean=1, gini_balance = True, make_corr=False, multi_epoch=True, combine=True)
        # c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = False, multi_epoch=True, combine=False) 
        # composites.append(c)
        # SN_Arrays.append(sn_arr)
        # if store_boots:
        #   boot_sn_arrays.append(boots)

        #use this for good composites
        # c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = True, verbose=True, multi_epoch=True, combine=True)
        # composites.append(c)
        # SN_Arrays.append(sn_arr)
        # if store_boots:
        #   boot_sn_arrays.append(boots)
            
        # c, sn_arr, boots = composite.main(query_strings[n], boot, medmean=1, gini_balance = True, verbose=True, multi_epoch=True, combine=True)
        composites.append(c)
        SN_Arrays.append(sn_arr)
        if store_boots:
            boot_sn_arrays.append(boots)

    for i, comp in enumerate(composites):
        dm15s = []
        kplot.plot_comp_and_all_spectra(comp, SN_Arrays[i], show_ivar=True)
        # for SN in SN_Arrays[i]:
        #   if SN.dm15_source is not None:
        #       dm15s.append(SN.dm15_source)
        #   else:
        #       dm15s.append(SN.dm15_from_fits)
        # plt.hist(dm15s)
        # plt.show()
    # comp2 = composites[1]
    # for sn in sn_arr:
    #   r = sa.measure_si_ratio(sn.wavelength[sn.x1:sn.x2], sn.flux[sn.x1:sn.x2])
    #   print r

    # for SN in SN_Arrays[0]:
    #   # if SN.name not in bad_vels and SN.source != 'swift_uv':
    #   if SN.name and SN.source != 'swift_uv':
    #       var = 1./SN.ivar
    #       vexp, SNR = sa.autosmooth(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], var_y = var[SN.x1:SN.x2])
    #       v, si_min_wave = sa.measure_velocity(SN.wavelength[SN.x1:SN.x2], SN.flux[SN.x1:SN.x2], 5900., 6300., vexp=vexp, plot=True)
    #       print SN.name, v, si_min_wave


    # for i, comp in enumerate(composites):
    #   comp.name = "Comp" + str(i)
    #   dm15 = np.round(np.nanmean(comp.dm15_array[comp.x1:comp.x2]),2)
    #   # r = sa.measure_si_ratio(comp.wavelength[comp.x1:comp.x2], comp.flux[comp.x1:comp.x2], vexp = .001, dm15=dm15)
    #   v_strong, si_min_wave = sa.measure_velocity(comp.wavelength[comp.x1:comp.x2],comp.flux[comp.x1:comp.x2], 5900., 6300.)
    #   print 'v = ', v_strong


    # vs = []
    # for b in boot_sn_arrays[0]:
    #   v_strong, si_min_wave = sa.measure_velocity(b.wavelength[b.x1:b.x2],b.flux[b.x1:b.x2], 5900., 6300.)
    #   vs.append(v_strong)
    # v_err = np.nanstd(vs)
    # print 'v_err = ', v_err


    # r = sa.measure_si_ratio(comp2.wavelength[comp2.x1:comp2.x2], comp2.flux[comp2.x1:comp2.x2], vexp = .001)
    # print r
    # set_min_num_spec(composites, 10)
    set_min_num_spec(composites, 2)
    # set_min_num_spec(composites, 20)
    normalize_comps(composites)
    # for comp in composites:
    #   sa.measure_si_velocity(comp)
    
    #plotting routines to visualize composites
    # set_min_num_spec(composites, 5)
    kplot.comparison_plot(composites)
    # scaled_plot(composites, min_num_show=2, min_num_scale=2)
    # stacked_plot(composites)

    # save_comps_to_files(composites)
    
        