=====================
Querying the Database
=====================

The database currently consists of two tables: *Spectra* and *Events*. These tables host the spectrum-specific and SN-specific metadata respectively. Currently the only way to interact with the database is via an SQL join on these tables. Right now, you should run your code from the /src directory. Our documentation will focus on the routines available in the ``kaepora.py``. Start by importing this module:

.. code-block:: python

    import kaepora as kpora

You can then define an array containing SQL queries and obtain spectra from the database. For example:

.. code-block:: python

    example_query = ["SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN where phase >= -1 and phase <= 1 and ((dm15_source < 1.8) or (dm15_from_fits < 1.8))"]
    spec_array = kpora.grab(example_query[0])

If you would like to remove atypical SNe Ia, SNe with flagged artifacts, and SNe with poor host reddening corrections use:

.. code-block:: python

    spec_array = kpora.grab(example_query[0], make_corr=True)

These spectra have been already been corrected for MW reddening. To correct these spectra for host-galaxy reddening (and exclude SNe with A\ :sub:`V` > 2) with a F99 reddening law use:

.. code_block:: python

    spec_array = kpora.host_dereddening(spec_array, cutoff=2.)

================
Spectrum Objects
================

``spec_array`` now contains an array of objects that contain our homogenized spectra and all of the spectrum- and SN-specific metadata. Currently these objects are made to represent single spectra, so objects generated from the same SNe will contain some redundant SN metadata. These spectra are normalized to their maximum flux. Basic information on these objects can be viewed with:

.. code_block:: python

    for spec in spec_array_dered:
        print spec.name, spec.filename, spec.source, spec.phase, spec.wavelength[spec.x1], spec.wavelength[spec.x2]

A spectrum and its variance can be plotted with:

.. code_block:: python

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,1)
    example_spec = spec_array_dered[20]
    ax[0].plot(example_spec.wavelength, example_spec.flux)
    ax[1].plot(example_spec.wavelength, 1/example_spec.ivar)
    plt.show()

Below we describe other attributes of these objects that are also queryable parameters of the database.

======
Schema
======

Spectral Attributes
===================

+-----------+------------+--------------------------------------------+--------+
| Attribute | SQL Format | Description                                | Type   |
+-----------+------------+--------------------------------------------+--------+
| name      | "SN"       | SN name                                    | String |
+-----------+------------+--------------------------------------------+--------+
| filename  | "Filename" | Filename from data source                  | String |
+-----------+------------+--------------------------------------------+--------+
| source    | "Source"   | Data source                                | String |
+-----------+------------+--------------------------------------------+--------+
| minwave   | "Minwave"  | Minimum wavelength of original spectrum    | float  |
+-----------+------------+--------------------------------------------+--------+
| maxwave   | "Maxwave"  | Maximum wavelength of original spectrum    | float  |
+-----------+------------+--------------------------------------------+--------+
| SNR       | "snr"      | Median S/N of the spectrum                 | float  |
+-----------+------------+--------------------------------------------+--------+
| mjd       | "MJD"      | Modified Julian Date of the spectrum       | float  |
+-----------+------------+--------------------------------------------+--------+
| phase     | "Phase"    | Rest-frame days from B-Band maximum        | float  |
+-----------+------------+--------------------------------------------+--------+
| ref       | "Ref"      | Bibtex code                                | String |
+-----------+------------+--------------------------------------------+--------+

SN Attributes
=============
These attributes contain the most metadata. We also include (but do not list) metadata from the results of several different light curve fits. If you would like to construct a query based on these metadata please contact me. 

+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| Attribute         | SQL Format          | Description                                                                          | Type   |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| redshift          | "Redshift"          | Redshift from NED                                                                    | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| mjd_max           | "MJD_max"           | Modified Julian date corresponding to the time of maximum-light                      | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| dm15_source       | "Dm15_source"       | Dm15 from the source survey                                                          | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| dm15_from_fits    | "Dm15_from_fits"    | Dm15 calculated from the polynomial relationship with a light-curve shape parameter  | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| e_dm15            | "e_dm15"            | Error in dm15                                                                        | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| av_25             | "Av_25"             | Estimated host galaxy extinction from an MLCS fit using R_v = 2.5                    | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| m_b_cfa           | "M_b_cfa"           | Absolute B-Band magnitude at maximum light from the CfA sample                       | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| m_b_cfa_err       | "M_b_cfa_err"       | Error in m_b_cfa                                                                     | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| b_minus_v_cfa     | "B_minus_V_cfa"     | B-V color at maximum light                                                           | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| b_minus_v_cfa_err | "B_minus_V_cfa_err" | Error in b_minus_v_cfa                                                               | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| v_at_max          | "V_at_max"          | Estimated velocity at maximum light                                                  | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| v_err             | "V_err"             | Error in v_at_max                                                                    | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| ned_host          | "NED_host"          | A simplified version of NED host galaxy morphology based on cross-listed objects     | String |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| carbon            | "Carbon_presence"   | 'A': carbon detected, 'F'; marginal carbon detected, 'N': no carbon detected         | String |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+
| hubble_res        | "Hubble_res"        | Hubble residual from a SALT fit                                                      | float  |
+-------------------+---------------------+--------------------------------------------------------------------------------------+--------+

You can view all attributes of the spectrum object with the code below:

.. code-block:: python

    spec_attributes = dir(spec_array[0])
    print spec_attributes