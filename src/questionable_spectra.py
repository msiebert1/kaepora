def bad_files():
		bad_files =        ['SN05hc_051018_r01_NTT_EM.dat', '2003du_20030501_4066_11015_00.dat',
							'2002er_20020901_3213_9175_00.dat', '2003du_20030503_4070_10994_00.dat',
							'2003du_20030512_4070_11015_00.dat', 'sn2006nz-20061124-ntt.dat',
							'sn2001eh-20010925.757-hst.flm', 'sn2001ep-20011102.871-hst.flm',
							'2005cf_20050601_3243_9720_00.dat', 'sn2006cm-20060529.46-fast.flm',
							'2002bo_20020321_3357_7725_00.dat', 'sn1994S-19940612.26-mmt.flm',
							'sn2003cg-20030329.27-mmt.flm', 'sn1995ac-19950929.27-fast.flm',
							'sn2007hj-20070903.28-fast.flm', '2000E_20000131_3274_7349_00.dat',
							'sn2006cj-20060521.29-fast.flm', 'sn2006oa-20061125.08-fast.flm',
							'sn2005cf-20050609.5-uvot-clip.flm', 'sn2006kf-20061030.385-ui.flm',
							'SN07bd_070412_b01_DUP_WF.dat', 'SN09ad_090223_b01_DUP_WF.dat',
							'SN05kc_051124_b01_DUP_MS.dat', 'sn2005eq-20051011.376-ui-corrected.flm',
							'sn2006et-20060919.345-ui.flm', 'sn2007cq-20070623.431-ui.flm',
							'sn1997bq-19970408.14-mmt.flm', 'sn2006lf-20061028.51-fast.flm',
							'sn2005eq-20051002.51-fast.flm', 'sn1995bd-19951223.34-fast.flm',
							'sn1998ab-19980403.38-fast.flm', 'sn1994M-19940612.22-mmt.flm',
							'2006X_20060209_3834_8139_00.dat', '2003du_20030508_4066_10997_00.dat',
							'SN05ke_051125_b01_T60_CS.dat', 'sn1994s-19940616-uoi.flm',
							'1994D_19940317_2999_10549_00.dat', '2005cf_20050608_3365_9997_00.dat',
							'2002dj_20020620_3210_9231_00.dat', 'sn2006lf-20061029.40-fast.flm',
							'sn2006lf-20061030.41-fast.flm', 'sn2006lf-20061031.42-fast.flm',
							'sn2006lf-20061101.40-fast.flm', 'sn2006lf-20061111.36-fast.flm',
							'sn2006lf-20061112.37-fast.flm','sn2006lf-20061113.36-fast.flm',
							'sn2006lf-20061115.42-fast.flm','sn2006lf-20061116.39-fast.flm',
							'sn2006lf-20061117.42-fast.flm','sn2006lf-20061119.35-fast.flm',
							'sn2006lf-20061122.35-fast.flm','sn2006lf-20061125.34-fast.flm',
							'sn2006lf-20061214.28-fast.flm','sn2006lf-20061216.30-fast.flm',
							'sn2006lf-20061226.26-fast.flm','sn2006lf-20061227.29-fast.flm',
							'sn1995al-19951114.52-fast.flm', '2002er_20020901_3213_9175_00.dat',
							'sn2003kf-20031216.37-fast.flm','SN05bg_050419_b01_DUP_WF.dat',
							'2003du_20040202_3182_9333_00.dat', '2002bo_20020328_3099_8801_00.dat',
							'2005cf_20050603_3721_8786_00.dat', 'SN06hx_061005_g01_NTT_EM.dat',
							'2003du_20030429_3428_9436_00.dat', '2000E_20000127_3213_7513_00.dat',
							'sn2006bt-20060427.629-br.flm','sn2003ic-20030919.40-fast.flm',
							'2000E_20000130_3274_7356_00.dat','sn2006cm-20060528.45-fast.flm',
							'sn2006nz-20061117.18-fast.flm','sn1996ai-19960620-uo.flm',
							'sn2006cm-20060528.424-ui.flm','2003du_20030506_4068_10991_00.dat',
							'sn2005A-20050111.04-ldss2.flm','sn2003iv-20031023.469-ui.flm',
							'SN07al_070314_b01_DUP_BC.dat', 'sn2006ke-20061024.397-ui.flm',
							'sn2002bf-20020307-os.flm','sn1996bk-19961015.10-fast.flm',
							'SN06mr_061113_r01_BAA_IM.dat','sn2003ic-20030927.38-fast.flm',
							'SN08hv_081220_b01_DUP_WF.dat','sn2002de-20020614.19-fast.flm',
							'sn2003ic-20030929.35-fast.flm','sn2006cj-20060529.26-fast.flm',
							'SN06mr_061116_b01_DUP_WF.dat','SN07al_070319_b01_DUP_BC.dat',
							'2003du_20030515_4969_9241_00.dat','SN07N_070131_b01_NTT_EM.dat',
							'sn2006ke-20061030.525-ui.flm','sn2003ic-20031002.25-fast.flm',
							'sn2005hf-20051027.27-fast.flm', 'SN06mr_061119_b01_DUP_WF.dat',
							'sn2005hf-20051028.22-fast.flm','SN06mr_061119_b01_HIL_BC.dat',
							'sn2001ic-20011211-ui.flm','sn2005mz-20060121.13-fast.flm',
							'sn2000cn-20000623.34-fast.flm','sn2007al-20070320.11-ldss3.flm',
							'sn2007ax-20070407.16-fast.flm','SN07on_071125_b01_CLA_MA.dat', 
							'SN08hv_081227_b01_DUP_WF.dat','sn2005hf-20051030.22-fast.flm',
							'sn2007al-20070321.33-fast.flm','sn1999gh-19991212.50-fast.flm',
							'sn2001N-20010202.47-fast.flm','sn2002bf-20020316.32-fast.flm',
							'SN06mr_061122_b01_DUP_WF.dat','sn2005hf-20051031.24-fast.flm',
							'sn1998bp-19980516.44-fast.flm','sn2003hu-20031002.13-fast.flm',
							'sn2005hf-20051101.22-fast.flm','SN06bd_060330_b01_DUP_WF.dat',
							'sn2007ci-20070609.244-ui-corrected.flm','SN07al_070326_b01_DUP_BC.dat',
							'sn1998bp-19980518.36-fast.flm','SN05kc_051208_b01_T60_CS.dat',
							'sn1999gh-19991217-ui.flm','sn2007bc-20070426.387-ui-corrected.flm',
							'sn2002do-20020704.39-fast.flm','SN07on_071203_b01_DUP_WF.dat',
							'sn2005mz-20060129.11-fast.flm','SN06gt_061013_b01_DUP_WF.dat',
							'sn1998bp-19980522-r.flm','sn2005hf-20051106.34-fast.flm',
							'sn2006cm-20060620.410-ui.flm','sn2000cn-20000704.32-fast.flm',
							'SN08R_080219feb08_b01_CLA_MA.dat','SN07jg_071016_b01_DUP_BC.dat',
							'2003du_20030530_4060_10974_00.dat','sn2006nz-20061213-ntt.dat',
							'sn2006te-20070126.326-ui.flm','SN08ia_090122_b01_CLA_LD.dat',
							'sn2001G-20010423.19-fast.flm','SN06X_060524_b01_DUP_BC.dat',
							'sn2006X-20060221.40-fast.flm','2006X_20060219_3731_8515_00.dat',
							'2006X_20060221_3981_8865_00.dat','sn1999cl-19990614.18-fast.flm',
							'sn2006X-20060222.41-fast.flm','sn2006x-20060222.413-ui.flm',
							'2006X_20060225_3734_8223_00.dat','sn2006X-20060225.36-fast.flm',
							'2006X_20060227_3918_8203_00.dat','sn2006X-20060227.44-fast.flm',
							'sn2006X-20060228.34-fast.flm','2002er_20020916_3336_8734_00.dat',
							'sn2006X-20060302.47-fast.flm','2006X_20060304_3783_8272_00.dat',
							'sn2006X-20060304.51-fast.flm','2006X_20060307_3861_8130_00.dat',
							'sn2006X-20060309.30-fast.flm','2002er_20020926_3489_8768_00.dat',
							'2002er_20021010_3560_9363_00.dat','1996X_19960613_2806_10203_00.dat',
							'sn2005m-20050201.556-ui.flm','2002er_20020906_3480_10263_00.dat',
							'sn2003cg-20030402.17-ldss2.flm']

		return bad_files
		#SN05hc_051018_r01_NTT_EM.dat very noisy
		#2003du_20030501_4066_11015_00.dat very large negative value
		#2002er_20020901_3213_9175_00.dat gap in spectrum
		#2003du_20030503_4070_10994_00.dat very strong emission line?
		#2003du_20030512_4070_11015_00.dat very negative values
		#sn2006nz-20061124-ntt.dat very negative values
		#sn2001eh-20010925.757-hst.flm some interpolated sections
		#sn2001ep-20011102.871-hst.flm causes problem not sure why
		#2005cf_20050601_3243_9720_00.dat causes problem not sure why
		#sn2006cm-20060529.46-fast.flm weird snr
		#2002bo_20020321_3357_7725_00.dat doesn't scale properly, flux very small at low wavelength
		#sn1994S-19940612.26-mmt.flm silicon line was interpolated
		#sn2003cg-20030329.27-mmt.flm SNR above 600 biases averaging
		#sn1995ac-19950929.27-fast.flm noisy
		#sn2007hj-20070903.28-fast.flm some interpolated sections
		#2000E_20000131_3274_7349_00.dat spectrum seems to be blueshifted
		#sn2006cj-20060521.29-fast.flm very noisy and some interpolation
		#sn2006oa-20061125.08-fast.flm some interpolated sections
		#sn2005cf-20050609.5-uvot-clip.flm high snr drags composite down
		#sn2006kf-20061030.385-ui.flm some interpolated sections
		#SN07bd_070412_b01_DUP_WF.dat some interpolated sections
		#SN09ad_090223_b01_DUP_WF.dat some interpolated sections
		#SN05kc_051124_b01_DUP_MS.dat some interpolated sections
		#sn2006et-20060919.345-ui.flm some interpolated sections
		#sn2007cq-20070623.431-ui.flm some interpolated sections
		#sn1997bq-19970408.14-mmt.flm telluric?
		#sn2006lf-20061028.51-fast.flm some interpolated sections, variance spectrum seems wrong
		#sn2005eq-20051002.51-fast.flm some interpolated sections
		#sn1995bd-19951223.34-fast.flm some interpolated sections
		#sn1998ab-19980403.38-fast.flm some interpolated sections
		#sn1994M-19940612.22-mmt.flm large interpolated section
		#2006X_20060209_3834_8139_00.dat host correction causes
		#2003du_20030508_4066_10997_00.dat very large negative value
		#SN05ke_051125_b01_T60_CS.dat some interpolated sections
		#sn1994s-19940616-uoi.flm some interpolated sections
		#1994D_19940317_2999_10549_00.dat not joined properly, telluric absorption
		#2005cf_20050608_3365_9997_00.dat telluric absorption
		#2002dj_20020620_3210_9231_00.dat telluric absorption
		#2002er_20020901_3213_9175_00.dat large gap
		#
		#All 2006lf cfa data seems to have large ivar (skews data)
		#sn2003kf-20031216.37-fast.flm seems to have large ivar
		#SN05bg_050419_b01_DUP_WF.dat ivar very large at large wavelengths
		#2003du_20040202_3182_9333_00.dat variance blows up
		#2002bo_20020328_3099_8801_00.dat large interpolated section
		#2005cf_20050603_3721_8786_00.dat large interpolated section
		#SN06hx_061005_g01_NTT_EM.dat slope is off after MW correction
		#2003du_20030429_3428_9436_00.dat telluric absorption
		#2000E_20000127_3213_7513_00.dat weird slope and noise
		#sn2006bt-20060427.629-br.flm slope is off after MW correction
		#sn2003ic-20030919.40-fast.flm slope is off after MW correction
		#2000E_20000130_3274_7356_00.dat slope is off after MW correction
		#sn2006cm-20060528.45-fast.flm slope is off after MW correction
		#sn2006nz-20061117.18-fast.flm slope is off after MW correction
		#sn1996ai-19960620-uo.flm large interpolated section
		#sn2006cm-20060528.424-ui.flm slope is off affter MW correction
		#2003du_20030506_4068_10991_00.dat large negative value
		#sn2005A-20050111.04-ldss2.flm slope is off after MW correction
		#sn2003iv-20031023.469-ui.flm slope is off after MW correction
		#SN07al_070314_b01_DUP_BC.dat slope is off after MW correction
		#sn2006ke-20061024.397-ui.flm slope is off after MW correction
		#sn2002bf-20020307-os.flm large blueshift?
		#sn1996bk-19961015.10-fast.flm slope is off after MW correction
		#SN06mr_061113_r01_BAA_IM.dat weird shape
		#sn2003ic-20030927.38-fast.flm slope is off after MW correction
		#SN08hv_081220_b01_DUP_WF.dat slope is off after MW correction
		#sn2002de-20020614.19-fast.flm slope is off after MW correction
		#sn2003ic-20030929.35-fast.flm slope is off after MW correction
		#sn2006cj-20060529.26-fast.flm slope is off after MW correction
		#SN06mr_061116_b01_DUP_WF.dat continuum seems off
		#SN07al_070319_b01_DUP_BC.dat continuum seems off
		#2003du_20030515_4969_9241_00.dat large blueshift?
		#SN07N_070131_b01_NTT_EM.dat continuum seems off
		#sn2006ke-20061030.525-ui.flm continuum seems off
		#sn2003ic-20031002.25-fast.flm slope is off after MW correction
		#sn2005hf-20051027.27-fast.flm slope is off after MW correction
		#SN06mr_061119_b01_DUP_WF.dat continuum seems off
		#sn2005hf-20051028.22-fast.flm slope is off after MW correction
		#SN06mr_061119_b01_HIL_BC.dat continuum seems off
		#sn2001ic-20011211-ui.flm almost entirely flat
		#sn2005mz-20060121.13-fast.flm slope is off after MW correction
		#sn2000cn-20000623.34-fast.flm slope is off after MW correction
		#sn2007al-20070320.11-ldss3.flm continuum seems off
		#sn2007ax-20070407.16-fast.flm slope is off after MW correction
		#SN07on_071125_b01_CLA_MA.dat continuum seems off, telluric contamination
		#SN08hv_081227_b01_DUP_WF.dat continuum seems off
		#sn2005hf-20051030.22-fast.flm spectrum is flat
		#sn2007al-20070321.33-fast.flm continuum seems off
		#sn1999gh-19991212.50-fast.flm continuum seems off
		#sn2001N-20010202.47-fast.flm slope is off of MW correction
		#sn2002bf-20020316.32-fast.flm slope is off of MW correction
		#SN06mr_061122_b01_DUP_WF.dat continuum seems off
		#sn2005hf-20051031.24-fast.flm spectrum is flat
		#sn1998bp-19980516.44-fast.flm continuum seems off
		#sn2003hu-20031002.13-fast.flm spectrum is flat
		#sn2005hf-20051101.22-fast.flm continuum seems off
		#SN06bd_060330_b01_DUP_WF.dat continuum seems off
		#sn2007ci-20070609.244-ui-corrected.flm continuum seems off
		#SN07al_070326_b01_DUP_BC.dat continuum seems off
		#sn1998bp-19980518.36-fast.flm continuum seems off
		#SN05kc_051208_b01_T60_CS.dat variance spectrum is wrong
		#sn1999gh-19991217-ui.flm continuum seems off 
		#sn2007bc-20070426.387-ui-corrected.flm continuum seems off
		#sn2002do-20020704.39-fast.flm spectrum is flat
		#SN07on_071203_b01_DUP_WF.dat continuum seems off
		#sn2005mz-20060129.11-fast.flm continuum seems off
		#SN06gt_061013_b01_DUP_WF.dat continuum seems off
		#sn1998bp-19980522-r.flm continuum seems off
		#sn2005hf-20051106.34-fast.flm spectrum is flat
		#sn2006cm-20060620.410-ui.flm continuum seems off
		#sn2000cn-20000704.32-fast.flm continuum seems off
		#SN08R_080219feb08_b01_CLA_MA.dat continuum seems off, telluric absorption
		#SN07jg_071016_b01_DUP_BC.dat very noisy, many interpolated sections
		#2003du_20030530_4060_10974_00.dat large negative values
		#sn2006nz-20061213-ntt.dat large negative values
		#sn2006te-20070126.326-ui.flm spectrum is flat
		#SN08ia_090122_b01_CLA_LD.dat redshifted?
		#sn2001G-20010423.19-fast.flm continuum seems off
		#SN06X_060524_b01_DUP_BC.dat doesnt scale properly
		#sn2006X-20060221.40-fast.flm blueshifted?
		#2006X_20060219_3731_8515_00.dat blueshifted?
		#2006X_20060221_3981_8865_00.dat blueshifted?
		#sn1999cl-19990614.18-fast.flm variance seems wrong
		#sn2006X-20060222.41-fast.flm blueshifted?
		#sn2006x-20060222.413-ui.flm blueshifted?
		#2006X_20060225_3734_8223_00.dat blueshifted?
		#sn2006X-20060225.36-fast.flm blueshifted?
		#2006X_20060227_3918_8203_00.dat blueshifted?
		#sn2006X-20060227.44-fast.flm blueshifted?
		#sn2006X-20060228.34-fast.flm blueshifted?
		#2002er_20020916_3336_8734_00.dat seems off near red edge
		#sn2006X-20060302.47-fast.flm blueshifted?
		#2006X_20060304_3783_8272_00.dat blueshifted?
		#sn2006X-20060304.51-fast.flm blueshifted?
		#2006X_20060307_3861_8130_00.dat blueshifted?
		#sn2006X-20060309.30-fast.flm blueshifted?
		#2002er_20020926_3489_8768_00.dat continuum seems off
		#2002er_20021010_3560_9363_00.dat continuum seems off
		#1996X_19960613_2806_10203_00.dat blueshifted?
		#sn2005m-20050201.556-ui.flm slope seems off
		#2002er_20020906_3480_10263_00.dat slope seems off
		#sn2003cg-20030402.17-ldss2.flm slope seems off

		# sn1994d-19940603.flm ??

		#sn2001az-20010430-ui.flm

		#Weird continuum
		#sn1997E-19970116.25-fast.flm   0.1221
		#sn1996bl-19961018.18-fast.flm  0.2945
		#sn1998ef-19981024.30-fast.flm  -0.0303
		#sn2003it-20031019.18-fast.flm  0.0932

		