;Create perigee tplot save files for Wygant


;SWEAP data is ~ 2.6 Hz rate.
;FIELDS data I have here in 1 min data and 4 Hz data.

path = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/solar_probe_routines/PSP_perigee_tplot_files_wygant/'
tplot_restore,filename=path + 'psp_fields_vars_jan2020.tplot'
tplot_restore,filename=path + 'psp_sweap_vars_jan2020.tplot'

;tplot_restore,filename=path + 'psp_fields_vars_apr2019.tplot'
;tplot_restore,filename=path + 'psp_sweap_vars_apr2019.tplot'
;tplot_restore,filename=path + 'psp_fields_vars_nov2018.tplot'
;tplot_restore,filename=path + 'psp_sweap_vars_nov2018.tplot'


;FIELDS data interplated to 1 min cadence ("_interp" quantities interpolated to the times of the flags variable)
;Ephem and the 1 min data are on basically the same time base, but with slight differences
vars = ['psp_fld_l2_mag_RTN_1min_interp',$               ;RTN 1 min mag
'psp_fld_l2_mag_RTN_1min_deltaT_interp',$        ;RTN 1 min mag dT
'psp_fld_l2_mag_SC_1min_interp',$               ;RTN 1 min mag
'spp_fld_ephem_spp_hertn_radial_distance_rs', $
'spp_fld_ephem_spp_hertn_radial_velocity', $
'spp_fld_ephem_spp_hertn_position',$
'spp_fld_ephem_spp_hertn_velocity',$
'spp_fld_ephem_spp_hertn_sc_x_vector',$
'spp_fld_ephem_spp_hertn_sc_y_vector',$
'spp_fld_ephem_spp_hertn_sc_z_vector',$
'spp_fld_ephem_eclipj2000_position', $
'spp_fld_ephem_eclipj2000_velocity',$
'psp_fld_l2_quality_flags']


;FIELDS data interpolated to 4 Hz
fields_4hzvars = ['psp_fld_l2_dfb_wf_scm_hg_sc_4Hz',$
'psp_fld_l2_dfb_wf_V1dc_4Hz',$
'psp_fld_l2_dfb_wf_V2dc_4Hz',$
'psp_fld_l2_dfb_wf_V3dc_4Hz',$
'psp_fld_l2_dfb_wf_V4dc_4Hz',$
'psp_fld_l2_dfb_wf_V5dc_4Hz',$
'psp_fld_l2_dfb_wf_V1dc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_V2dc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_V3dc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_V4dc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_V5dc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_scm_hg_sc_4Hz',$
'psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_dVdc_sc_4Hz',$
'psp_fld_l2_dfb_wf_dVdc_sensor_4Hz',$
'psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_scm_lg_sc_4Hz',$
'psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz',$
'psp_fld_l2_mag_RTN_4_Sa_per_Cyc',$       ;~4 Hz RTN mag
'psp_fld_l2_mag_RTN_4_Sa_per_Cyc_deltaT',$
'psp_fld_l2_mag_SC_4_Sa_per_Cyc',$       ;~4 Hz RTN mag
'psp_fld_l2_mag_SC_4_Sa_per_Cyc_deltaT']



;************************************
;***NEED TO CHECK THIS PART FOR EACH PERIGEE PASS
;*************************************

;Interpolate the SWEAP data to the 4 Hz time base of some of the FIELDS data
;varinterp = 'psp_fld_l2_dfb_wf_V1dc_4Hz'

;****NOTE: for the Jan-Feb 2020 perigee pass the SWEAP data goes on much longer
;than the FIELDS data. Therefore I can't interpolate to the FIELDS data.
;Instead, create an artificial 4Hz time base.
t0 = '2020-01-16'
t1 = '2020-02-13'
ndays = floor((time_double(t1) - time_double(t0))/86400)

nelem_day = 86400.*4.
tbase = dindgen(nelem_day*ndays)/4. + time_double(t0)

store_data,'timebase',tbase,tbase
varinterp = 'timebase'

;************************************
;*************************************


tinterpol_mxn,'NP_MOMENT',varinterp,newname='NP_MOMENT_4Hz'
tinterpol_mxn,'WP_MOMENT',varinterp,newname='WP_MOMENT_4Hz'
tinterpol_mxn,'VP_MOMENT_SC',varinterp,newname='VP_MOMENT_SC_4Hz'
tinterpol_mxn,'VP_MOMENT_RTN',varinterp,newname='VP_MOMENT_RTN_4Hz'
tinterpol_mxn,'SC_POS_HCI',varinterp,newname='SC_POS_HCI_4Hz'
tinterpol_mxn,'SC_VEL_HCI',varinterp,newname='SC_VEL_HCI_4Hz'
tinterpol_mxn,'CARR_LATITUDE',varinterp,newname='CARR_LATITUDE_4Hz'
tinterpol_mxn,'CARR_LONGITUDE',varinterp,newname='CARR_LONGITUDE_4Hz'

store_data,['NP_MOMENT','WP_MOMENT','VP_MOMENT_SC','VP_MOMENT_RTN','SC_POS_HCI','SC_VEL_HCI','CARR_LATITUDE','CARR_LONGITUDE'],/del


;Now because I'm not using the FIELDS data to interpolate to 4Hz, I need to interpolate
;the 4Hz FIELDS data to the timebase 4Hz data so they line up exactly with the SWEAP data.
for i=0,n_elements(fields_4hzvars)-1 do tinterpol_mxn,fields_4hzvars[i],varinterp,/overwrite


copy_data,'psp_fld_l2_mag_RTN_4_Sa_per_Cyc','psp_fld_l2_mag_RTN_4Hz'
copy_data,'psp_fld_l2_mag_RTN_4_Sa_per_Cyc_deltaT','psp_fld_l2_mag_RTN_deltaT_4Hz'
copy_data,'psp_fld_l2_mag_SC_4_Sa_per_Cyc','psp_fld_l2_mag_SC_4Hz'
copy_data,'psp_fld_l2_mag_SC_4_Sa_per_Cyc_deltaT','psp_fld_l2_mag_SC_deltaT_4Hz'

store_data,'*4_Sa_per_Cyc*',/del







tplot,['psp_fld_l2_dfb_wf_V5dc_4Hz','NP_MOMENT_4Hz']



tplot_save,'*',filename='~/Desktop/fields_sweap_perigee_jan2020'
