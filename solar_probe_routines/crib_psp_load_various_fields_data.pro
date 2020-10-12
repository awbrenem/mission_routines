;Cribsheet for loading various PSP L2 data quantities.
;Saves as a tplot variable that's used by crib_psp_combine_sweap_fields_perigee_wygant.pro

rbsp_efw_init

setenv, 'PSP_STAGING_DIR=/Users/aaronbreneman/data/spp/' ; <- replace this
setenv, 'USER=wygant' ; <- replace this
setenv, 'PSP_STAGING_PW=flds@psp' ; <- replace this

;;Perihelion 1 (Nov 5th, 2018)
;t0 = '2018-11-01'
;t1 = '2018-11-13'
;filename = '~/Desktop/psp_fields_vars_nov2018'
;Perihelion 2 (April 4th, 2019)
;t0 = '2019-03-30'
;t1 = '2019-04-13'
;filename = '~/Desktop/psp_fields_vars_apr2019'
;Perihelion 3 (Sept 1st, 2019)
;t0 = '2019-08-26'
;t1 = '2019-09-10'
;filename = '~/Desktop/psp_fields_vars_sept2019'
;Perihelion 4 (Jan 29, 2020)

;t0 = '2020-01-17'
;t1 = '2020-01-19'

t0 = '2020-01-16'
t1 = '2020-02-13'
filename = '~/Desktop/psp_fields_vars_jan2020'



varinterp = 'psp_fld_l2_mag_RTN_4_Sa_per_Cyc'   ;interpolate (non 1-min) data to this variable


ndays = floor((time_double(t1) - time_double(t0))/86400)
timespan,t0,ndays,/days


;-----------------------------------------------------------
;Load L2 FIELDS data
;-----------------------------------------------------------

;type = 'mag_RTN'  ;normal mag data
;type = 'mag_RTN_4_Sa_per_Cyc'   for downsampled data
;type = 'mag_RTN_1min'   for even coarser downsampling
;You could also use 'SC' in place of RTN coordinates.







;spp_fld_load,type = 'rfs_hfr'
;;39 psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2
;;47 psp_fld_l2_rfs_hfr_auto_peaks_ch0_V1V2
;store_data,['psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2_gain','psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2_nsum','psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2_hl','psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4','psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4_gain','psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4_nsum','psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4_hl','psp_fld_l2_rfs_hfr_auto_peaks_ch0_V1V2_gain'],/del
;store_data,['psp_fld_l2_rfs_hfr_auto_peaks_ch0_V1V2_nsum','psp_fld_l2_rfs_hfr_auto_peaks_ch0_V1V2_hl','psp_fld_l2_rfs_hfr_auto_peaks_ch1_V3V4','psp_fld_l2_rfs_hfr_auto_peaks_ch1_V3V4_gain','psp_fld_l2_rfs_hfr_auto_peaks_ch1_V3V4_nsum','psp_fld_l2_rfs_hfr_auto_peaks_ch1_V3V4_hl','psp_fld_l2_rfs_hfr_cross_im_V1V2_V3V4','psp_fld_l2_rfs_hfr_cross_im_V1V2_V3V4_gain','psp_fld_l2_rfs_hfr_cross_im_V1V2_V3V4_nsum'],/del
;store_data,['psp_fld_l2_rfs_hfr_cross_im_V1V2_V3V4_hl','psp_fld_l2_rfs_hfr_cross_re_V1V2_V3V4','psp_fld_l2_rfs_hfr_cross_re_V1V2_V3V4_gain','psp_fld_l2_rfs_hfr_cross_re_V1V2_V3V4_nsum','psp_fld_l2_rfs_hfr_cross_re_V1V2_V3v4_hl','psp_fld_l2_rfs_hfr_coher_V1V2_V3V4','psp_fld_l2_rfs_hfr_coher_V1V2_V3V4_gain'],/del
;store_data,['psp_fld_l2_rfs_hfr_coher_V1V2_V3V4_nsum','psp_fld_l2_rfs_hfr_coher_V1V2_V3V4_hl','psp_fld_l2_rfs_hfr_phase_V1V2_V3V4','psp_fld_l2_rfs_hfr_phase_V1V2_V3V4_gain','psp_fld_l2_rfs_hfr_phase_V1V2_V3V4_nsum','psp_fld_l2_rfs_hfr_phase_V1V2_V3V4_hl','psp_fld_l2_rfs_hfr_averages','psp_fld_l2_rfs_hfr_peaks','psp_fld_l2_rfs_hfr_ch0','psp_fld_l2_rfs_hfr_ch0_string','psp_fld_l2_rfs_hfr_ch1','psp_fld_l2_rfs_hfr_ch1_string','psp_fld_l2_rfs_hfr_nsum','psp_fld_l2_rfs_hfr_gain','psp_fld_l2_rfs_hfr_hl'],/del



spp_fld_load,type = 'mag_RTN_4_Sa_per_Cyc'
get_data,'psp_fld_l2_mag_RTN_4_Sa_per_Cyc',data=d
tdiff = d.x - shift(d.x,1)
tdiff[0] = tdiff[1]
store_data,'psp_fld_l2_mag_RTN_4_Sa_per_Cyc_deltaT',d.x,tdiff    ; & ylim,'mag_RTN_4_Sa_per_Cyc_deltaT',0,80
tplot,['psp_fld_l2_mag_RTN_4_Sa_per_Cyc','psp_fld_l2_mag_RTN_4_Sa_per_Cyc_deltaT']

spp_fld_load,type = 'mag_SC_4_Sa_per_Cyc'
get_data,'psp_fld_l2_mag_SC_4_Sa_per_Cyc',data=d
tdiff = d.x - shift(d.x,1)
tdiff[0] = tdiff[1]
store_data,'psp_fld_l2_mag_SC_4_Sa_per_Cyc_deltaT',d.x,tdiff    ; & ylim,'mag_RTN_4_Sa_per_Cyc_deltaT',0,80
tplot,['psp_fld_l2_mag_SC_4_Sa_per_Cyc','psp_fld_l2_mag_SC_4_Sa_per_Cyc_deltaT']



;The following are high-res files and can be difficult to load. So, break them
;up by loading single day at a time, downsampling data, then moving to next day.

;names = strmid(time_string(tnew),0,10)

tnew = strmid(time_string(time_double(t0) + 86400*indgen(ndays)),0,10)

get_data,varinterp,times,d

for i=0,ndays-1 do begin

  goo = where((times ge time_double(tnew[i])) and (times le time_double(tnew[i])+86400d))
  ttmp = times[goo]

  timespan,time_double(tnew[i]),1,/days

  spp_fld_load,type = 'dfb_wf_dvdc'  ;~300 Samples/sec
  tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sc',ttmp,newname='psp_fld_l2_dfb_wf_dVdc_sc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sensor',ttmp,newname='psp_fld_l2_dfb_wf_dVdc_sensor_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz_' + tnew[i]
  ;del the high res versions
  store_data,['psp_fld_l2_dfb_wf_dVdc_sc','psp_fld_l2_dfb_wf_dVdc_sensor','psp_fld_l2_dfb_wf_dVdc_sample_rate'],/del



  spp_fld_load,type = 'dfb_wf_scm' ;~300 Samples/sec mag data
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sensor',ttmp,newname='psp_fld_l2_dfb_wf_scm_hg_sensor_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sc',ttmp,newname='psp_fld_l2_dfb_wf_scm_hg_sc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_lg_sensor',ttmp,newname='psp_fld_l2_dfb_wf_scm_lg_sensor_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_lg_sc',ttmp,newname='psp_fld_l2_dfb_wf_scm_lg_sc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_lg_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz_' + tnew[i]
  store_data,['psp_fld_l2_dfb_wf_scm_hg_sensor','psp_fld_l2_dfb_wf_scm_hg_sc','psp_fld_l2_dfb_wf_scm_hg_sample_rate','psp_fld_l2_dfb_wf_scm_lg_sensor','psp_fld_l2_dfb_wf_scm_lg_sc','psp_fld_l2_dfb_wf_scm_lg_sample_rate'],/del




  spp_fld_load,type = 'dfb_wf_vdc' ;up to ~150 Samples/sec electric field data
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V1dc',ttmp,newname='psp_fld_l2_dfb_wf_V1dc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V2dc',ttmp,newname='psp_fld_l2_dfb_wf_V2dc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V3dc',ttmp,newname='psp_fld_l2_dfb_wf_V3dc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V4dc',ttmp,newname='psp_fld_l2_dfb_wf_V4dc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V5dc',ttmp,newname='psp_fld_l2_dfb_wf_V5dc_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V1dc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_V1dc_sample_rate_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V2dc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_V2dc_sample_rate_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V3dc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_V3dc_sample_rate_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V4dc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_V4dc_sample_rate_4Hz_' + tnew[i]
  tinterpol_mxn,'psp_fld_l2_dfb_wf_V5dc_sample_rate',ttmp,newname='psp_fld_l2_dfb_wf_V5dc_sample_rate_4Hz_' + tnew[i]
  store_data,['psp_fld_l2_dfb_wf_V?dc','psp_fld_l2_dfb_wf_V?dc_sample_rate'],/del


endfor




;Combine all of the single-day tplot variables from above.
tfin = [0d]
valfin = [[0.],[2]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz_'+tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d]
valfin = [[0.],[2]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_dVdc_sensor_4Hz_'+tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_dVdc_sensor_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d]
valfin = [[0.],[2]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_dVdc_sc_4Hz_'+tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_dVdc_sc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif
;***********

tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_hg_sensor_4Hz_'+tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y] & $
  if is_struct(d) then freqv = d.v
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_hg_sensor_4Hz',tfin,valfin,freqv,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_hg_sc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y] & $
  if is_struct(d) then freqv = d.v
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_hg_sc_4Hz',tfin,valfin,freqv,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_lg_sensor_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y] & $
  if is_struct(d) then freqv = d.v
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_lg_sensor_4Hz',tfin,valfin,freqv,dlim=dlim,lim=lim
endif


tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_lg_sc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y] & $
  if is_struct(d) then freqv = d.v
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_lg_sc_4Hz',tfin,valfin,freqv,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [[0],[0],[0]]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y] & $
  if is_struct(d) then freqv = d.v
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem,*]
  store_data,'psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz',tfin,valfin,freqv,dlim=dlim,lim=lim
endif

;************************************


tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V1dc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V1dc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V2dc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V2dc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V3dc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V3dc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V4dc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V4dc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V5dc_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V5dc_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V1dc_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V1dc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V2dc_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V2dc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V3dc_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V3dc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V4dc_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V4dc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif

tfin = [0d] & valfin = [0]
for i=0,ndays-1 do begin $
  get_data,'psp_fld_l2_dfb_wf_V5dc_sample_rate_4Hz_' + tnew[i],data=d,dlim=dlim,lim=lim & $
  if is_struct(d) then tfin = [tfin,d.x] & $
  if is_struct(d) then valfin = [valfin,d.y]
endfor
if n_elements(tfin) gt 1 then begin
  nelem = n_elements(tfin)-1
  tfin = tfin[1:nelem]
  valfin = valfin[1:nelem]
  store_data,'psp_fld_l2_dfb_wf_V5dc_sample_rate_4Hz',tfin,valfin,dlim=dlim,lim=lim
endif


;ylim,['psp_fld_l2_dfb_wf_V?dc'],-1,1



;Reset the timespan
timespan,t0,ndays,/days



;ephemeris
spp_fld_make_or_retrieve_cdf, 'ephem_eclipj2000', /load
spp_fld_make_or_retrieve_cdf, 'ephem_spp_hertn', /load
tplot, 'spp_fld_ephem_' + [$
'spp_hertn_radial_distance_rs', $
'spp_hertn_radial_velocity', $
'eclipj2000_position', $
'eclipj2000_velocity']



;NOT TOO USEFUL
;spp_fld_load,type = 'f2_100bps' ;HSK-type stuff and FBK data
;spp_fld_load,type = 'mag_VSO' ;Mostly nonexistent
;spp_fld_load,type = 'dfb_ac_bpf'
;spp_fld_load,type = 'dfb_dc_bpf'





;----------------------------------------------------------
;FINAL VARIABLES
;----------------------------------------------------------


;Interpret hires data to 4 sec cadence
;varinterp = 'psp_fld_l2_mag_RTN_4_Sa_per_Cyc'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sc',varinterp,newname='psp_fld_l2_dfb_wf_scm_hg_sc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V1dc',varinterp,newname='psp_fld_l2_dfb_wf_V1dc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V2dc',varinterp,newname='psp_fld_l2_dfb_wf_V2dc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V3dc',varinterp,newname='psp_fld_l2_dfb_wf_V3dc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V4dc',varinterp,newname='psp_fld_l2_dfb_wf_V4dc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V5dc',varinterp,newname='psp_fld_l2_dfb_wf_V5dc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V1dc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_V1dc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V2dc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_V2dc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V3dc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_V3dc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V4dc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_V4dc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_V5dc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_V5dc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sc',varinterp,newname='psp_fld_l2_dfb_wf_scm_hg_sc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_hg_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sc',varinterp,newname='psp_fld_l2_dfb_wf_dVdc_sc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sensor',varinterp,newname='psp_fld_l2_dfb_wf_dVdc_sensor_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_dVdc_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_lg_sc',varinterp,newname='psp_fld_l2_dfb_wf_scm_lg_sc_4Hz'
;tinterpol_mxn,'psp_fld_l2_dfb_wf_scm_lg_sample_rate',varinterp,newname='psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz'


spp_fld_load,type = 'mag_RTN_1min'
;TPLOT(356):   38 psp_fld_l2_mag_RTN_1min
;TPLOT(356):   37 psp_fld_l2_quality_flags
spp_fld_load,type = 'mag_SC_1min'


;Interpolate 1min data to flag cadence (the two are extremely close but not exact)
varinterp2 = 'psp_fld_l2_quality_flags'
tinterpol_mxn,'psp_fld_l2_mag_RTN_1min',varinterp2,newname='psp_fld_l2_mag_RTN_1min_interp'
tinterpol_mxn,'psp_fld_l2_mag_RTN_1min_deltaT',varinterp2,newname='psp_fld_l2_mag_RTN_1min_deltaT_interp'
tinterpol_mxn,'psp_fld_l2_mag_SC_1min',varinterp2,newname='psp_fld_l2_mag_SC_1min_interp'


get_data,'psp_fld_l2_mag_RTN_1min_interp',data=d
tdiff = d.x - shift(d.x,1)
tdiff[0] = tdiff[1]
store_data,'psp_fld_l2_mag_RTN_1min_deltaT_interp',d.x,tdiff & ylim,'psp_fld_l2_mag_RTN_1min_deltaT_interp',0,80
tplot,['psp_fld_l2_mag_RTN_1min_deltaT_interp','psp_fld_l2_mag_RTN_1min_interp']





vars = ['psp_fld_l2_mag_RTN_1min_interp',$               ;RTN 1 min mag
'psp_fld_l2_mag_RTN_1min_deltaT_interp',$        ;RTN 1 min mag dT
'psp_fld_l2_mag_SC_1min_interp',$               ;RTN 1 min mag
'psp_fld_l2_dfb_wf_V?dc_4Hz',$                ;V1-V5  (~3 to 150 Samples/sec)
'psp_fld_l2_dfb_wf_V?dc_sample_rate_4Hz',$
'psp_fld_l2_mag_RTN_4_Sa_per_Cyc',$       ;~4 Hz RTN mag
'psp_fld_l2_mag_RTN_4_Sa_per_Cyc_deltaT',$
'psp_fld_l2_mag_SC_4_Sa_per_Cyc',$       ;~4 Hz RTN mag
'psp_fld_l2_mag_SC_4_Sa_per_Cyc_deltaT',$
'psp_fld_l2_dfb_wf_dVdc_sc_4Hz',$   ;dVX and dVY (~70 to 300 Samples/sec)
'psp_fld_l2_dfb_wf_dVdc_sensor_4Hz',$   ;volts dv12, dv34 (~70 to 300 Samples/sec)
'psp_fld_l2_dfb_wf_dVdc_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_scm_hg_sc_4Hz',$    ;~70 to 300 Samples/sec
'psp_fld_l2_dfb_wf_scm_hg_sample_rate_4Hz',$
'psp_fld_l2_dfb_wf_scm_lg_sc_4Hz',$     ;~150 Samples/sec
'psp_fld_l2_dfb_wf_scm_lg_sample_rate_4Hz',$
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



stop


tplot_save,vars,filename=filename
end
