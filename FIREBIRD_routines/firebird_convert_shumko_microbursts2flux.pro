;For a single date, turn the microburst values loaded from load_firebird_microburst_list.pro (which are signal/background ratios for each channel) into counts using:
;   counts = sig*sqrt(A+1) + A,
;where A is the running average.
;Once in counts, these can be turned into fluxes using calibrations from
;  firebird_get_calibration_counts2flux.pro.

;Example:
; timespan,'2019-10-04'
; fb = '4'  ;FIREBIRD flight unit 4
; firebird_load_data,fb
; dettime = 3.  ;sec

; split_vec,'fu4_fb_col_hires_counts'
; rbsp_detrend,'fu4_fb_col_hires_counts_0',dettime
; tname = 'fu4_fb_col_hires_counts_0_smoothed'


; ub_list = firebird_load_shumko_microburst_list(fb,filename='FU4_microbursts_bw=0.75sec.csv')


; firebird_convert_shumko_microbursts2flux, ub_list, tname


pro firebird_convert_shumko_microbursts2flux, ub_list, tname 




  ;Extract the date of interest
  tr = time_string(timerange())
  date_curr = strmid(tr[0],0,10)

  datetmp = strarr(n_elements(ub_list.time))
  for i=0, n_elements(ub_list.time)-1 do datetmp[i] = strmid(ub_list.time[i],0,10)


  sctmp = sc
  cal = firebird_get_calibration_counts2flux(datetime,sctmp)
  ;chn = strtrim(cal.CHANNEL_USED_FOR_SURVEY_CALIBRATION - 1,2)

  flux_conv = 1/((cal.cadence/1000.) * (cal.energy_range_collimated[0,1] - cal.energy_range_collimated[0,0]) * cal.g_factor_collimated[0])
  ; flux = counts/cadence/energy_width/geometric_factor



  good = where(datetmp eq date_curr)
  ub_counts = fltarr(n_elements(good))
  ub_flux = fltarr(n_elements(good))

  if good[0] ne -1 then begin
    for i=0, n_elements(good)-1 do begin
      sig = ub_list.sig_ch1[good[i]]
      A = tsample(,time_double(ub_list.time[good[i]]), time=tm)

      ub_counts[i] = sig*sqrt(A+1) + A
      ub_flux[i] = flux_conv * ub_counts[i]

      ;firebird_get_calibration_sig2counts,sig,A
    endfor






end 
