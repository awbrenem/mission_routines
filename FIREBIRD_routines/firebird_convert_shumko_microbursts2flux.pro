;--Calibrate to flux and time-correct Mike Shumko's microbursts. 

;For a single date, turn the microburst values loaded from load_firebird_microburst_list.pro (which are signal/background ratios for each channel from column detector) into counts using:
;   counts = sig*sqrt(A+1) + A,
;where A is the running average.
;Once in counts, these can be turned into fluxes using calibrations from
;  firebird_get_calibration_counts2flux.pro.


;Variables: 
;   fb -> which FIREBIRD flight unit ('3' or '4') 
;   ub_list -> microburst list from firebird_load_shumko_microburst_list.pro 
;   tname -> name of tplot variable containing the smoothed background counts for energy channel of interest;
;     For example, the lowest energy channel counts from firebird_load_data.pro  
;   e_channel -> integer energy channel of collimated detector to use. Exact energy values change with time (see Arlo Johnson 2020 paper), but are generally: 
;     0: 220-283
;     1: 283-383
;     2: 383-520
;     3: 520-720
;     4: 720-985
;     5: >985



;Example:
; timespan,'2019-10-04'
; fb = '4'  ;FIREBIRD flight unit 4
; firebird_load_data,fb
; smooth_time = 3.  ;sec to smooth data over 
;
; e_channel = 0.   ;250 keV channel (lowest FB energy bin) 
; split_vec,'fu4_fb_col_hires_counts'
; rbsp_detrend,'fu4_fb_col_hires_counts_0',smooth_time
; tname = 'fu4_fb_col_hires_counts_0_smoothed'
;
;
; ub_list = firebird_load_shumko_microburst_list(fb,filename='FU4_microbursts_bw=0.75sec.csv')
; ub_flux = firebird_convert_shumko_microbursts2flux(fb, ub_list, tname, e_channel)


function firebird_convert_shumko_microbursts2flux, fb, ub_list, tname, tname_timecorrection, e_channel 



  sctmp = fb


  ;Extract the date of interest
  tr = time_string(timerange())
  date_curr = strmid(tr[0],0,10)




  ;Extract date string from date-time string
  datetmp = strarr(n_elements(ub_list.time))
  for i=0, n_elements(ub_list.time)-1 do datetmp[i] = strmid(ub_list.time[i],0,10)



  ;Get conversion factor for counts to flux  (flux = counts/cadence/energy_width/geometric_factor)
  cal = firebird_get_calibration_counts2flux(date_curr, sctmp)
  flux_conv = 1/((cal.cadence/1000.) * (cal.energy_range_collimated[e_channel,1] - cal.energy_range_collimated[e_channel,0]) * cal.g_factor_collimated[e_channel])





  case e_channel of
     0: sig_vals = ub_list.sig_ch1
     1: sig_vals = ub_list.sig_ch2
     2: sig_vals = ub_list.sig_ch3
     3: sig_vals = ub_list.sig_ch4
     4: sig_vals = ub_list.sig_ch5
     5: sig_vals = ub_list.sig_ch6
  endcase


  good = where(datetmp eq date_curr)
  ub_counts = fltarr(n_elements(good))
  ub_flux = fltarr(n_elements(good))
  ub_times = dblarr(n_elements(good))


  if good[0] ne -1 then begin
    for i=0, n_elements(good)-1 do begin
      A = tsample(tname,time_double(ub_list.time[good[i]]), time=tm)
      tcorr = tsample(tname_timecorrection,time_double(ub_list.time[good[i]]))
      ub_counts[i] = sig_vals[good[i]]*sqrt(A+1) + A
      ub_flux[i] = flux_conv * ub_counts[i]
      ub_times[i] = tm + tcorr
    endfor



    return,{ub_times:ub_times, ub_flux:ub_flux} 
  endif else return, {ub_times:!values.f_nan, ub_flux:!values.f_nan} 


end 
