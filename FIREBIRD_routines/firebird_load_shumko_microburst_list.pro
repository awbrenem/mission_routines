;Load the list that Mike Shumko's FIREBIRD microburst id code ("microburst_detection" package)
;spits out. These are individual microbursts pulled from the FIREBIRD flux data. 
;These microbursts have been time-corrected and flux corrected (see test_firebird_time_correction_flux_calibration.pro)

;Microbursts with the following flags have been removed:
;'time_gap' 
;'saturated'
;'n_zeros'


;NOTES:
;(1) Mike Shumko's code returns microburst counts/sec, which I calibrate here into flux using firebird_get_calibration_counts2flux.pro.
;(2) I've verified from Mike that he ONLY uses the collimated detector. 
;(3) Mike (email Jan 19, 2022) included three variables to test for "bad" microburst ids. 
;n_zeros: represents the number of time stamps with 0 counts around the microburst (within 5 seconds). 
; This is a good indicator of severe saturation---look at plots where n_zeros > 5 or 10---so I may even trust it more than the dropout flag.
;saturated: 0 or 1 
;time_gap: 0 or 1




;fb = '3' or '4'  corresponding to FIREBIRD FU3 or FU4

;----------------------------------------------------------------------------------------------------


function firebird_load_shumko_microburst_list,fb,filename=fntmp


    ;Set the max allowed value of the flag variable "n_zeros".
    max_n_zeros = 1.

    paths = get_project_paths()


    if not keyword_set(fntmp) then fntmp = 'FU'+fb+'_microburst_catalog_02.csv'



    ft = [7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
    fn = ['time','lat','lon','alt','mcilwainL','MLT','kp','counts_s_0','counts_s_1','counts_s_2','counts_s_3','counts_s_4','counts_s_5','sig_0','sig_1','sig_2','sig_3','sig_4','sig_5','time_gap','saturated','n_zeros']
    floc = [0, 27, 46, 66, 84,102,120,125,143,162,180,199,203,207,225,244,263,283,287,308,310,312]
    fg = indgen(22)

    template = {version:1.,$
      datastart:1L,$
      delimiter:44B,$
      missingvalue:!values.f_nan,$
      commentsymbol:'',$
      fieldcount:22L,$
      fieldtypes:ft,$
      fieldnames:fn,$
      fieldlocations:floc,$
      fieldgroups:fg}



    vals = read_ascii(paths.SHUMKO_MICROBURST_DETECTION+fntmp,template=template)


    ;----------------------------------------------------------------------------------------------------
    ;Remove all "bad" microbursts. Note that some of the remaining ones may also not be microbursts, but
    ;the bad ones are almost always bad. 
    ;----------------------------------------------------------------------------------------------------


    bad = where((vals.time_gap) or (vals.saturated) or (vals.n_zeros gt max_n_zeros))
    vals.time[bad] = !values.f_nan
    good = where(finite(vals.time) ne 0.)    
    times = time_double(vals.time[good])



    ;--------------------------
    ;Time correct the data 
    ;--------------------------

    firebird_load_time_corrections


    get_data,'FU'+fb+'_tcorr',data=tcorr

    for i=0,n_elements(times)-1 do begin
       goo = where(times[i] le tcorr.x)
       if goo[0] ne -1 then times[i] += tcorr.y[goo[0]]
    endfor


    ;Corrected times
    times = time_string(times,prec=3)





    ;day occurrence of each microburst
    dates = time_string(times,tformat='YYYY-MM-DD')


    ;------------------------------
    ;Convert counts to flux (only differential channels)
    flux = fltarr(n_elements(vals.counts_s_0[good]),5)


    cadence = 1. ;sec  (Mike's counts are actually counts/sec)


    for i=0.,n_elements(vals.counts_s_0[good])-1 do begin
      x = firebird_get_calibration_counts2flux(dates[i],fb)
      energy_width = x.energy_range_collimated[*,1] - x.energy_range_collimated[*,0]


      if is_struct(x) then begin 
        ;To calibrate from counts to flux (only for differential channels):
        ;(see header of firebird_get_calibration_counts2flux.pro)
        flux[i,0] = vals.counts_s_0[good[i]]/cadence/energy_width[0]/x.g_factor_collimated[0]
        flux[i,1] = vals.counts_s_1[good[i]]/cadence/energy_width[1]/x.g_factor_collimated[1]
        flux[i,2] = vals.counts_s_2[good[i]]/cadence/energy_width[2]/x.g_factor_collimated[2]
        flux[i,3] = vals.counts_s_3[good[i]]/cadence/energy_width[3]/x.g_factor_collimated[3]
        flux[i,4] = vals.counts_s_4[good[i]]/cadence/energy_width[4]/x.g_factor_collimated[4]
      endif else begin
        flux[i,*] = !values.f_nan
      endelse
      
    endfor





    vals_fin = {TIME:times,LAT:vals.lat[good],LON:vals.lon[good],ALT:vals.alt[good],MCILWAINL:vals.mcilwainl[good],MLT:vals.mlt[good],KP:vals.kp[good],$
    flux_0:flux[*,0],flux_1:flux[*,1],flux_2:flux[*,2],flux_3:flux[*,3],flux_4:flux[*,4],$
    COUNTS_S_0:vals.counts_s_0[good],COUNTS_S_1:vals.counts_s_1[good],COUNTS_S_2:vals.counts_s_2[good],COUNTS_S_3:vals.counts_s_3[good],$
    COUNTS_S_4:vals.counts_s_4[good],COUNTS_S_5:vals.counts_s_5[good],$
    SIG_0:vals.sig_0[good],SIG_1:vals.sig_1[good],SIG_2:vals.sig_2[good],SIG_3:vals.sig_3[good],SIG_4:vals.sig_4[good],SIG_5:vals.sig_5[good]}


    return,vals_fin

end