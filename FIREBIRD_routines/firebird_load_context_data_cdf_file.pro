;Load FIREBIRD context (survey) CDF files (can load multiple days) that Aaron created at
;http://rbsp.space.umn.edu/firebird/FU?/YYYY/
;from the FIREBIRD context files (each of these spans an entire campaign, which is why I broke them up into daily files)
;
;NOTE: This code both time-corrects the data and applies the proper calibration of counts to flux 
;(For testing see firebird_test_time_correction_flux_calibration.pro)
;
;
;
;e.g. FU3_context_20150204_v01.cdf

;These files are created with
;firebird_create_cdf_files


;These have the tplot variables:
;1 D0  (counts/6sec, not flux) - from integral detector
;2 D1  (counts/6sec, not flux) - from a single differential detector channel
;3 Alt ("Altitude of spacecraft as calculated by STK")
;4 Count_Time_Correction ("Time correction for count data expressed as (Ground Time - Spacecraft Time). Should be added to reported time stamp. Applies only to count and flux data, ephemeris does not need to be corrected.")
;5 Flag ("Data flags from the CRC; 1 is good data, 0 is potentially bad data. False negatives (0's) are possible, false positives (1's) are not.")
;6 Lat ("Latitude of spacecraft as calculated by STK")
;7 Lon ("Longitude of spacecraft as calculated by STK")
;8 Loss_cone_type (0=open; 1=trapped; 2=BLC)
;9 MLT
;10 McIlwainL ("Mcllwain L shell parameter for particles mirroring at the spacecraft using the T89 external model and IGRF inernal model")
;11 kp


;Note on counts vs flux:
;Calibration is done by calling firebird_get_calibration_counts2flux.pro 




;Usage
;timespan,'2016-06-11',4,/days
;firebird_load_context_data_cdf_file,'4'

;options are '3' or '4'



pro firebird_load_context_data_cdf_file,cubesat,$
  file_fail=file_fail   ;indicates if no file was returned



  tr = timerange()
  ndays_load = floor((tr[1]-tr[0])/86400)

  

  for i=0,ndays_load -1 do begin

    date = time_string(tr[0]+i*86400,/date_only,tformat='YYYYMMDD')
  
  
    year = strmid(date,0,4)
    mn = strmid(date,4,2)
    dy = strmid(date,6,2)
  
  
    ;file to load
    fn = 'FU' + cubesat + '_context_'+year+mn+dy+'_v01.cdf'
    ;Grab local path to save data
    homedir = (file_search('~',/expand_tilde))[0]+'/'
    folder = homedir+'data/firebird/FU'+cubesat+'/' + year + '/'
    url = 'http://rbsp.space.umn.edu/firebird/'     ;FU?/YYYY/
    path = 'FU' + cubesat+'/'+year+'/'+fn
    file_loaded = spd_download(remote_file=url+path,$
                  local_path=folder,$
                  local_file=folder+fn,$
                  /last_version)
  
  
  
    cdf2tplot,file_loaded,tplotnames=tntmp 
  

;    if tn[0] eq '' then begin 
;      print,'-----------------------------------------------------------'
;      print,'***********************************************************'
;      print,'NO CONTEXT CDF FILE EXISTS FOR THIS DAY....RETURNING'
;      print,'***********************************************************'
;      print,'-----------------------------------------------------------'
;      file_fail = 1
;      return
;    endif
  
  
  
    if tntmp[0] ne '' then begin 
  
  
      ;--------------------------------------------------------------
      ;Change counts to flux (see header for firebird_get_calibration_counts2flux.pro)
      ;--------------------------------------------------------------
    
    
      cal = firebird_get_calibration_counts2flux(date,cubesat)
    
    
      ;array reference for energy channels and geometric factor
      index = cal.channel_used_for_survey_calibration-1
    
      ;Only campaign 1 doesn't have a good calibration for survey channel (b/c mostly non-functioning surface detectors used)
      if cal.channel_type_for_survey_data eq 'collimated' then gfactor = cal.g_factor_collimated[index]
      if cal.channel_type_for_survey_data eq 'surface'    then gfactor = !values.f_nan
    
      if cal.channel_type_for_survey_data eq 'collimated' then channel_energies = cal.energy_range_collimated[index,*]
      if cal.channel_type_for_survey_data eq 'surface'    then channel_energies = !values.f_nan
    
    
      if finite(channel_energies[0]) then dE = channel_energies[1] - channel_energies[0] else dE = !values.f_nan   ;keV
    
    
    
      get_data,'D1',data=d1  ;Only differential channel (counts/6sec)
    
      ;Divide out the integration time for the counts channel
      cadence = 6.  ;sec   (counts are actually counts/6sec)
    
      flux = d1.y/(dE*gfactor*cadence)   ;#/s-cm2-sr-keV  
    
    
    
    
    
      ;Apply time correction and rename channels
      get_data,'Count_Time_Correction',data=tc 
    
    
      ;time-corrected times
      time = d1.x+tc.y
    
    
    
    
      store_data,'flux_context_FU'+cubesat,time,flux
      options,'flux_context_FU'+cubesat,'ytitle','Context flux:!CFrom '+ cal.channel_type_for_survey_data +' channel' + strtrim(cal.channel_used_for_survey_calibration,2)
      options,'flux_context_FU'+cubesat,'ysubtitle','[#/keV-sr-cm2-s]'
    
      store_data,'counts_context_FU'+cubesat,time,d1.y
      options,'counts_context_FU'+cubesat,'ytitle','Context counts!CFrom '+ cal.channel_type_for_survey_data +' channel' + strtrim(cal.channel_used_for_survey_calibration,2)
      options,'counts_context_FU'+cubesat,'ysubtitle','[counts/6sec]'
    
      get_data,'D0',data=d0  ;Only integral channel (counts/6sec)
      store_data,'counts_context_integral_channel_FU'+cubesat,d0.x+tc.y,d0.y
      options,'counts_context_integral_channel_FU'+cubesat,'ytitle','Context counts!Cintegral channel'
      options,'counts_context_integral_channel_FU'+cubesat,'ysubtitle','[counts/6sec]'
    
      store_data,'channel_used_for_survey_calibration'+cubesat,d0.x+tc.y,replicate(cal.channel_used_for_survey_calibration,n_elements(d0.x))
      ylim,'channel_used_for_survey_calibration'+cubesat,0,10
    
    
    
      ;Correct the times on the following:
      get_data,'Flag',data=dd
      store_data,'Flag',time,dd.y
      get_data,'McIlwainL',data=dd
      store_data,'McIlwainL',time,dd.y
      get_data,'Loss_cone_type',data=dd
      store_data,'Loss_cone_type',time,dd.y
    
    
    
      ylim,'flux_context_FU'+cubesat,0.1,1000,1
      ylim,'Flag',0,2 
      ylim,'McIlwainL',0,12
      ylim,'Loss_cone_type',0,3
    
      options,['MLT','McIlwainL','kp','Alt','Flag','Loss_cone_type'],'panel_size',0.5
      options,'flux_context_FU'+cubesat,'psym',-4
    
    ;  tplot,['flux_context_FU'+cubesat,$
    ;        'MLT',$
    ;        'McIlwainL',$
    ;        'Alt',$
    ;        'Flag',$
    ;        'Loss_cone_type']              
    
    
    
      ;final tplot names    
      tnfinvars = ['flux_context_FU'+cubesat,'counts_context_FU'+cubesat,'counts_context_integral_channel_FU'+cubesat,$
        'channel_used_for_survey_calibration'+cubesat,'Alt','Count_Time_Correction','Flag','Lat','Lon',$
        'Loss_cone_type','MLT','McIlwainL','kp']

    
    
      if ndays_load gt 1 then begin
    
    
        ;Rename to temporary variable
        if i eq 0 then for j=0,n_elements(tnfinvars)-1 do copy_data,tnfinvars[j],tnfinvars[j]+'_fin'
    
    
        ;If both renamed and original, then combine
        if i gt 0 then begin
          for j=0, n_elements(tnfinvars)-1 do begin
            get_data,tnfinvars[j],data=tntmp
            get_data,tnfinvars[j]+'_fin',data=tnfin
            store_data,tnfinvars[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]}  
          endfor  ;for each tplot variable

        endif  ;i>0
        store_data,tnfinvars,/del
      endif  ;if more than 1 day to load
    endif  ;for file exists for current day
    
  
  endfor  ;for each day to load
  
  
  
  
  ;Final rename of variables
  if ndays_load gt 1 then begin
    for j=0, n_elements(tnfinvars)-1 do copy_data,tnfinvars[j]+'_fin',tnfinvars[j]
    store_data,tnfinvars+'_fin',/del
  endif
  

 
  
  
  store_data,['D0','D1'],/del

end
