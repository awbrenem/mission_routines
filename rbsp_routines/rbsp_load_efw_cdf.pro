;Load final RBSP EFW data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;Returns tplot variables


;L2 data types to load:
;---e-highres-uvw
;---e-spinfit-mgse
;---e-spinfit-mgse_both_booms
;---esvy_despun
;---fbk
;---spec
;---vsvy-highres


;L3 data types to load 
;N/A - only a single type



;e.g. load L2 data
;timespan,'2012-11-01'
;sc = 'rbspa'
;lvl = 'l2'
;type = 'fbk'

;e.g. load L3 data
;timespan,'2012-11-01'
;sc = 'rbspa'
;lvl = 'l3'
;type = N/A



pro rbsp_load_efw_cdf,sc,lvl,type



  tr = timerange()


  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'
  


  dtst = time_string(tr[0],tformat='YYYY-MM-DD')
  yr = strmid(dtst,0,4)
  mm = strmid(dtst,5,2)
  dd = strmid(dtst,8,2)
  
 
  
  if lvl eq 'l2' then begin  
    remote_path = url +sc+'/'+lvl+'/efw/'+type+'/'+yr+'/'
    remote_file = sc+'_efw-'+lvl+'_'+type+'_'+yr+mm+dd+'_v??.cdf'
    local_path = 'data/rbsp/'+sc+'/'+lvl+'/efw/'+type+'/'+yr+'/'
  endif
  if lvl eq 'l3' then begin
    remote_path = url +sc+'/'+lvl+'/efw/'+yr+'/'
    remote_file = sc+'_efw-'+lvl+'_'+yr+mm+dd+'_v??.cdf'
    local_path = 'data/rbsp/'+sc+'/'+lvl+'/efw/'+yr+'/'
  endif



  ;Download file and turn into tplot variables
  file = spd_download(remote_file=remote_file,remote_path=remote_path,$
          local_path=local_path,/last_version)
  cdf2tplot,file



end
