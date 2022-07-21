;Load final RBSP EFW data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;Returns tplot variables
;NOTE: can load multiple days of data


;L2 data "types" to load:
;---e-highres-uvw
;---e-spinfit-mgse
;---esvy_despun
;---fbk
;---spec
;---vsvy-highres


;L3 data "types" to load 
;N/A - only a single type


;---------
;e.g. load L2 data
;timespan,'2012-11-01',2,/days
;sc = 'a'
;lvl = 'l2'
;type = 'fbk'

;---------
;e.g. load L3 data
;timespan,'2012-11-01',2,/days
;sc = 'a'
;lvl = 'l3'
;type = N/A



pro rbsp_load_efw_cdf,sc,lvl,type,$
  paths=paths



  ;Fix a very common user input error
  if type eq 'vsvy-hires' then type = 'vsvy-highres'
  if type eq 'e-hires-uvw' then type = 'e-highres-uvw'



  sc2 = 'rbsp'+sc

  tr = timerange()

 
  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'
  


  ndays_load = floor((tr[1]-tr[0])/86400)



  for i=0,ndays_load -1 do begin 


    dtst = time_string(tr[0]+i*86400,tformat='YYYY-MM-DD')
    yr = strmid(dtst,0,4)
    mm = strmid(dtst,5,2)
    dd = strmid(dtst,8,2)
  
  
    
    
    if lvl eq 'l2' then begin    
      ;The following need to be adjusted. The folder is "vsvy-highres" but the filename has "vsvy-hires"
      ;Similarly, "e-hires-uvw" --> "e-hires-uvw"
      type2 = type
      if type eq 'vsvy-highres' then type2 = 'vsvy-hires'
      if type eq 'e-highres-uvw' then type2 = 'e-hires-uvw'
      remote_path = url +sc2+'/'+lvl+'/efw/'+type+'/'+yr+'/'
      remote_file = sc2+'_efw-'+lvl+'_'+type2+'_'+yr+mm+dd+'_v??.cdf'
      local_path = '/Users/abrenema/data/rbsp/'+sc2+'/'+lvl+'/efw/'+type+'/'+yr+'/'
    endif
    if lvl eq 'l3' then begin
      remote_path = url +sc2+'/'+lvl+'/efw/'+yr+'/'
      remote_file = sc2+'_efw-'+lvl+'_'+yr+mm+dd+'_v??.cdf'
      local_path = '/Users/abrenema/data/rbsp/'+sc2+'/'+lvl+'/efw/'+yr+'/'
    endif
  
  
  
    ;Download file and turn into tplot variables
    file = spd_download(remote_file=remote_file,remote_path=remote_path,$
            local_path=local_path,/last_version)
   
    cdf2tplot,file,varnames=tnames
  
  

    if ndays_load gt 1 then begin
      
      
      ;Rename 
      if i eq 0 then for j=0,n_elements(tnames)-1 do copy_data,tnames[j],tnames[j]+'_fin'      
      
      
      ;If both renamed and original, then combine
      if i gt 0 then begin
        for j=0, n_elements(tnames)-1 do begin


sz = size(tnames[j],/n_dimensions)

;print,'****' + tnames[j]

          get_data,tnames[j],data=tntmp
          get_data,tnames[j]+'_fin',data=tnfin
          
          sz = size(tntmp.y,/n_dimensions)

;print,'***before size'
;help,tnfin
          if sz eq 1 or sz eq 2 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]}
          if sz eq 3 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y],v:tnfin.v}

;get_data,tnames[j]+'_fin',data=ttt
;print,'***after size'
;help,ttt
;stop
          
        endfor  ;for each tplot variable
      endif  ;i>0 
      store_data,tnames,/del
    endif  ;if more than 1 day to load



  endfor  ;for each day to load


  ;Final rename of variables
  if ndays_load gt 1 then begin
    for j=0, n_elements(tnames)-1 do copy_data,tnames[j]+'_fin',tnames[j] 
    store_data,tnames+'_fin',/del
  endif



  paths = {remote_path:remote_path,remote_file:remote_file,local_path:local_path}

end
