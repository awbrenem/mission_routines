;Load final RBSP EMFISIS data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;Returns tplot variables
;Can load multiple days of data


;-------------------------
;L2 data "types" to load:
;(NOTE: for L2, type variable needs to combine these (e.g. type='hfr/spectra-burst')
;-------------------------

;hfr
;   spectra-burst
;   spectra-merged (NOT CURRENTLY WORKING - cdf2tplot CAN'T TURN THESE INTO TPLOT VARIABLES)
;   spectra   (WORKS)
;   waveform  (NOT CURRENTLY WORKING - weird data format)
;magnetometer
;   uvw (NOT CURRENTLY WORKING - LOCAL FILENAME IS SLIGHTLY INCORRECT)
;wfr
;   spectral-matrix-diagonal-merged (NOT CURRENTLY WORKING)
;   spectral-matrix-diagonal (NOT CURRENTLY WORKING)
;   spectral-matrix (NOT CURRENTLY WORKING)
;   waveform-continuous-burst (NOT CURRENTLY WORKING)
;   waveform  (NOT CURRENTLY WORKING - weird data format)

;-------------------------
;L3 data "types" to load 
;NOTE: since there's only "magnetometer" data, type='1sec', '4sec', or 'hires'
;with subcategories gei, geo, gse, gsm, sm 
;e.g. type = '4sec/sm'
;-------------------------

;magnetometer
;   1sec
;   4sec
;   hires
      ;gei
      ;geo
      ;gse
      ;gsm
      ;sm



;---------
;e.g. load L2 data
;timespan,'2012-11-01',2,/days
;sc = 'a'
;lvl = 'l2'
;type = 'hfr/spectra'

;---------
;e.g. load L3 data
;timespan,'2012-11-01',2,/days
;sc = 'b'
;lvl = 'l3'
;type = '4sec/geo'

;---------
;Example file names:
;rbsp-a_hfr-waveform_emfisis-l2_20120901_v1.2.9.cdf
;rbsp-a_hfr-spectra-burst_emfisis-l2_20120901_v1.2.6.cdf
;rbsp-a_housekeeping_emfisis-l2_20120104_v1.2.1.cdf  
;rbsp-a_magnetometer_uvw_emfisis-l2_20120830_v1.7.1.cdf
;rbsp-a_magnetometer_hires-gei_emfisis-l3_20120830_v1.7.1.cdf
;-------------------------------------------------------------



pro rbsp_load_emfisis_cdf,sc,lvl,type,$
  paths=paths




  sc2 = 'rbsp-'+sc
  scpath = strmid(sc2,0,4)+strmid(sc2,5,1)


  tr = timerange()
  ndays_load = floor((tr[1]-tr[0])/86400)


  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'
  


  for i=0,ndays_load -1 do begin


    dtst = time_string(tr[0]+i*86400,tformat='YYYY-MM-DD')
    yr = strmid(dtst,0,4)
    mm = strmid(dtst,5,2)
    dd = strmid(dtst,8,2)
    
   
    typefn = strsplit(type,'/',/extract)
  
  
    
    if lvl eq 'l2' then begin  
      remote_path = url +scpath+'/'+lvl+'/emfisis/'+type+'/'+yr+'/'
      remote_file = sc2+'_'+typefn[0]+'-'+typefn[1]+'_emfisis-l2_'+yr+mm+dd+'_v?.?.?.cdf
      local_path = '/Users/abrenema/data/rbsp/'+scpath+'/'+lvl+'/emfisis/'+type+'/'+yr+'/'
    endif
    if lvl eq 'l3' then begin
      remote_path = url +scpath+'/'+lvl+'/emfisis/magnetometer/'+type+'/'+yr+'/'
      remote_file = sc2 + '_magnetometer_'+typefn[0]+'-'+typefn[1]+'_emfisis-l3_'+yr+mm+dd+'_v?.?.?.cdf'
      local_path = '/Users/abrenema/data/rbsp/'+scpath+'/'+lvl+'/emfisis/magnetometer/'+type+'/'+yr+'/'
    endif
  
  
  
    ;Download file and turn into tplot variables
    file = spd_download(remote_file=remote_file,remote_path=remote_path,$
            local_path=local_path,/last_version)

    cdf2tplot,file,varnames=tnames


    ;----------------------------------------------------------------------------
    ;Some of the EMFISIS data has a format that's not compatible with tplot. Fix here
    ;----------------------------------------------------------------------------

    if type eq 'hfr/spectra' then begin
      get_data,'HFR_Spectra',data=dd,dlim=dlim,lim=lim
      store_data,'HFR_Spectra',dd.x,dd.y,reform(dd.v),dlim=dlim,lim=lim
    endif


   
    if ndays_load gt 1 then begin


      ;Rename
      if i eq 0 then for j=0,n_elements(tnames)-1 do copy_data,tnames[j],tnames[j]+'_fin'


      ;If both renamed and original, then combine
      if i gt 0 then begin
        for j=0, n_elements(tnames)-1 do begin


          sz = size(tnames[j],/n_dimensions)


          get_data,tnames[j],data=tntmp
          get_data,tnames[j]+'_fin',data=tnfin

          sz = size(tntmp.y,/n_dimensions)


          ;***This is clunkly but works for l3 data and l2 hfr/spectra, which are currently the only two types I can load
          if lvl eq 'l2' then begin
            if sz eq 1 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]},dlim=dlim,lim=lim
            if sz eq 2 or sz eq 3 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y],v:tnfin.v},dlim=dlim,lim=lim
          endif
          if lvl eq 'l3' then begin
            if sz eq 1 or sz eq 2 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]},dlim=dlim,lim=lim
            if sz eq 3 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y],v:tnfin.v},dlim=dlim,lim=lim
          endif


        endfor  ;for each tplot variable
      endif  ;i>0
      store_data,tnames,/del
    endif ;if more than 1 day to load
  endfor  ;for each day to load




  ;Final rename of variables
  if ndays_load gt 1 then begin
    for j=0, n_elements(tnames)-1 do copy_data,tnames[j]+'_fin',tnames[j]
    store_data,tnames+'_fin',/del
  endif



  paths = {remote_path:remote_path,remote_file:remote_file,local_path:local_path}



end
