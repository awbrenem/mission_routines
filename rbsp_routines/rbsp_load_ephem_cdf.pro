;rbsp_load_ephem_cdf.pro

;Loads the definitive RBSP ephemeris data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;e.g. /pub/data/rbsp/rbspa/ephemeris/ect-mag-ephem/cdf/def-1min-t89d/2013/rbsp-a_mag-ephem_def-1min-t89d_20130101_v01.cdf
;Can load multiple days of EFW data


;Returns tplot variables


;***EFW variables returned
;rbsp?_r_gse
;rbsp?_v_gse
;rbsp?_q_uvw2gse
;rbsp?_wsc_gse
;rbsp?_mlt
;rbsp?_mlat
;rbsp?_lshell
;rbsp?_sphase_ssha
;rbsp?_spin_period
;rbsp?_spin_phase


;***ECT variables returned (some of them)
;CD=centered dipole
;ED=eccentric dipole
;Pfn=???
;Pfs=???

;MLT variables: CDMAG_MLT, EDMAG_MLT, Pfn_CD_MLT, Pfn_ED_MLT, Pfs_CD_MLT, Pfs_ED_MLT
;L-shell variables: Lsimple, Lm_eq, Lstar, L
;MLAT variables: EDMAG_MLAT, MlatFromBoverBeq



;ECT "types" to load:
;---def-1min-op77q  (rbsp-a_mag-ephem_def-1min-op77q_20120831_v01.cdf)
;---def-1min-t89d   (rbsp-a_mag-ephem_def-1min-t89d_20130101_v01.cdf)
;---def-1min-t89q   (rbsp-a_mag-ephem_def-1min-t89q_20130101_v01.cdf)
;---def-5min-ts04d  (rbsp-a_mag-ephem_def-5min-ts04d_20130101_v01.cdf)



;---------
;e.g. loading EFW spice ephem
;timespan,'2012-11-01',2,/days
;sc = 'a'
;src = 'efw'
;rbsp_load_ephem_cdf,sc
;---------
;e.g. loading ECT definitive mag ephem
;timespan,'2012-11-01'
;sc = 'a'
;type = 'def-1min-t89d'
;src = 'ect'
;rbsp_load_ephem_cdf,sc,type,source='ect'



;pro rbsp_load_ephem_cdf,sc,type,$
;  paths=paths,source=src

timespan,'2012-11-01',2,/days
sc = 'a'
src = 'efw'



  if ~keyword_set(src) then src = 'efw'
  
  if src eq 'ect' then sc2 = 'rbsp-'+sc else sc2 = 'rbsp'+sc
  sc3 = 'rbsp'+sc

  tr = timerange()


  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'


  ndays_load = floor((tr[1]-tr[0])/86400)



  for i=0,ndays_load -1 do begin


    dtst = time_string(tr[0]+i*86400,tformat='YYYY-MM-DD')
    yr = strmid(dtst,0,4)
    mm = strmid(dtst,5,2)
    dd = strmid(dtst,8,2)
  
  
    if src eq 'ect' then begin
      remote_path = url +sc3+'/ephemeris/ect-mag-ephem/cdf/'+type+'/'+yr+'/'
      remote_file = sc2+'_mag-ephem_'+type+'_'+yr+mm+dd+'_v??.cdf'
      local_path = '/Users/abrenema/data/rbsp/'+sc3+'/ephemeris/ect-mag-ephem/cdf/'+type+'/'+yr+'/'
    endif
    if src eq 'efw' then begin
      remote_path = url +sc3+'/ephemeris/efw-ephem/'+yr+'/'
      remote_file = sc2+'_spice_products_'+yr+mm+dd+'_v??.cdf'
      local_path = '/Users/abrenema/data/rbsp/'+sc3+'/ephemeris/efw-ephem/'+yr+'/'
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


          get_data,tnames[j],data=tntmp
          get_data,tnames[j]+'_fin',data=tnfin


          sz = size(tntmp.y,/n_dimensions)
          if sz eq 1 or sz eq 2 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]}
          if sz eq 3 then store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y],v:tnfin.v}


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
