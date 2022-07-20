;rbsp_load_ephem_cdf.pro

;Loads the definitive RBSP ephemeris data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;e.g. /pub/data/rbsp/rbspa/ephemeris/ect-mag-ephem/cdf/def-1min-t89d/2013/rbsp-a_mag-ephem_def-1min-t89d_20130101_v01.cdf

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
;timespan,'2012-11-01'
;sc = 'a'
;source = 'efw'
;---------
;e.g. loading ECT definitive mag ephem
;timespan,'2012-11-01'
;sc = 'a'
;type = 'def-1min-t89d'
;source = 'ect'


pro rbsp_load_ephem_cdf,sc,type,$
  paths=paths,source=src


  if ~keyword_set(src) then src = 'efw'
  
  if src eq 'ect' then sc2 = 'rbsp-'+sc else sc2 = 'rbsp'+sc
  sc3 = 'rbsp'+sc

  tr = timerange()


  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'



  dtst = time_string(tr[0],tformat='YYYY-MM-DD')
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

  cdf2tplot,file


  paths = {remote_path:remote_path,remote_file:remote_file,local_path:local_path}

  

end
