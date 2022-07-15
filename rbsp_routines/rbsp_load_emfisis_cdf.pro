;Load final RBSP EMFISIS data from https://spdf.gsfc.nasa.gov/pub/data/rbsp/
;Returns tplot variables


;-------------------------
;L2 data "types" to load:
;(NOTE: for L2, type variable needs to combine these (e.g. type='hfr/spectra-burst')
;-------------------------

;hfr
;   spectra-burst
;   spectra-merged
;   spectra
;   waveform
;magnetometer
;   uvw
;wfr
;   spectral-matrix-diagonal-merged
;   spectral-matrix-diagonal
;   spectral-matrix
;   waveform-continuous-burst
;   waveform

;-------------------------
;L3 data "types" to load 
;NOTE: since there's only "magnetometer" data, type='1sec', '4sec', or 'hires'
;with subcategories gei, geo, gse, gsm, sm 
;e.g. type = '4sec-sm'
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
;timespan,'2012-11-01'
;sc = 'rbsp-a'
;lvl = 'l2'
;type = 'hfr/waveform'

;---------
;e.g. load L3 data
;timespan,'2012-11-01'
;sc = 'rbsp-a'
;lvl = 'l3'
;type = '4sec/geo'

;rbsp-a_hfr-waveform_emfisis-l2_20120901_v1.2.9.cdf
;rbsp-a_hfr-spectra-burst_emfisis-l2_20120901_v1.2.6.cdf
;rbsp-a_housekeeping_emfisis-l2_20120104_v1.2.1.cdf  
;rbsp-a_magnetometer_uvw_emfisis-l2_20120830_v1.7.1.cdf
;rbsp-a_magnetometer_hires-gei_emfisis-l3_20120830_v1.7.1.cdf



pro rbsp_load_emfisis_cdf,sc,lvl,type


  scpath = strmid(sc,0,4)+strmid(sc,5,1)

  tr = timerange()


  url = 'https://spdf.gsfc.nasa.gov/pub/data/rbsp/'
  


  dtst = time_string(tr[0],tformat='YYYY-MM-DD')
  yr = strmid(dtst,0,4)
  mm = strmid(dtst,5,2)
  dd = strmid(dtst,8,2)
  
 
  typefn = strsplit(type,'/',/extract)


  
  if lvl eq 'l2' then begin  
    remote_path = url +scpath+'/'+lvl+'/emfisis/'+type+'/'+yr+'/'
    remote_file = sc+'_'+typefn[0]+'-'+typefn[1]+'_emfisis-l2_'+yr+mm+dd+'_v?.?.?.cdf
    local_path = 'data/rbsp/'+scpath+'/'+lvl+'/emfisis/'+type+'/'+yr+'/'
  endif
  if lvl eq 'l3' then begin
    remote_path = url +scpath+'/'+lvl+'/emfisis/magnetometer/'+type+'/'+yr+'/'
    remote_file = sc + '_magnetometer_'+typefn[0]+'-'+typefn[1]+'_emfisis-l3_'+yr+mm+dd+'_v?.?.?.cdf'
    local_path = 'data/rbsp/'+scpath+'/'+lvl+'/emfisis/magnetometer/'+type+'/'+yr+'/'
  endif



  ;Download file and turn into tplot variables
  file = spd_download(remote_file=remote_file,remote_path=remote_path,$
          local_path=local_path,/last_version)
  cdf2tplot,file


stop
end
