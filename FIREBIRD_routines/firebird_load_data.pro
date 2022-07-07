;+
; NAME: firebird_load_data
;
; SYNTAX:
;
; PURPOSE: Fetches/loads FIREBIRD official hires data and stores as tplot variables
;
; Usage: timespan,'2019-01-24'
;        firebird_load_data,'3'
;
; INPUT: N/A
;
; OUTPUT: N/A
;
; KEYWORDS:
;   cubesat = '3' or '4'
;   fileexists --> returns a 1 if the file has been found
; HISTORY:
;	Created Dec 2017, Aaron Breneman
; NOTES:
;
; VERSION:
;
;-

;Test
;2016-01-20T19:43:35.301500
;mlt_fb = 10.413636009088213
;lon_fb = -121.6999989029282
;lat_fb = 56.788301648640505



pro firebird_load_data,cubesat,plot=plot,fileexists=fileexists

  ;default. Can change later if file not found
  fileexists = 1

  if cubesat eq '3' then cubesat = 'FU_3' else cubesat = 'FU_4'
  cubesat2 = strmid(cubesat,0,2) + strmid(cubesat,3,1)

  tr = timerange()
  date = time_string(tr[0],/date_only,tformat='YYYYMMDD')
  yyyy = strmid(date,0,4)
  mm = strmid(date,4,2)
  dd = strmid(date,6,2)

  type = 'hires'
  type2 = 'Hires'

  url = 'http://solar.physics.montana.edu/FIREBIRD_II/Data/' + cubesat + '/' + type + '/'

  fn = cubesat2 + '_' + type2 + '_' + yyyy+'-'+mm+'-'+dd+'_L2.txt'

  dprint,dlevel=3,verbose=verbose,relpathnames,/phelp

  ;local_path = '/Users/aaronbreneman/Desktop/Research/RBSP_Firebird_microburst_conjunctions_all/firebird/'
  local_path = '/Users/abrenema/data/firebird/' + cubesat2 + '/' + yyyy + '/'

  files = spd_download(remote_path=url,remote_file=fn,$
  local_path=local_path,$
  /last_version)

stop

;---------------------------------------------------------------
;Read file
  ft = [7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3]
  fl = [0,27,45,64,84,105,109,113,117,121,125,129,133,137,142,147,151,155,159,163,167,171,175,179,183,187,205,223,242,261,280,299,304,308]
  fg = indgen(n_elements(fl))
  fns = ['time','col_flux1','col_flux2','col_flux3','col_flux4','col_flux5',$
  'col_flux6','sur_flux1','sur_flux2','sur_flux3','sur_flux4','sur_flux5','sur_flux6','col_counts1',$
  'col_counts2','col_counts3','col_counts4','col_counts5','col_counts6','sur_counts1',$
  'sur_counts2','sur_counts3','sur_counts4','sur_counts5','sur_counts6','count_time_correction',$
  'mcilwainl','lat','lon','alt','mlt','kp','flag','loss_cone_type']


  template = {VERSION: 1.0,$
              DATASTART: 125L,$
              DELIMITER: 32b,$
              MISSINGVALUE: !values.f_nan,$
              COMMENTSYMBOL: '',$
              FIELDCOUNT:34L,$
              FIELDTYPES:ft,$
              FIELDNAMES:fns,$
              FIELDLOCATIONS:fl,$
              FIELDGROUPS:fg}



  ;Check to see if file has been downloaded.
  file_exists = FILE_TEST(local_path+fn)
  if ~file_exists then begin
    fileexists = 0
    print,'NO FIREBIRD DATA FOR THIS DATE....RETURNING'
    return
  endif





  data = read_ascii(local_path+fn,template=template)

  time = time_double(data.time)

  csstr = strlowcase(cubesat2)

  store_data,csstr+'_fb_col_hires_flux',time,double([[data.col_flux1],[data.col_flux2],[data.col_flux3],[data.col_flux4],[data.col_flux5],[data.col_flux6]])
  store_data,csstr+'_fb_col_hires_counts',time,double([[data.col_counts1],[data.col_counts2],[data.col_counts3],[data.col_counts4],[data.col_counts5],[data.col_counts6]])
  store_data,csstr+'_fb_sur_hires_flux',time,double([[data.sur_flux1],[data.sur_flux2],[data.sur_flux3],[data.sur_flux4],[data.sur_flux5],[data.sur_flux6]])
  store_data,csstr+'_fb_sur_hires_counts',time,double([[data.sur_counts1],[data.sur_counts2],[data.sur_counts3],[data.sur_counts4],[data.sur_counts5],[data.sur_counts6]])

  store_data,csstr+'_fb_geolat_from_hiresfile',time,double(data.lat) ;Geographic lat
  store_data,csstr+'_fb_geolon_from_hiresfile',time,double(data.lon) ;Geographic long
  store_data,csstr+'_fb_alt_from_hiresfile',time,double(data.alt)
  store_data,csstr+'_fb_mlt_from_hiresfile',time,double(data.mlt)
  store_data,csstr+'_fb_mcilwainL_from_hiresfile',time,double(data.mcilwainl)



  ylim,csstr+'_fb_col_hires_flux',0.1,1000,1

;Campaign 21 energy bins for FU_3
;201.2 - 250.2 keV","250.2 - 299.2 keV","299.2 - 348.2 keV","348.2 - 446.2keV","446.2 - 1055.2 keV"
;"ENERGY_WIDTHS": ["49.0 keV","49.0 keV","49.0 keV","98.0 keV","609.0 keV"],
;Campaign 21 energy bins for FU_4
;201.0 - 246.5 keV","246.5 - 301.1 keV","301.1 - 346.6 keV","346.6 - 446.7 keV","446.7 - 983.6 keV
;"ENERGY_WIDTHS": ["45.5 keV","54.6 keV","45.5 keV","100.1 keV","536.9 keV"]


;  probe = 'a'
;  rbsp_load_efw_spec,probe=probe,type='calibrated',/pt
;  rbsp_load_efw_spec,probe=probe,type='calibrated',/pt

;  get_data,'fu3_fb_col_hires_flux_0',data=d
;  period = d.x[1]-d.x[0]
;  rbsp_detrend,csstr+'_fb_col_hires_flux_?',2*period
;  rbsp_detrend,csstr+'_fb_sur_hires_flux_?',2*period

  options,csstr+'_fb_col_hires_flux_?','psym',-5
  options,csstr+'_fb_sur_hires_flux_?','psym',-5

  if KEYWORD_SET(plot) then begin
    rbsp_efw_init
    split_vec,csstr+'_fb_col_hires_flux'
    ylim,csstr+'_fb_col_hires_flux_?',0,0,0
    tplot,['rbsp'+probe+'_efw_64_spec[0:4]',csstr+'_fb_col_hires_flux_[0-5]']
    stop
    split_vec,csstr+'_fb_sur_hires_flux'
    ylim,csstr+'_fb_sur_hires_flux_?',0,0,0
;    tplot,csstr+'_fb_sur_hires_flux_[0-5]'
    tplot,['rbsp'+probe+'_efw_64_spec[0:4]',csstr+'_fb_sur_hires_flux_[0-5]']
    stop
  endif

;upper = srt/2.
;rbsp_decimate,'fu3_fb_col_hires_flux_0',upper=upper,level=lvl,newname='tst'
;rbsp_decimate,'fu3_fb_col_hires_flux_0',level=1,newname='tst'




end
