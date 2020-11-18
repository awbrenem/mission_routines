;Load FIREBIRD context (survey) CDF files that Aaron created at
;http://rbsp.space.umn.edu/firebird/FU?/YYYY/


;e.g. FU3_context_20150204_v01.cdf

;These have the tplot variables:
;1 D0  (counts, not flux)
;2 D1  (counts, not flux)
;3 Alt
;4 Count_Time_Correction
;5 Flag
;6 Lat
;7 Lon
;8 Loss_cone_type
;9 MLT
;10 McIlwainL
;11 kp

;These files are created with
;firebird_create_cdf_files


;Usage
;timespan,'2018-08-15'
;firebird_load_context_data_cdf_file,'FU3'

;options are FU3 or FU4



pro firebird_load_context_data_cdf_file,sc

  rbsp_efw_init

;  sc = 'FU3'
;  datetime = '2018-08-15'
;  timespan,datetime

  tr = timerange()
  datetime = strmid(time_string(tr[0]),0,10)


  year = strmid(datetime,0,4)
  mn = strmid(datetime,5,2)
  dy = strmid(datetime,8,2)

  fn = sc + '_context_'+year+mn+dy+'_v01.cdf'


  ;Grab path to save data 
  spawn,'pwd',tmp
  strvals = strsplit(tmp,'/',/extract)

;local folder
;  folder = '/Users/aaronbreneman/data/firebird/'+sc+'/'+year+'/'
  folder = '/'+strvals[0]+'/'+strvals[1]+'/data/firebird/'+sc+'/' + year + '/'



;  url = 'http://rbsp.space.umn.edu/kersten/data/firebird/'     ;FU?/YYYY/
  url = 'http://rbsp.space.umn.edu/firebird/'     ;FU?/YYYY/


  path = sc+'/'+year+'/'+fn


  file_loaded = spd_download(remote_file=url+path,$
                local_path=folder,$
                /last_version)



  cdf2tplot,file_loaded


end
