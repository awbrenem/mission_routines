
;Load Efield fit files at

;http://rbsp.space.umn.edu/rbsp_efw/test_perigee_correction/maneuver_fit/


;e.g. 	rbspb_fit_e_mgse_2018_v01.cdf


;These have the tplot variables:
; rbspb_dey 
; rbspb_dez

;Each of these variables has 3 components:
;[n,0]: old (data after running spin fit)
;[n,1]: fit  
;[n,2]: new (=old-fit). 


;These CDF files were created by Sheng Tian, 2020.


pro rbsp_load_efit_cdf_file,sc;,testing=testing


;    sc = 'a'
 ;   testing = 1
  ;  date = '2013-07-15'
   ; timespan,date
  rbsp_efw_init

  tr = timerange()

    ;Extract year 
 
  datetime = strmid(time_string(tr[0]),0,10)


  year = strmid(datetime,0,4)
;  mn = strmid(datetime,5,2)
;  dy = strmid(datetime,8,2)


  fn = 'rbsp'+sc+'_fit_e_mgse_'+year+'_v01.cdf'

  if ~keyword_set(folder) then folder = !rbsp_efw.local_data_dir + $
                                       'rbsp' + strlowcase(sc[0]) + path_sep() + $
                                       'e_fit' + path_sep() + $
                                       year + path_sep()


    
    url = 'http://rbsp.space.umn.edu/rbsp_efw/test_perigee_correction/maneuver_fit/'



  file_loaded = spd_download(remote_file=url+fn,$
                local_path=folder,$
                /last_version)

  cdf2tplot,file_loaded

;  if ~KEYWORD_SET(testing) then $
 ;   cdf2tplot,file_loaded else $
  ;  tplot_restore,filename='~/Desktop/rbspa_test_perigee_correction_2013_0715_v01.tplot'



end
