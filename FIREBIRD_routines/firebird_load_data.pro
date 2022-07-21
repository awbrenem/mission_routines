;+
; NAME: firebird_load_data
; 
;
; SYNTAX:
;
; PURPOSE: Fetches/loads (multiple days of) FIREBIRD official hires data and stores as tplot variables.
; 
; NOTE: the data are time-corrected AND have the proper flux calibration applied within this routine
;(For testing see firebird_test_time_correction_flux_calibration.pro)
;CAVEAT: I use the integral channel flux values from the FIREBIRD file, b/c I don't know how to properly calibrate these.
;
;
; Usage: timespan,'2015-06-11',2,/days
;        firebird_load_data,'4'
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
  ndays_load = floor((tr[1]-tr[0])/86400)




  for i=0,ndays_load -1 do begin

  
    date = time_string(tr[0]+i*86400,/date_only,tformat='YYYYMMDD')
    yyyy = strmid(date,0,4)
    mm = strmid(date,4,2)
    dd = strmid(date,6,2)
  
    type = 'hires'
    type2 = 'Hires'
  
    url = 'http://solar.physics.montana.edu/FIREBIRD_II/Data/' + cubesat + '/' + type + '/'
  
    fn = cubesat2 + '_' + type2 + '_' + yyyy+'-'+mm+'-'+dd+'_L2.txt'
  
    dprint,dlevel=3,verbose=verbose,relpathnames,/phelp
  
  
    local_path = '/Users/abrenema/data/firebird/' + cubesat2 + '/' + yyyy + '/'
  
    files = spd_download(remote_path=url,remote_file=fn,$
    local_path=local_path,$
    /last_version)
  
  
  
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
  
  
    ;-------------------------------------------------------------  
    ;Add the time correction to the hires data 
    ;-------------------------------------------------------------
  
    time = time_double(data.time) + double(data.count_time_correction)
    
    ;-------------------------------------------------------------  
    ;Properly calibrate counts to flux. 
    ;NOTE: we won't return the flux values in this file b/c they are determined with the incorrect
    ;geometric factor.
    ;-------------------------------------------------------------
  
  
  
    x = firebird_get_calibration_counts2flux(time_string(tr[0]+i*86400,tformat='YYYY-MM-DD'),strmid(cubesat2,2,1))
  
  
    dE = x.ENERGY_RANGE_COLLIMATED[*,1] - x.ENERGY_RANGE_COLLIMATED[*,0]
    col_flux1 = data.col_counts1/((x.cadence/1000.)*dE[0]*x.G_FACTOR_COLLIMATED[0])
    col_flux2 = data.col_counts2/((x.cadence/1000.)*dE[1]*x.G_FACTOR_COLLIMATED[1])
    col_flux3 = data.col_counts3/((x.cadence/1000.)*dE[2]*x.G_FACTOR_COLLIMATED[2])
    col_flux4 = data.col_counts4/((x.cadence/1000.)*dE[3]*x.G_FACTOR_COLLIMATED[3])
    col_flux5 = data.col_counts5/((x.cadence/1000.)*dE[4]*x.G_FACTOR_COLLIMATED[4])
  
    dE = x.ENERGY_RANGE_SURFACE[*,1] - x.ENERGY_RANGE_SURFACE[*,0]
    sur_flux1 = data.sur_counts1/((x.cadence/1000.)*dE[0]*x.G_FACTOR_SURFACE[0])
    sur_flux2 = data.sur_counts2/((x.cadence/1000.)*dE[1]*x.G_FACTOR_SURFACE[1])
    sur_flux3 = data.sur_counts3/((x.cadence/1000.)*dE[2]*x.G_FACTOR_SURFACE[2])
    sur_flux4 = data.sur_counts4/((x.cadence/1000.)*dE[3]*x.G_FACTOR_SURFACE[3])
    sur_flux5 = data.sur_counts5/((x.cadence/1000.)*dE[4]*x.G_FACTOR_SURFACE[4])
  
  
  
    csstr = strlowcase(cubesat2)
  
    store_data,csstr+'_fb_col_hires_flux',time,double([[col_flux1],[col_flux2],[col_flux3],[col_flux4],[col_flux5],[data.col_flux6]])
    store_data,csstr+'_fb_col_hires_counts',time,double([[data.col_counts1],[data.col_counts2],[data.col_counts3],[data.col_counts4],[data.col_counts5],[data.col_counts6]])
    store_data,csstr+'_fb_sur_hires_flux',time,double([[sur_flux1],[sur_flux2],[sur_flux3],[sur_flux4],[sur_flux5],[data.sur_flux6]])
    store_data,csstr+'_fb_sur_hires_counts',time,double([[data.sur_counts1],[data.sur_counts2],[data.sur_counts3],[data.sur_counts4],[data.sur_counts5],[data.sur_counts6]])
  
    store_data,csstr+'_fb_geolat_from_hiresfile',time,double(data.lat) ;Geographic lat
    store_data,csstr+'_fb_geolon_from_hiresfile',time,double(data.lon) ;Geographic long
    store_data,csstr+'_fb_alt_from_hiresfile',time,double(data.alt)
    store_data,csstr+'_fb_mlt_from_hiresfile',time,double(data.mlt)
    store_data,csstr+'_fb_mcilwainL_from_hiresfile',time,double(data.mcilwainl)
  
  
  
    ylim,csstr+'_fb_col_hires_flux',0.1,1000,1
    options,csstr+'_fb_col_hires_flux_?','psym',-5
    options,csstr+'_fb_sur_hires_flux_?','psym',-5
  


    ;final names of the tplot variables (don't delete)
    tnames = csstr + '_fb_'+ ['col_hires_flux','col_hires_counts','sur_hires_flux','sur_hires_counts','geolat_from_hiresfile','geolon_from_hiresfile','alt_from_hiresfile','mlt_from_hiresfile','mcilwainL_from_hiresfile']


    if ndays_load gt 1 then begin


      ;Rename to temporary variable
      if i eq 0 then for j=0,n_elements(tnames)-1 do copy_data,tnames[j],tnames[j]+'_fin'


      ;If both renamed and original, then combine
      if i gt 0 then begin
        for j=0, n_elements(tnames)-1 do begin


          get_data,tnames[j],data=tntmp
          get_data,tnames[j]+'_fin',data=tnfin
          store_data,tnames[j]+'_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]}



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




end
