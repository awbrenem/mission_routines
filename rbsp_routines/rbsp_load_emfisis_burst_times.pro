;+
; NAME: RBSP_LOAD_EMFISIS_BURST_TIMES
;
; SYNTAX:
;	rbsp_load_emfisis_burst_times,probe='a b'
;
; PURPOSE: Loads EMFISIS Burst data availability. Creates a tplot
;	variable indicating available times.
;
;	date = 'yyyy-mm-dd'
;	probe = 'a' or 'b'
;
; VERSION:
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2018-12-06 09:25:24 -0800 (Thu, 06 Dec 2018) $
;   $LastChangedRevision: 26261 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/rbsp_load_emfisis_burst_times.pro $
;
;-

pro rbsp_load_emfisis_burst_times,probe




	remote_data_dir = 'https://emfisis.physics.uiowa.edu/events/'
	subdir = 'rbsp-'+probe+'/burst/'
	local_path = '/Users/abrenema/data/rbsp/emfisis/burst_availability/rbsp-'+probe+'/burst/'


	tr = timerange()
	ndays_load = floor((tr[1]-tr[0])/86400)


	for i=0,ndays_load -1 do begin

	  dtst = time_string(tr[0]+i*86400,tformat='YYYY-MM-DD')
  	day = strmid(dtst,0,4)+strmid(dtst,5,2)+strmid(dtst,8,2)
  
  
 
  
  	;grab the online EMFISIS list of burst times and read in
  	fn = 'rbsp-'+probe+'_burst_times_'+day+'.txt'
  	files = spd_download(remote_path=remote_data_dir+subdir,remote_file=fn,$
    local_path=local_path,$
    /last_version)
  
  	resulttst = FILE_TEST(files)

    if resulttst then begin
    
    	openr,lun,local_path+fn,/get_lun
    	jnk = ''
    	readf,lun,jnk
    
    	lines = strarr(50000.)
    	q=0.
    	while not eof(lun) do begin
    		readf,lun,jnk
    		lines[q] = jnk
    		q++
    	endwhile
    	close,lun
    	free_lun,lun
    
    
    
    	;remove blank strings
    	goo = where(lines eq '')
    	lines = lines[0:goo[0]-1]
    
    	t0 = dblarr(n_elements(lines))
    	t1 = dblarr(n_elements(lines))
    	type = strarr(n_elements(lines))
    
    	;Extract the times
    	for q=0,n_elements(lines)-1 do begin
    		tmp = strsplit(lines[q],',',/extract)
    		t0[q] = time_double(tmp[0])
    		t1[q] = time_double(tmp[1])
    		type[q] = tmp[3]
    	endfor
    
    
    	;create artificial array of times.
;      timestmp = time_double(dtst) + dindgen(100.*86400.)/99.
      timestmp = time_double(dtst) + dindgen(86400.)
    	valstmp = fltarr(n_elements(timestmp))
    
    
    	for b=0,n_elements(t0)-1 do begin $
    			goo = where((timestmp ge t0[b]) and (timestmp le t1[b])) & $
    			if goo[0] ne -1 then valstmp[goo] = 1.
    	endfor
    
    	store_data,'rbsp'+probe+'_emfisis_burst',timestmp,valstmp





    	if ndays_load gt 1 then begin

    	  ;Rename
        if i eq 0 then copy_data,'rbsp'+probe+'_emfisis_burst','rbsp'+probe+'_emfisis_burst_fin'


    	  ;If both renamed and original, then combine
    	  if i gt 0 then begin

  	      get_data,'rbsp'+probe+'_emfisis_burst',data=tntmp
  	      get_data,'rbsp'+probe+'_emfisis_burst_fin',data=tnfin
  	      store_data,'rbsp'+probe+'_emfisis_burst_fin',data={x:[tnfin.x,tntmp.x],y:[tnfin.y,tntmp.y]}

    	  endif  ;i>0
    	  store_data,'rbsp'+probe+'_emfisis_burst',/del
    	endif  ;if more than 1 day to load



    endif ;for a daily file to load
   endfor  ;for each day to load


   ;Final rename of variables
   if ndays_load gt 1 then begin
     ;for j=0, n_elements(tnames)-1 do copy_data,tnames[j]+'_fin',tnames[j]
     copy_data,'rbsp'+probe+'_emfisis_burst_fin','rbsp'+probe+'_emfisis_burst'
     store_data,'rbsp'+probe+'_emfisis_burst_fin',/del
   endif



   ylim,'rbsp'+probe+'_emfisis_burst',0,2




  
end
