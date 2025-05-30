;+
; Scan L1 esvy and vsvy data for timg tag jump of +/-1 sec.
;
; This file should be used by EFW to batch process the data over the whole mission.
;-

pro rbsp_efw_phasef_find_l1_time_tag_irregularity, probe=probe, filename=file, errmsg=errmsg


;---Check inputs.
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif
    
    out_dir = fgetpath(file)
    if file_test(out_dir,/directory) eq 0 then file_mkdir, out_dir
    if file_test(file) eq 1 then file_delete, file  ; overwrite old files.
    ginfo = dictionary($
        'TITLE', 'For RBSP EFW L1 data, times of potential time tag shift', $
        'TEXT', 'Generated by Sheng Tian at the University of Minnesota' )
    cdf_save_setting, ginfo, filename=file
    
    log_file = join_path([out_dir,file_basename(file)+'.log'])
    if file_test(log_file) eq 1 then file_delete, log_file
    ftouch, log_file


;---Constants and settings.
    secofday = 86400d
    prefix = 'rbsp'+probe+'_'
    rbspx = 'rbsp'+probe


;---The L1 data are supposed to be available on the local disk.
    local_root = join_path([default_local_root(),'data','rbsp'])
    version = 'v*'
    foreach l1_type, ['esvy','vsvy'] do begin
        lprmsg, 'Processing L1 data: '+l1_type+' ...', log_file
        lprmsg, '', log_file
        
    ;---Find the mission time range and prepare the files for each day.
        time_range = phasef_get_valid_range(l1_type+'_l1', probe=probe)
        days = make_bins(time_range+[0,-1]*secofday, secofday)
        nday = n_elements(days)
        
        cdf_files = strarr(n_elements(days))
        foreach day, days, day_id do begin
            base_name = rbspx+'_l1_'+l1_type+'_%Y%m%d_'+version+'.cdf'
            local_path = [local_root,rbspx,'l1',l1_type,'%Y']
            local_file = apply_time_to_pattern(join_path([local_path,base_name]), day)
            files = file_search(local_file)
            cdf_files[day_id] = files[-1]   ; Get the highest version, '' if no data for the day.
        endforeach
        
        
    ;---Find time tag gap of 1 sec.
        current_times = list()
        previous_times = list()
        wanted_dtime = 1.
        ; The buffer to save the time tags.
        current_buffer = []
        previous_buffer = []
        foreach cdf_file, cdf_files, day_id do begin
            lprmsg, '', log_file
            lprmsg, 'Processing '+time_string(days[day_id],tformat='YYYY_MMDD')+' ...', log_file
            lprmsg, 'File of day: '+cdf_file+' ...', log_file
            
            ; The time tag of the current day.
            if file_test(cdf_file) eq 0 then begin
                current_buffer = []
            endif else begin
                current_buffer = convert_time(cdf_read_var('epoch', filename=cdf_file), from='epoch16', to='unix')
            endelse
            
            times = [previous_buffer,current_buffer]
            ntime = n_elements(times)
            if ntime gt 2 then begin
                dtime = round(times[1:ntime-1]-times[0:ntime-2])
                index = where(abs(dtime) eq wanted_dtime, count)
                if count ne 0 then begin
                    lprmsg, 'Found '+string(count,format='(I0)')+' times: ', log_file
                    for ii=0,count-1 do begin
                        lprmsg, strjoin(time_string(times[index[ii]:index[ii]+1],tformat='YYYY-MM-DD/hh:mm:ss.ffffff'),'    '), log_file
                    endfor
                    current_times.add, times[index+1], /extract
                    previous_times.add, times[index], /extract
                endif
            endif
            
            ; Prepare for the next day.
            previous_buffer = current_buffer
        endforeach
        
    ;---Save to CDF.
        current_times = current_times.toarray()
        previous_times = previous_times.toarray()
        index = uniq(current_times, sort(current_times))
        current_times = current_times[index]
        previous_times = previous_times[index]
    
        if n_elements(current_times) ne 0 then begin
            time_var = prefix+l1_type+'_current_times'
            data = current_times
            settings = dictionary($
                'FIELDNAM', 'current_unix_time', $
                'UNITS', 'sec', $
                'VAR_TYPE', 'data' )
            cdf_save_var, time_var, value=data, filename=file, cdf_type='CDF_DOUBLE'
            cdf_save_setting, settings, var=time_var, filename=file
        endif
        
        if n_elements(previous_times) ne 0 then begin
            time_var = prefix+l1_type+'_previous_times'
            data = current_times
            settings = dictionary($
                'FIELDNAM', 'previous_unix_time', $
                'UNITS', 'sec', $
                'VAR_TYPE', 'data' )
            cdf_save_var, time_var, value=data, filename=file, cdf_type='CDF_DOUBLE'
            cdf_save_setting, settings, var=time_var, filename=file
        endif        
    endforeach
    
end





;log_file = join_path([homedir(),'rbsp_l1_time_tag_irregularity.txt'])
;if file_test(log_file) eq 1 then file_delete, log_file
;ftouch, log_file
;tab = '    '
;tformat = 'YYYY-MM-DD/hh:mm:ss.ffffff'
;foreach probe, ['a','b'] do begin
;    prefix = 'rbsp'+probe+'_'
;    file = join_path([homedir(),'rbsp'+probe+'_l1_time_tag_irregularity.cdf'])
;    times = convert_time(cdf_read_var('epoch', filename=file),from='epoch',to='unix')
;    pre_times = convert_time(cdf_read_var('epoch_previous', filename=file),from='epoch',to='unix')
;    rbx = 'RB-'+strupcase(probe)
;    dtimes = times-pre_times
;    foreach time, times, ii do begin
;        msg = rbx+tab+time_string(time,tformat=tformat)+tab+time_string(pre_times[ii],tformat=tformat)+tab+string(dtimes[ii])
;        lprmsg, msg, log_file
;    endforeach
;endforeach
;stop


foreach probe, ['a','b'] do begin
    file = join_path([homedir(),'rbsp'+probe+'_l1_time_tag_irregularity.cdf'])
    rbsp_efw_phasef_find_l1_time_tag_irregularity, probe=probe, filename=file, errmsg=errmsg
endforeach
end