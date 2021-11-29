;+
; Read flags.
;-

pro rbsp_efw_read_flags, time, probe=probe, errmsg=errmsg

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'data','rbsp'])
    if n_elements(remote_root) eq 0 then remote_root = join_path([rbsp_efw_phasef_get_server()])
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(datatype) eq 0 then datatype = 'flags'

;---Init settings.
    type_dispatch = hash()
    valid_range = (probe eq 'a')? time_double(['2012-09-08','2019-02-23']): time_double(['2012-09-08','2019-07-17'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_efw_flags_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'flags','%Y']
    remote_path = [remote_root,rbspx,'flags',rbspx,'%Y']

    request = dictionary($
        'pattern', dictionary($
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,'']), $
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list())

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            rbsp_efw_read_flags_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

    if n_elements(files) eq 0 then begin
        errmsg = 'Failed to find file ...'
        return
    endif
    cdf2tplot, files
    flag_names = cdf_read_var('flag_labels', filename=files[0])
    prefix = 'rbsp'+probe+'_'
    var = prefix+'efw_flags'
    add_setting, var, /smart, dictionary($
        'display_type', 'stack', $
        'yrange', [-0.2,1.2], $
        'labels', flag_names )

end


time_range = time_double(['2013-01-01','2013-01-02'])
probes = ['a']
foreach probe, probes do $
    rbsp_efw_read_flags, time_range, probe=probe

end
