;+
; Read RBSP EFW data.
;
; time. A time or a time range in ut time. Set time to find files
;   automatically, or set files to read data in them directly.
; id=. A string sets the data type to read. Check supported ids by setting
;   print_datatype.
; print_datatype=. A boolean. Set to print all supported ids.
; probe=. A string set the probe to read data for.
; local_root=. A string to set the local root directory.
; remote_root=. A string to set the remote root directory.
; local_files=. A string or an array of N full file names. Set to fine
;   tuning the files to read data from.
; file_times=. An array of N times. Set to fine tuning the times of the files.
;-
pro rbsp_read_efw, time, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root

    compile_opt idl2
    on_error, 0
    errmsg = ''


;---Check inputs.
    sync_threshold = 86400d*120
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'data','rbsp'])
    if n_elements(remote_root) eq 0 then remote_root = 'https://cdaweb.sci.gsfc.nasa.gov/pub/data/rbsp'
    if n_elements(version) eq 0 then version = 'v[0-9]{2}'


;---Init settings.
    type_dispatch = hash()
    rbspx = 'rbsp'+probe
    ; Level 1 burst 1.
    foreach key, ['vb1','mscb1'] do begin
        base_name = rbspx+'_l1_'+key+'_%Y%m%d_'+version+'.cdf'
        local_path = [local_root,rbspx,'efw','l1',key,'%Y']
        remote_path = [remote_root,rbspx,'l1','efw',key,'%Y']
        type_dispatch['l1%'+key] = dictionary($
            'pattern', dictionary($
                'local_file', join_path([local_path,base_name]), $
                'local_index_file', join_path([local_path,default_index_file(/sync)]), $
                'remote_file', join_path([remote_path,base_name]), $
                'remote_index_file', join_path([remote_path,''])), $
            'sync_threshold', sync_threshold, $
            'cadence', 'day', $
            'extension', fgetext(base_name), $
            'var_list', list($
                dictionary($
                    'in_vars', [key], $
                    'out_vars', [rbspx+'_efw_'+key], $
                    'time_var_name', 'epoch', $
                    'time_var_type', 'epoch16')))
    endforeach

    base_name = rbspx+'_l1_esvy_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efw','l1','esvy','%Y']
    remote_path = ['http://themis.ssl.berkeley.edu/data/rbsp', $
        rbspx,'l1','esvy','%Y']
    type_dispatch['l1%esvy'] = dictionary($
        'pattern', dictionary($
        'local_file', join_path([local_path,base_name]), $
        'local_index_file', join_path([local_path,default_index_file(/sync)]), $
        'remote_file', join_path([remote_path,base_name]), $
        'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['esvy'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))

    base_name = rbspx+'_l1_vsvy_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'l1','vsvy','%Y']
    remote_path = ['http://themis.ssl.berkeley.edu/data/rbsp', $
        rbspx,'l1','vsvy','%Y']
    type_dispatch['l1%vsvy'] = dictionary($
        'pattern', dictionary($
        'local_file', join_path([local_path,base_name]), $
        'local_index_file', join_path([local_path,default_index_file(/sync)]), $
        'remote_file', join_path([remote_path,base_name]), $
        'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['vsvy'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))

    ; Level 2.
    base_name = rbspx+'_efw-l2_e-hires-uvw_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efw','l2','e-highres-uvw','%Y']
    remote_path = [remote_root,rbspx,'l2','efw','e-highres-uvw','%Y']
    type_dispatch['l2%uvw'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)]), $
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['e_hires_uvw'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))
    base_name = rbspx+'_efw-l2_vsvy-hires_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efw','l2','vsvy-highres','%Y']
    remote_path = [remote_root,rbspx,'l2','efw','vsvy-highres','%Y']
    type_dispatch['l2%vsvy-highres'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)]), $
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['vsvy'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))
    type_dispatch['l2%vsvy-highres2'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)]), $
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['vsvy'], $
                'time_var_name', 'epoch_v', $
                'time_var_type', 'epoch16')))
    base_name = rbspx+'_efw-l2_e-spinfit-mgse_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efw','l2','e-spinfit-mgse','%Y']
    remote_path = [remote_root,rbspx,'l2','efw','e-spinfit-mgse','%Y']
    type_dispatch['l2%spinfit'] = dictionary($
        'pattern', dictionary($
        'local_file', join_path([local_path,base_name]), $
        'local_index_file', join_path([local_path,default_index_file(/sync)]), $
        'remote_file', join_path([remote_path,base_name]), $
        'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['efield_spinfit_mgse'], $
                'out_vars', [rbspx+'_e_mgse'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))

    ; Level 3.
    base_name = rbspx+'_efw-l3_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'efw','l3','%Y']
    remote_path = [remote_root,rbspx,'l3','efw','%Y']
    type_dispatch['l3%efw'] = dictionary($
        'pattern', dictionary($
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file(/sync)]), $
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,''])), $
        'sync_threshold', sync_threshold, $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', ['efield_inertial_frame_mgse'], $
                'out_vars', [rbspx+'_e_mgse'], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch16')))

    if keyword_set(print_datatype) then begin
        print, 'Suported data type: '
        ids = type_dispatch.keys()
        foreach id, ids do print, '  * '+id
        return
    endif


;---Dispatch patterns.
    if n_elements(datatype) eq 0 then begin
        errmsg = handle_error('No input datatype ...')
        return
    endif
    if not type_dispatch.haskey(datatype) then begin
        errmsg = handle_error('Do not support type '+datatype+' yet ...')
        return
    endif
    request = type_dispatch[datatype]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)

;---Read data from files and save to memory.
    read_files, time, files=files, request=request

end


rbsp_read_efw, /print_datatype
time = time_double(['2013-06-07/04:52','2013-06-07/05:02'])
rbsp_read_efw, time, probe='b', id='l3%efw'
end
