;+
; Read density.
;-

pro rbsp_efw_phasef_read_density, time, probe=probe, boom_pairs=boom_pairs, dmin=dmin, dmax=dmax


    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'data','rbsp'])
    if n_elements(remote_root) eq 0 then remote_root = join_path([rbsp_efw_phasef_get_server()])
    if n_elements(version) eq 0 then version = 'v01'
    if n_elements(datatype) eq 0 then datatype = 'density'
    if n_elements(boom_pairs) eq 0 then boom_pairs = ['12','34','13','14','23','24']
    if n_elements(dmin) eq 0 then dmin = 10d
    if n_elements(dmax) eq 0 then dmax = 3000d


;---Init settings.
    type_dispatch = hash()
    valid_range = (probe eq 'a')? time_double(['2012-09-08','2019-10-15']): time_double(['2012-09-08','2019-07-17'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_efw_density_uh_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'density_uh','%Y']
    remote_path = [remote_root,'density_uh',rbspx,'%Y']

    type_dispatch['density'] = dictionary($
        'pattern', dictionary($
            'remote_file', join_path([remote_path,base_name]), $
            'remote_index_file', join_path([remote_path,'']), $
            'local_file', join_path([local_path,base_name]), $
            'local_index_file', join_path([local_path,default_index_file()])), $
        'valid_range', time_double(valid_range), $
        'cadence', 'day', $
        'extension', fgetext(base_name), $
        'var_list', list($
            dictionary($
                'in_vars', rbspx+['_density_'+boom_pairs], $
                'time_var_name', 'epoch', $
                'time_var_type', 'epoch')))

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
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            rbsp_efw_phasef_read_density_gen_file, file_time, probe=probe, filename=local_file
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif


;---Read data from files and save to memory.
    read_files, time, files=files, request=request


;---Apply flags.
    rbsp_efw_read_flags, time, probe=probe
    prefix = 'rbsp'+probe+'_'
    get_data, prefix+'density_'+boom_pairs[0], common_times
    flag_var = prefix+'efw_flags'
    interp_time, flag_var, common_times
    flags = get_var_data(flag_var)
    flag_names = get_setting(flag_var, 'labels')
    fillval = !values.f_nan
    foreach boom_pair, boom_pairs do begin
        var = prefix+'density_'+boom_pair
        density = get_var_data(var)
        foreach flag_type, ['charging','charging_extreme'] do begin
            index = where(flag_names eq flag_type+'_'+boom_pair, count)
            if count ne 0 then begin
                index = where(flags[*,index] eq 1, count)
                if count ne 0 then density[index] = fillval
            endif
        endforeach

        index = lazy_where(density, '][', [dmin,dmax], count=count)
        if count ne 0 then density[index] = fillval

        store_data, var, common_times, density, limits={ylog:1, ytitle:'(cm!U-3)', labels:'Density '+boom_pair}
    endforeach

end


time_range = time_double(['2014-08-28','2014-08-29'])
probe = 'a'
rbsp_efw_phasef_read_density, time_range, probe=probe
end
