;+
; Read Themis EFI data.
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

pro themis_read_efi, time, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root, $
    coordinate=coord

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 86400d*120
    if n_elements(probe) eq 0 then probe = 'x'
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'data','themis'])
    if n_elements(remote_root) eq 0 then remote_root = 'https://cdaweb.sci.gsfc.nasa.gov/pub/data/themis'
    if n_elements(version) eq 0 then version = 'v[0-9]{2}'
    if n_elements(coord) eq 0 then coord = 'gsm'

;---Init settings.
    type_dispatch = hash()
    thx = 'th'+probe
    ; Level 2.
    base_name = thx+'_l2_efi_%Y%m%d_'+version+'.cdf'
    local_path = [local_root,thx,'l2','efi','%Y']
    remote_path = [remote_root,thx,'l2','efi','%Y']
    ; efs, 3sec resolution.
    ; eff, 1/8 sec resolution.
    foreach key, ['efs','eff'] do begin
        type_dispatch['l2%'+key] = dictionary($
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
                    'in_vars', thx+'_'+key+'_dot0_'+coord, $
                    'time_var_name', thx+'_'+key+'_dot0_time', $
                    'time_var_type', 'unix')))
    endforeach

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

themis_read_efi, /print_datatype
time = time_double(['2013-10-30/23:00','2013-10-31/06:00'])
themis_read_efi, time, id='l2%efs', probe='d'
end
