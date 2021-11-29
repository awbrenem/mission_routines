;+
; Read EFW L2 data.
;
; datatype can be 'e-hires-uvw','e-spinfit-mgse','esvy_despun','fbk','spec','vsvy-hires'
;-

pro rbsp_efw_read_l2, tr, probe=probe, datatype=datatype, trange=trange, $
    level=level, verbose=verbose, downloadonly=downloadonly, $
     cdf_data=cdf_data,get_support_data=get_support_data, $
     tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
     varformat=varformat, valid_names = valid_names, files=files, $
     type=type, _extra = _extra

    rbsp_efw_init
    vb = keyword_set(verbose) ? verbose : 0
    vb = vb > !rbsp_efw.verbose

    if n_elements(probe) eq 0 then probe = 'a'
    if n_elements(version) eq 0 then version = 'v*'
    if n_elements(trange) ne 0 then time_range = trange
    if n_elements(tr) ne 0 then time_range = tr
    if n_elements(time_range) eq 0 then time_range = timerange()
    if size(time_range[0],/type) eq 7 then time_range = time_double(time_range)
    if n_elements(datatype) eq 0 then begin
        dprint, 'No datatype ...', verbose=vb
        return
    endif
    data_types = ['e-hires-uvw','e-spinfit-mgse','esvy_despun','fbk','spec','vsvy-hires']
    index = where(data_types eq datatype[0], count)
    if count eq 0 then begin
        dprint, 'Invalid datatype: '+datatype[0]+' ...', verbose=vb
        return
    endif


    rbspx = 'rbsp'+probe
    base = rbspx+'_efw-l2_'+datatype+'_YYYYMMDD_'+version+'.cdf'
    local_root = !rbsp_efw.local_data_dir
    local_path = [local_root,rbspx,'efw','l2',datatype,'YYYY',base]
    remote_root = rbsp_efw_remote_root()
    remote_path = [remote_root,rbspx,rbsp_efw_remote_sub_dirs(level='l2',datatype=datatype),base]
    local_files = file_dailynames(file_format=join_path(local_path), trange=time_range)
    remote_files = file_dailynames(file_format=join_path(remote_path), trange=time_range)


    nfile = n_elements(local_files)
    for file_id=0,nfile-1 do begin
        url = remote_files[file_id]
        spd_download_expand, url, last_version=1, $
            ssl_verify_peer=ssl_verify_peer, ssl_verify_host=ssl_verify_host, _extra=_extra
        base = file_basename(url)
        local_file = join_path([file_dirname(local_files[file_id]),base])
        tmp = spd_download_file(url=url, filename=local_file)
        local_files[file_id] = local_file
    endfor


    suffix = ''
    prefix = rbspx+'_efw_'
    cdf2tplot, file=local_files, all=0, prefix=prefix, suffix=suffix, verbose=vb, $
        tplotnames=tns, /convert_int1_to_int2, get_support_data=0
        
        
    if datatype eq 'spec' then begin
        vars = tnames(prefix+'spec*')
        options, vars, 'ylog', 1
        options, vars, 'ytitle', 'Freq [Hz]'
        foreach var, vars do begin
            get_data, var, dlimit=dlim
            options, var, 'ztitle', '['+dlim.ysubtitle+']'
            options, var, 'ysubtitle', ''
        endforeach
    endif else if datatype eq 'fbk' then begin
        vars = tnames(prefix+'fbk*')
        options, vars, 'ylog', 1
        options, vars, 'ytitle', 'Freq!Ccenter!C[Hz]'
        foreach var, vars do begin
            get_data, var, dlimit=dlim
            options, var, 'ztitle', dlim.ysubtitle
            options, var, 'ysubtitle', ''
        endforeach

    endif

end


; Set the time and probe for loading data.
time_range = ['2013-01-01','2013-01-03']
probe = 'b'

; Load the L2 data.
rbsp_efw_read_l2, time_range, probe=probe, datatype='esvy_despun'
rbsp_efw_read_l2, time_range, probe=probe, datatype='vsvy-hires'

prefix = 'rbsp'+probe+'_efw_'
vars = prefix+[$

    ; The E field in mGSE.
    'efield_mgse', $

    ; The single-ended boom potential.
    'vsvy', $

    ; The flags.
    'efw_qual' ]

; Plot the variables.
tplot, vars, trange=time_range
end
