;+
; Read EFW burst data. Save data in rbspx_efw_mscb1_<coord>.
;
; datatype can be 'mscb1','mscb2'.
;
; coord=. A lower case string sets the coordinate of the output data. 'mgse' by default.
;-


pro rbsp_efw_read_burst_bfield, tr, probe=probe, datatype=datatype, trange=trange, $
    level=level, verbose=verbose, downloadonly=downloadonly, $
     cdf_data=cdf_data,get_support_data=get_support_data, $
     tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
     varformat=varformat, valid_names = valid_names, files=files, $
     type=type, _extra = _extra

    rbsp_efw_init
    vb = keyword_set(verbose) ? verbose : 0
    vb = vb > !rbsp_efw.verbose

    if n_elements(coord) eq 0 then coord = 'mgse'
    if n_elements(probe) eq 0 then probe = 'a'
    if n_elements(version) eq 0 then version = 'v*'
    if n_elements(trange) ne 0 then time_range = trange
    if n_elements(tr) ne 0 then time_range = tr
    if n_elements(time_range) eq 0 then time_range = timerange()
    if size(time_range[0],/type) eq 7 then time_range = time_double(time_range)
    timespan, time_range[0], total(time_range*[-1,1]), /seconds
    if n_elements(datatype) eq 0 then datatype = 'mscb1'
    data_types = ['mscb1','mscb2']
    index = where(data_types eq datatype[0], count)
    if count eq 0 then begin
        dprint, 'Invalid datatype: '+datatype[0]+' ...', verbose=vb
        return
    endif


    rbspx = 'rbsp'+probe
    base = rbspx+'_l1_'+datatype+'_YYYYMMDD_hhmm_'+version+'.cdf'
    local_root = !rbsp_efw.local_data_dir
    local_path = [local_root,rbspx,'efw','l1',datatype+'split','YYYY',base]
    remote_root = rbsp_efw_remote_root()
    remote_path = [remote_root,rbspx,rbsp_efw_remote_sub_dirs(level='l1',datatype=datatype+'_split'),base]
    resolution = 15*60d ; sec.
    local_files = file_dailynames(file_format=join_path(local_path), trange=time_range, resolution=resolution)
    remote_files = file_dailynames(file_format=join_path(remote_path), trange=time_range, resolution=resolution)


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


;---Calbirate B field.
    in_var = prefix+datatype
    store_data, in_var, dlimits={data_att:{units:'ADC'}}
    rbsp_efw_cal_waveform, probe=probe, datatype=datatype, trange=time_range


;---Convert vector from UVW to wanted coord.
    rgb = [6,4,2]
    xyz = ['x','y','z']
    get_data, in_var, times, vec, limits=lim
    vec = cotran(vec, times, 'uvw2'+coord[0], probe=probe)
    out_var = in_var+'_'+coord[0]
    store_data, out_var, times, vec, limits={$
        ytitle:prefix+'burst_bfield!C[nT]', $
        labels:strupcase(coord)+' B'+xyz, $
        colors:rgb }

end


; Set the time and probe for loading data.
time_range = ['2013-06-10/05:57:20','2013-06-10/05:59:40']
probe = 'b'

; Load the burst data.
rbsp_efw_read_burst_efield, time_range, probe=probe, datatype='vb1', coord='mgse'
rbsp_efw_read_burst_bfield, time_range, probe=probe, datatype='mscb1', coord='mgse'

prefix = 'rbsp'+probe+'_efw_'
vars = prefix+[$

    ; The E and B fields in mgse.
    ['eb1','mscb1']+'_mgse' ]

; Plot the variables.
tplot, vars, trange=time_range
end
