;+
; Read EFW burst data. Save data in rbspx_efw_<xxx>.
;
; datatype can be 'vb1','vb2','mscb1','mscb2'.
;-

pro rbsp_efw_read_burst, tr, probe=probe, datatype=datatype, trange=trange, $
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
    timespan, time_range[0], total(time_range*[-1,1]), /seconds
    if n_elements(datatype) eq 0 then datatype = 'vb1'
    data_types = ['vb1','vb2','mscb1','mscb2','eb1','eb2']
    index = where(data_types eq datatype[0], count)
    if count eq 0 then begin
        dprint, 'Invalid datatype: '+datatype[0]+' ...', verbose=vb
        return
    endif


    rbspx = 'rbsp'+probe
    base = rbspx+'_efw_l1_'+datatype+'-split_YYYYMMDDThhmm00_'+version+'.cdf'
    local_root = !rbsp_efw.local_data_dir
    local_path = [local_root,rbspx,'efw','l1',datatype+'-split','YYYY',base]
    remote_root = rbsp_efw_remote_root()
    remote_path = [remote_root,rbspx,rbsp_efw_remote_sub_dirs(level='l1',datatype=datatype+'-split_v03'),base]
    resolution = 15*60d ; sec.
    local_files = file_dailynames(file_format=join_path(local_path), trange=time_range, resolution=resolution)
    remote_files = file_dailynames(file_format=join_path(remote_path), trange=time_range, resolution=resolution)

    if strmid(datatype,0,/reverse) eq '2' then begin
        base = rbspx+'_l1_'+datatype+'_YYYYMMDD_'+version+'.cdf'
        local_root = !rbsp_efw.local_data_dir
        local_path = [local_root,rbspx,'efw','l1',datatype,'YYYY',base]
        remote_root = 'http://themis.ssl.berkeley.edu/data/rbsp/rbspa/l1/'
        remote_path = [datatype,'YYYY',base]
        local_files = file_dailynames(file_format=join_path(local_path), trange=time_range)
        remote_files = remote_root+file_dailynames(file_format=join_path(remote_path), trange=time_range)
    endif


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

    v_var = prefix+datatype
    store_data, v_var, dlimits={data_att:{units:'ADC'}}
    get_data, v_var, data=dd
    if size(dd,/type) ne 8 then message, 'No data ...'
    rbsp_efw_cal_waveform, probe=probe, datatype=datatype, trange=minmax(dd.x)


end


; Set the time and probe for loading data.
time_range = ['2013-06-10/05:57:20','2013-06-10/05:59:40']
probe = 'b'

; Load the spinfit data.
rbsp_efw_read_burst, time_range, probe=probe, datatype='vb1'

prefix = 'rbsp'+probe+'_efw_'
vars = prefix+[ 'vb1' ]

; Plot the variables.
tplot, vars, trange=time_range

end
