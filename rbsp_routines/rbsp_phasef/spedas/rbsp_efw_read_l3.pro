;+
; Read EFW L3 data.
;-

pro rbsp_efw_read_l3, tr, probe=probe, datatype=datatype, trange=trange, $
    level=level, verbose=verbose, downloadonly=downloadonly, $
     cdf_data=cdf_data,get_support_data=get_support_data, $
     tplotnames=tns, make_multi_tplotvar=make_multi_tplotvar, $
     varformat=varformat, valid_names = valid_names, files=files, $
     type=type, _extra = _extra

    rbsp_efw_init
    vb = keyword_set(verbose) ? verbose : 0
    vb = vb > !rbsp_efw.verbose

    if n_elements(probe) eq 0 then probe = 'a'
    if n_elements(version) eq 0 then version = 'v04'
    if n_elements(trange) ne 0 then time_range = trange
    if n_elements(tr) ne 0 then time_range = tr
    if n_elements(time_range) eq 0 then time_range = timerange()
    if size(time_range[0],/type) eq 7 then time_range = time_double(time_range)
    

    rbspx = 'rbsp'+probe
    base = rbspx+'_efw-l3_YYYYMMDD_'+version+'.cdf'
    local_root = !rbsp_efw.local_data_dir
    local_path = [local_root,rbspx,'efw','l3','YYYY',base]
    remote_root = rbsp_efw_remote_root()
    remote_path = [remote_root,rbspx,rbsp_efw_remote_sub_dirs(level='l3'),base]
    local_files = file_dailynames(file_format=join_path(local_path), trange=time_range)
    remote_files = file_dailynames(file_format=join_path(remote_path), trange=time_range)


    nfile = n_elements(local_files)
    for file_id=0,nfile-1 do begin
        tmp = spd_download_file(url=remote_files[file_id], filename=local_files[file_id])
;        remote_path = file_dirname(remote_files[file_id])
;        remote_file = file_basename(remote_files[file_id])
;        local_path = file_dirname(local_files[file_id])
;        local_file = file_basename(local_files[file_id])
;        tmp = spd_download(remote_file=remote_file, remote_path=remote_path, $
;            local_file=local_file, local_path=local_path)
    endfor



    suffix = ''
    prefix = rbspx+'_efw_'
    cdf2tplot, file=local_files, all=0, prefix=prefix, suffix=suffix, verbose=vb, $
        tplotnames=tns, /convert_int1_to_int2, get_support_data=0

end



; Set the time and probe for loading data.
time_range = ['2013-01-01','2013-01-03']
probe = 'a'

; Load the spinfit data.
rbsp_efw_read_l3, time_range, probe=probe

prefix = 'rbsp'+probe+'_efw_'
vars = prefix+[$

    ; The spinfit E field with E_spinaxis = 0.
    'efield_in_corotation_frame_spinfit_mgse', $

    ; The spinfit E field with E_spinaxis calculated from E dot B = 0,
    ; when B is away from the spin plane by >15 deg.
    'efield_in_corotation_frame_spinfit_edotb_mgse', $

    ; The spacecraft potential.
    'spacecraft_potential', $

    ; The EFW density calibrated according to the upper-hybrid line.
    'density', $

    ; The ephemeris data.
    'position_gse', 'velocity_gse', 'mlt', 'mlat', 'lshell', 'orbit_num' ]

; Plot the variables.
tplot, vars, trange=time_range
end
