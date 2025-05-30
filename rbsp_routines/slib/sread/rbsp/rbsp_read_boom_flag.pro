;+
; Read RBSP boom flag, where 1 is for working boom.
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
; version=. A string to set specific version of files. By default, the
;   program finds the files of the highest version.
;-

pro rbsp_read_boom_flag_gen_file, date, probe=probe, filename=data_file, errmsg=errmsg, local_root=local_root
;---Internal, do not check inputs.
    local_root = join_path([default_local_root(),'data','rbsp'])
    rbspx = 'rbsp'+probe

;---Load Vsvy.
    max_valid_v = 200.
    secofday = constant('secofday')
    date_time_range = date+[0,secofday]
    time_range = date_time_range+[-1,1]*600.
    change_time = time_double('2016-02-28')
    if product(time_range-change_time) le 0 then begin      ; during change format.
        datatype = 'l2%vsvy-highres'
        rbsp_read_efw, [time_range[0],change_time], probe=probe, id=datatype
        get_data, 'vsvy', tmp1, tmp2
        datatype = 'l2%vsvy-highres2'
        rbsp_read_efw, [change_time,time_range[1]], probe=probe, id=datatype
        get_data, 'vsvy', tmp3, tmp4
        store_data, 'vsvy', [tmp1,tmp3], [tmp2,tmp4]
    endif else if time_range[1] lt change_time then begin   ; all before change format.
        datatype = 'l2%vsvy-highres'
        rbsp_read_efw, time_range, probe=probe, id=datatype
    endif else begin                                        ; all after change format.
        datatype = 'l2%vsvy-highres2'
        rbsp_read_efw, time_range, probe=probe, id=datatype
    endelse
    

    efield_time_step = 1d/16
    spinfit_time_step = 10d
    flag_time_step = 60.
    spin_period = rbsp_info('spin_period')
    smooth_width = round(spin_period/efield_time_step)*2    ; about 20-24 sec.
    prefix = 'rbsp'+probe+'_'
    vsvy_var = prefix+'vsvy'
    rename_var, 'vsvy', to=vsvy_var
    get_data, vsvy_var, times, vsvy
    if n_elements(vsvy) le 6 then begin
        errmsg = handle_error('No Vsvy data ...')
        return
    endif
    nboom = 4
    v_colors = sgcolor(['red','green','blue','black'])
    v_labels = 'V'+string(findgen(nboom)+1,format='(I0)')
    vsvy = vsvy[*,0:nboom-1]
    index = where(abs(vsvy) ge max_valid_v, count)
    if count ne 0 then vsvy[index] = !values.f_nan
    ; smooth to remove oscilations around perigee.
    for ii=0, nboom-1 do vsvy[*,ii] = smooth(vsvy[*,ii], smooth_width, /edge_zero, /nan)
    store_data, vsvy_var, times, vsvy, limits={$
        ytitle: '(V)', $
        colors: v_colors, $
        labels: v_labels}
    highres_times = make_bins(time_range, efield_time_step)
    interp_time, vsvy_var, highres_times
    vsvy = get_var_data(vsvy_var, limits=lim)
    index = where(finite(snorm(vsvy),/nan), count)
    fillval = !values.f_nan
    pad_time = 30.  ; sec.
    if count ne 0 then begin
        nan_times = highres_times[time_to_range(index,time_step=1)]
        nnan_time = n_elements(nan_times)*0.5
        for section_id=0,nnan_time-1 do begin
            index = lazy_where(highres_times, '[]', nan_times[section_id,*]+[-1,1]*pad_time, count=count)
            if count eq 0 then continue
            vsvy[index,*] = fillval
        endfor
        store_data, vsvy_var, highres_times, vsvy
    endif


;---Low-res version.
    vsvy_lowres_var = prefix+'vsvy_lowres'
    nhighres_time = n_elements(highres_times)
    time_index = smkarthm(0,nhighres_time,spinfit_time_step/efield_time_step,'dx')
    lowres_times = highres_times[time_index]
    vsvy = vsvy[time_index,*]
    store_data, vsvy_lowres_var, lowres_times, vsvy, limits=lim


;---Get the best-estimated Vsc.
    max_good_v = 150.  ; By inspecting monthly plots, this seems to be a good threshold for bad |V|.
    index = where(abs(vsvy) ge max_good_v, count)
    if count ne 0 then vsvy[index] = !values.f_nan
    nlowres_time = n_elements(lowres_times)
    vsc_median = fltarr(nlowres_time)+!values.f_nan
    for ii=0,nlowres_time-1 do begin
        the_vs = reform(vsvy[ii,*])
        vsc_median[ii] = median(the_vs)
    endfor
    vsc_median_var = prefix+'vsc_median'
    store_data, vsc_median_var, lowres_times, vsc_median, limits={$
        ytitle: '(V)', $
        labels: 'Vsc median'}


;---Get Vx-Vsc.
    vsvy = get_var_data(vsvy_var, times=highres_times)
    vsc = get_var_data(vsc_median_var, at=highres_times)
    for ii=0,nboom-1 do begin
        id_str = string(ii+1,format='(I0)')
        dv = vsvy[*,ii]-vsc
        dv0 = smooth(dv, smooth_width, /edge_mirror, /nan)
        dv0 = dv0[time_index]
        store_data, prefix+'dv0_'+id_str, lowres_times, dv0, limits={$
            ytitle:'(V)', $
            labels:'(V'+id_str+'-Vsc)_BG', $
            ystyle: 1, $
            yrange: [-1,1]*5}
    endfor



;---Remove SDT, eclipse.
    rbsp_read_eclipse_flag, time_range, probe=probe
    rbsp_read_sdt_flag, time_range, probe=probe

    flag_vars = prefix+['eclipse','sdt']+'_flag'
    flag_time_step = 60.    ; sec.
    get_data, flag_vars[0], flag_times
    nflag_time = n_elements(flag_times)
    other_flags = intarr(nflag_time)
    foreach flag_var, flag_vars do other_flags += get_var_data(flag_var)
    index = where(other_flags eq 1, count)
    flag_time_ranges = (count eq 0)? !null: time_to_range(flag_times[index], time_step=flag_time_step)
    nflag_time_range = n_elements(flag_time_ranges)*0.5


;---Boom flag.
    max_valid_dv0 = 5.  ; V.
    pad_time = 120.
    for ii=0,nboom-1 do begin
        id_str = string(ii+1,format='(I0)')
        dv0 = get_var_data(prefix+'dv0_'+id_str)

    ;---Tell from dv0.
        probe_flags = intarr(nlowres_time)
        index = where(abs(dv0) lt max_valid_dv0, count)
        if count ne 0 then probe_flags[index] = 1
        
    ;---Data gap if fine?
        index = where(finite(dv0,/nan), count)
        if count ne 0 then probe_flags[index] = 1
        
   
    ;---Mask times when other flags are 1.
        for jj=0, nflag_time_range-1 do begin
            index = lazy_where(lowres_times, '[]', flag_time_ranges[jj,*]+[-1,1]*pad_time, count=count)
            if count eq 0 then continue
            probe_flags[index,*] = 0
        endfor

        flag_var = prefix+'v'+id_str+'_flag'
        store_data, flag_var, lowres_times, probe_flags, limits={$
            ytitle: '(#)', $
            labels: 'V'+id_str+' flag!C  1: good', $
            ystyle: 1, $
            yrange: [0,1]+[-1,1]*0.2, $
            ytickv: [0,1], $
            panel_size: 0.4, $
            yticks: 1, $
            yminor: 0}
    endfor


;;---Fix false flag due to the data gap around day change.
;    vsc_median = get_var_data(prefix+'vsc_median')
;    index = where(finite(vsc_median,/nan), count)
;    boom_flags = intarr(nlowres_time,nboom)
;    for ii=0,nboom-1 do begin
;        boom_flags[*,ii] = get_var_data(prefix+'v'+id_str+'_flag'
;    endfor
;    
;    if count ne 0 then begin
;        nan_times = lowres_times[time_to_range(index,time_step=1)]
;        nnan_time = n_elements(nan_times)*0.5
;        max_dt = 60.
;        foreach test_time, date_time_range do begin
;            section_index = where(abs(nan_times[*,0]-test_time) le max_dt or $
;                abs(nan_times[*,0]-test_time) le max_dt, nsection)
;            if nsection eq 0 then continue
;            for ii=0,nsection-1 do begin
;                index = lazy_where(lowres_times, '[]', nan_times[section_index[ii],*], count)
;                if count eq 0 then continue
;            endfor
;            stop
;        endforeach
;    endif



;---Uniform time.
    time_index = lazy_where(lowres_times, '[]', date_time_range)
    common_times = lowres_times[time_index]
    ncommon_time = n_elements(common_times)
    
    ; Boom flag.
    boom_flags = intarr(ncommon_time,nboom)
    for ii=0, nboom-1 do begin
        id_str = string(ii+1,format='(I0)')
        boom_flags[*,ii] = (get_var_data(prefix+'v'+id_str+'_flag'))[time_index]
    endfor
    store_data, prefix+'boom_flag', common_times, boom_flags
    
    ; Other vars.
    foreach var, prefix+['vsc_median'] do begin
        data = get_var_data(var)
        store_data, var, common_times, data[time_index,*]
    endforeach




;---Write to file.
    odir = fgetpath(data_file)
    if file_test(odir,/directory) eq 0 then file_mkdir, odir
    if file_test(data_file) eq 1 then file_delete, data_file  ; overwrite old files.

    settings = dictionary($
        'title', 'RBSP flag to tell if the 4 spin-plane booms work or not',$
        'text', 'Generated by Sheng Tian at the University of Minnesota' )
    cdf_save_setting, settings, filename=data_file

    ; utsec.
    time_var = 'ut_flag'
    tdat = common_times     ; value=xxx will make xxx undefined.
    settings = dictionary($
        'unit', 'sec', $
        'time_var_type', 'unix', $
        'var_type', 'support_data')
    cdf_save_var, time_var, value=tdat, filename=data_file
    cdf_save_setting, settings, filename=data_file, varname=time_var

    ; Vsc_median.
    vsc_var = 'vsc_median'
    vsc_median = get_var_data(prefix+'vsc_median')
    settings = dictionary($
        'fieldnam', 'Vsc median', $
        'display_type', 'scalar', $
        'unit', 'V', $
        'short_name', 'V!DSC!N', $
        'depend_0', time_var)
    cdf_save_var, vsc_var, value=vsc_median, filename=data_file
    cdf_save_setting, settings, filename=data_file, varname=vsc_var

    ; boom_flag.
    boom_flag_var = 'boom_flag'
    settings = dictionary($
        'fieldnam', 'Spin plane boom flag: 1 for working', $
        'ytitle', '(#)', $
        'labels', v_labels, $
        'yrange', [-0.2,1.2], $
        'ytickv', [0,1], $
        'yticks', 1, $
        'yminor', 0, $
        'ystyle', 1, $
        'colors', v_colors, $
        'depend_0', time_var)
    boom_flags = get_var_data(prefix+'boom_flag')
    cdf_save_var, boom_flag_var, value=boom_flags, filename=data_file
    cdf_save_setting, settings, filename=data_file, varname=boom_flag_var

end

pro rbsp_read_boom_flag, time, id=datatype, probe=probe, $
    print_datatype=print_datatype, errmsg=errmsg, $
    local_files=files, file_times=file_times, version=version, $
    local_root=local_root, remote_root=remote_root

    compile_opt idl2
    on_error, 0
    errmsg = ''

;---Check inputs.
    sync_threshold = 0
    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('no probe ...')
        return
    endif
    index = where(probe eq ['a','b'])
    if index[0] eq -1 then begin
        errmsg = handle_error('invalid probe ...')
        return
    endif
    if n_elements(local_root) eq 0 then local_root = join_path([default_local_root(),'data','rbsp'])
    if n_elements(remote_root) eq 0 then remote_root = join_path([rbsp_efw_phasef_get_server()])
    if n_elements(version) eq 0 then version = 'v01'

;---Init settings.
    type_dispatch = hash()
    valid_range = (probe eq 'a')? $
        time_double(['2012-09-08','2019-10-15']): $
        time_double(['2012-09-08','2019-07-17'])
    rbspx = 'rbsp'+probe
    base_name = rbspx+'_boom_flag_%Y_%m%d_'+version+'.cdf'
    local_path = [local_root,rbspx,'boom_flag','%Y']
    remote_path = [remote_root,'boom_flag',rbspx,'%Y']


    type_dispatch['all'] = dictionary($
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
                'in_vars', ['boom_flag','vsc_median'], $
                'out_vars', rbspx+'_'+['boom_flag','vsc_median'], $
                'time_var_name', 'ut_flag', $
                'time_var_type', 'unix')))

    if keyword_set(print_datatype) then begin
        print, 'Suported data type: '
        ids = type_dispatch.keys()
        foreach id, ids do print, '  * '+id
        return
    endif

;---Dispatch patterns.
    datatype = 'all'
    request = type_dispatch[datatype]

;---Find files, read variables, and store them in memory.
    files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
        file_times=file_times, time=time, nonexist_files=nonexist_files)
    if n_elements(nonexist_files) ne 0 then begin
        foreach file, request.nonexist_files do begin
            file_time = file.file_time
            local_file = file.local_file
            rbsp_read_boom_flag_gen_file, file_time, probe=probe, filename=local_file, local_root=local_root
        endforeach
        files = prepare_files(request=request, errmsg=errmsg, local_files=files, $
            file_times=file_times, time=time, nonexist_files=nonexist_files)
    endif

;---Read data from files and save to memory.
    read_files, time, files=files, request=request

end



;probe = 'a'
;date = time_double('2014-06-19')    ; Vsvy miss data on that day.
;
;probe = 'b'
;date = time_double('2013-01-27')    ; Flag on at the beginning of the day.
;date = time_double('2014-01-09')    ; Small data gaps.
;date = time_double('2014-08-27')    ; Charging.
;date = time_double('2014-08-29')    ; Flag on in a small chunck.
;
;data_file = join_path([homedir(),'test.cdf'])
;rbsp_read_boom_flag_gen_file, date, probe=probe, filename=data_file
;stop


; Run through the whole mission.
stop
secofday = constant('secofday')
foreach probe, ['a','b'] do begin
    valid_time_range = rbsp_info('efw_l2_data_range', probe=probe)
    valid_time_range = time_double(['2016-02-27','2016-02-29'])
    dates = make_bins(valid_time_range, secofday)
    foreach date, dates do begin
        date_time_range = date+[0,secofday]
        rbsp_read_boom_flag, date_time_range, probe=probe
    endforeach
endforeach
end
