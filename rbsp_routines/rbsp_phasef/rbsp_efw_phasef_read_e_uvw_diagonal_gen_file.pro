;+
; Adopted from rbsp_efw_spinfit_vxb_subtract_crib
;
;-

pro rbsp_efw_phasef_read_e_uvw_diagonal_gen_file, time, probe=probe, filename=file

;---Check inputs.
    if n_elements(file) eq 0 then begin
        errmsg = handle_error('No output file ...')
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = handle_error('No input probe ...')
        return
    endif

    if n_elements(time) eq 0 then begin
        errmsg = handle_error('No input time ...')
        return
    endif


;---Constants and settings.
    secofday = 86400d
    errmsg = ''


;---Load E UVW.
    rbspx = 'rbsp'+probe
    prefix = 'rbsp'+probe+'_'
    date = time[0]-(time[0] mod secofday)
    time_range = date+[0,secofday]
    tr = time_range+[-1,1]*60d
    timespan, tr[0], total(tr*[-1,1]), /seconds
    rbsp_load_efw_waveform, probe=probe, datatype='vsvy', coord='uvw', noclean=1


;---Manually rule out some times.
    mask_list = list()
    mask_list.add, dictionary($
        'probe', 'a', $
        'time_range', time_double(['2017-04-13/23:00','2017-04-14/00:01']))
    mask_list.add, dictionary($
        'probe', 'b', $
        'time_range', time_double(['2015-06-12/10:00','2015-06-12/10:40']))

    l1_efw_var = prefix+'efw_vsvy'
    get_data, l1_efw_var, times, vsvy
    foreach mask, mask_list do begin
        if mask.probe ne probe then continue
        index = lazy_where(times, '[]', mask.time_range, count=count)
        if count ne 0 then times[index] = !values.f_nan
    endforeach
    index = where(finite(times))
    times = times[index]
    vsvy = vsvy[index,*]
    store_data, l1_efw_var, times, vsvy


;---Implement the time correction, if necessary.
    rbsp_efw_read_l1_time_tag_correction, probe=probe
    get_data, prefix+'l1_time_tag_correction', start_times, time_ranges, corrections
    nsection = n_elements(corrections)
    if n_elements(time_ranges) le 1 then nsection = 0
    l1_efw_var = 'rbsp'+probe+'_efw_vsvy'
    get_data, l1_efw_var, times, vsvy
    var_updated = 0
    for ii=0, nsection-1 do begin
        tmp = where(times ge time_ranges[ii,0] and times le time_ranges[ii,1], count)
        if count eq 0 then continue
        var_updated = 1
        time_step = sdatarate(times)
        dtimes = times[1:-1]-times[0:-2]
        bad_index = where(abs(dtimes) ge 0.5, count)
        if count eq 0 then bad_index = !null    ; Can have no bad index but the whole day is shifted.
        if min(times) ge time_ranges[ii,0] then i0 = 0 else begin
            foreach index, bad_index do begin
                if round(dtimes[index]) ne -round(corrections[ii]) then continue
                if abs(times[index+1]-time_ranges[ii,0]) ge time_step then continue
                i0 = index+1
            endforeach
        endelse
        if max(times) le time_ranges[ii,1] then i1 = n_elements(times) else begin
            foreach index, bad_index do begin
                if round(dtimes[index]) ne round(corrections[ii]) then continue
                if abs(times[index+1]-time_ranges[ii,1]) ge time_step then continue
                i1 = index+1
            endforeach
        endelse
        times[i0:i1-1] += corrections[ii]
    endfor
    if var_updated then store_data, l1_efw_var, times, vsvy


;---Leap second corrections.
    rbsp_efw_read_l1_time_tag_leap_second, probe=probe
    get_data, prefix+'l1_time_tag_leap_second', start_times
    nsection = n_elements(start_times)
    get_data, l1_efw_var, times, vsvy
    var_updated = 0
    for ii=0, nsection-1 do begin
        tmp = where(time[0] le start_times[ii] and times[-1] ge start_times[ii], count)
        if count eq 0 then continue
        var_updated = 1
        dtimes = times[1:-1]-times[0:-2]
        bad_index = where(abs(dtimes) ge 0.5, count)
        for jj=0, count-1 do begin
            i0 = bad_index[jj]
            tmp = times[0:i0]
            index = where(tmp ge times[i0+1], count)
            if count eq 0 then continue
            tmp[index] = !values.d_nan
            times[0:i0] = tmp
        endfor
        index = where(finite(times), count)
        if count eq 0 then stop    ; Something is wrong, must have some finite data.
        times = times[index]
        vsvy = vsvy[index,*]
    endfor
    if var_updated then store_data, l1_efw_var, times, vsvy

    dtimes = round(times[1:-1]-times[0:-2])
    bad_index = where(abs(dtimes) eq 1, count)
    if count ne 0 then stop


;---Get rid of non-monotonic times, which sometimes show up.
    get_data, l1_efw_var, times, vsvy
    index = uniq(times, sort(times))
    times = times[index]
    vsvy = vsvy[index,*]
    store_data, l1_efw_var, times, vsvy


;---Construct E combo.
    ; Get boom length.
    cp0 = rbsp_efw_get_cal_params(tr[0])
    cp = (probe eq 'a')? cp0.a: cp0.b
    boom_length = cp.boom_length
    boom_shorting_factor = cp.boom_shorting_factor

    boom_length_adj = sqrt(2)*boom_length[0]*0.5
    pairs = ['13','14','23','24']
    foreach pair, pairs do begin
        idx = [strmid(pair,0,1),strmid(pair,1,1)]-1
        edata = (vsvy[*,idx[0]]-vsvy[*,idx[1]])*(1000./boom_length_adj)
        e_var = prefix+'efw_esvy_'+pair
        store_data, e_var, times, edata, $
            limits={units:'mV/m', coord:'uvw', labels:pair}
    endforeach



;---Remove DC offset.
    spin_period = 11d   ; the rough number works fine, no need to get the accurate number
    dt = median(times[1:-1]-times[0:-2])
    width = spin_period/dt
    e_vars = prefix+'efw_esvy_'+pairs
    foreach e_var, e_vars do begin
        get_data, e_var, times, edata
        offset1 = smooth(edata, width, /nan, /edge_zero)
        offset2 = smooth(offset1, width, /nan, /edge_zero)
        edata -= offset2
        store_data, e_var, times, float(edata)
    endforeach
    ; Trim data to the wanted time range.
    index = where(times ge time_range[0] and times lt time_range[1], count)
    if count ne 0 then begin
        foreach e_var, e_vars do begin
            get_data, e_var, times, edata
            store_data, e_var, times[index], edata[index]
        endforeach
    endif


;---Save data.
    path = fgetpath(file)
    if file_test(path,/directory) eq 0 then file_mkdir, path
    data_file = file
    if file_test(data_file) eq 1 then file_delete, data_file  ; overwrite old files.

    ginfo = dictionary($
        'TITLE', 'RBSP EFW calibrated E combo in the spacecraft frame', $
        'TEXT', 'Generated by Sheng Tian at the University of Minnesota, adopted from rbsp_efw_spinfit_vxb_crib' )
    cdf_save_setting, ginfo, filename=file
    save_vars = e_vars
    stplot2cdf, save_vars, istp=1, filename=file, time_var='epoch'

end


;stop
;probes = ['a','b']
;root_dir = join_path([default_local_root(),'data','rbsp'])
;foreach probe, probes do begin
;    prefix = 'rbsp'+probe+'_'
;    rbspx = 'rbsp'+probe
;    time_range = (probe eq 'a')? time_double(['2012-09-08','2019-10-14']): time_double(['2012-09-08','2019-07-16'])
;    time_range = time_range>time_double('2016-01-01')
;    days = make_bins(time_range, constant('secofday'))
;    foreach day, days do begin
;        str_year = time_string(day,tformat='YYYY')
;        path = join_path([root_dir,rbspx,'e_uvw_diagonal',str_year])
;        base = prefix+'efw_e_combo_'+time_string(day,tformat='YYYY_MMDD')+'_v01.cdf'
;        file = join_path([path,base])
;        rbsp_efw_phasef_read_e_uvw_diagonal_gen_file, day, probe=probe, filename=file
;    endforeach
;endforeach


probe = 'a'
time = time_double(['2016-01-01'])
;time = time_double(['2013-01-01'])
probe = 'b'
time = time_double('2015-07-01')
file = join_path([homedir(),'test.cdf'])
if file_test(file) eq 1 then file_delete, file
rbsp_efw_phasef_read_e_uvw_diagonal_gen_file, time, probe=probe, filename=file
end
