;+
; Generate E after wake correction spinfit for 12 and 34 pairs, based on the calibrated E UVW.
;-

pro rbsp_phasef_read_e_wake_spinfit_gen_file, time, probe=probe, filename=file

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


;---Do spinfit.
    ; Get the time range.
    secofday = 86400d
    date = time-(time mod secofday)
    tr = date+[0,secofday]
    rbspx = 'rbsp'+probe
    prefix = 'rbsp'+probe+'_'
    timespan, tr[0], total(tr*[-1,1]), /seconds

    rbsp_load_spice_cdf_file, probe

    ; Restore the calibrated E UVW data.
    rbsp_efw_phasef_read_wake_flag, tr, probe=probe, id='euv'
    spin_axis_var = rbspx+'_spinaxis_direction_gse'
    e_uvw_var = rbspx+'_e_uvw'
    ; Add UVW to dlim to tell spinfit the coord.
    get_data, prefix+'eu_fixed', uts, eu
    get_data, prefix+'ev_fixed', uts, ev
    nut = n_elements(uts)
    e_uvw = [[eu],[ev],[fltarr(nut)]]
    store_data, e_uvw_var, uts, e_uvw
    data_att = {coord_sys:'uvw'}
    dlim = {data_att:data_att}
    store_data, e_uvw_var, dlimits=dlim

    spinfit_var = e_uvw_var+'_spinfit'
    pairs = ['12','34']
    foreach pair, pairs, pair_id do begin
        rbsp_spinfit, e_uvw_var, plane_dim=pair_id
        dsc_var = prefix+'e_wake_spinfit_v'+pair+'_dsc'
        the_var = prefix+'e_wake_spinfit_mgse_v'+pair

        ; Transform the spinfit data from DSC to MGSE (AARON'S UPDATED VERSION which is very fast and gives same result as old method)
        copy_data, spinfit_var, dsc_var
        ;tinterpol_mxn, spin_axis_var, dsc_var,/quadratic. This is done in rbsp_efw_dsc_to_mgse.
        rbsp_efw_dsc_to_mgse, probe, dsc_var, spin_axis_var
        mgse_var = dsc_var+'_mgse'
        copy_data, mgse_var, the_var
    endforeach


;---Save data.
    path = fgetpath(file)
    if file_test(path,/directory) eq 0 then file_mkdir, path
    data_file = file
    if file_test(data_file) eq 1 then file_delete, data_file  ; overwrite old files.

    ginfo = dictionary($
        'TITLE', 'RBSP EFW E after wake correction spinfit in the corotation frame', $
        'TEXT', 'Generated by Sheng Tian at the University of Minnesota, adopted from rbsp_efw_spinfit_vxb_crib' )
    cdf_save_setting, ginfo, filename=file
    save_vars = prefix+'e_wake_spinfit_mgse_v'+pairs
    foreach save_var, save_vars do begin
        get_data, save_var, times, e_uvw
        store_data, save_var, times, float(e_uvw), limits={units:'mV/m', coord:'mgse'}
    endforeach
    stplot2cdf, save_vars, istp=1, filename=file, time_var='epoch'

end



;time = time_double('2014-08-01')
;probe = 'a'
;file = join_path([homedir(),'test.cdf'])
;if file_test(file) eq 1 then file_delete, file
;tic
;rbsp_phasef_read_e_wake_spinfit_gen_file, time, probe=probe, filename=file
;toc
;stop

probes = ['b']
root_dir = join_path([default_local_root(),'data','rbsp'])
foreach probe, probes do begin
    prefix = 'rbsp'+probe+'_'
    rbspx = 'rbsp'+probe
    time_range = (probe eq 'a')? time_double(['2012-09-08','2019-10-14']): time_double(['2012-09-08','2019-07-16'])
    days = make_bins(time_range, constant('secofday'))
    foreach day, days do begin
        str_year = time_string(day,tformat='YYYY')
        path = join_path([root_dir,rbspx,'e_wake_spinfit_mgse',str_year])
        base = prefix+'efw_e_wake_spinfit_mgse_'+time_string(day,tformat='YYYY_MMDD')+'_v01.cdf'
        file = join_path([path,base])
        rbsp_phasef_read_e_wake_spinfit_gen_file, day, probe=probe, filename=file
    endforeach
endforeach


stop



end