;+
; Save position variables.
;-

pro rbsp_efw_phasef_save_pos_var_to_file, date, probe=probe, filename=file, $
    errmsg=errmsg, log_file=log_file, saved_vars=new_vars


    errmsg = ''

;---Check input.
    if n_elements(file) eq 0 then begin
        errmsg = 'cdf file is not set ...'
        lprmsg, errmsg, log_file
        return
    endif
    if file_test(file) eq 0 then begin
        errmsg = 'file does not exist ...'
        lprmsg, errmsg, log_file
        return
    endif

    if n_elements(probe) eq 0 then begin
        errmsg = 'No input probe ...'
        lprmsg, errmsg, log_file
        return
    endif
    if probe ne 'a' and probe ne 'b' then begin
        errmsg = 'Invalid probe: '+probe+' ...'
        lprmsg, errmsg, log_file
        return
    endif
    prefix = 'rbsp'+probe+'_'
    rbspx = 'rbsp'+probe

    data_type = 'pos_var'
    valid_range = rbsp_efw_phasef_get_valid_range(data_type, probe=probe)
    if n_elements(date) eq 0 then begin
        errmsg = 'No input date ...'
        lprmsg, errmsg, log_file
        return
    endif
    if size(date,/type) eq 7 then date = time_double(date)
    if product(date-valid_range) gt 0 then begin
        errmsg = 'Input date: '+time_string(date,tformat='YYYY-MM-DD')+' is out of valid range ...'
        lprmsg, errmsg, log_file
        return
    endif


;---Check data availability.
    secofday = 86400d
    time_range = date+[0,secofday]
    timespan, date, secofday, /second

    ; Needs rbspx_r_gse, rbspx_v_gse, rbspx_mlt, rbspx_mlat, rbspx_lshell.
    ;rbsp_efw_phasef_read_pos_var, date, probe=probe, errmsg=errmsg, log_file=log_file
    tplot_vars = prefix+['r_gse','v_gse','mlt','mlat','lshell']
    foreach var, tplot_vars do begin
        if check_if_update(var) then begin
            errmsg = 'Var does not exist: '+var+' ...'
            lprmsg, errmsg, log_file
            return
        endif
    endforeach


;---Save to CDF.
    new_vars = ['pos_gse','vel_gse','mlt','mlat','lshell']
    foreach new_var, new_vars, var_id do begin
        tplot_var = tplot_vars[var_id]
        new_var = new_vars[var_id]
        get_data, tplot_var, times, data
        ntime = n_elements(times)
        
        if ~cdf_has_var('pos_gse', filename=file) then begin
            cdf_rename_var, 'position_gse', to='pos_gse', filename=file
        endif
        if ~cdf_has_var('vel_gse', filename=file) then begin
            cdf_rename_var, 'velocity_gse', to='vel_gse', filename=file
        endif

        ; Check if need to write epoch.
        vatts = cdf_read_setting(new_var, filename=file)
        time_var = vatts['DEPEND_0']
        epochs = cdf_read_var(time_var, filename=file)
        if n_elements(epochs) eq 0 then begin
            ;epochs = tplot_time_to_epoch(times, epoch16=1)
            epochs = convert_time(times, from='unix', to='epoch16')
            cdf_save_data, time_var, value=epochs, filename=file
        endif
        if n_elements(epochs) ne n_elements(times) then message, 'Inconsistent epoch and data ...'

        ; Save data.
        cdf_save_data, new_var, value=transpose(data), filename=file

        ; Update labeling.
        wanted_label = phasef_get_labeling(new_var)
        foreach the_key, ['VAR_NOTES','UNITS','LABL_PTR_1','CATDESC'] do begin
            if ~wanted_label.haskey(the_key) then continue
            vatts[the_key] = wanted_label[the_key]
        endforeach
        cdf_save_setting, vatts, filename=file, varname=new_var

        the_key = 'labels'
        if wanted_label.haskey(the_key) then begin
            label_var = vatts['LABL_PTR_1']
            if cdf_has_var(label_var, filename=file) then begin
                label_vatts = cdf_read_setting(label_var, filename=file)
            endif else label_vatts = dictionary()
            cdf_save_var, label_var, filename=file, $
                value=transpose(wanted_label[the_key])
            if n_elements(label_vatts) ne 0 then begin
                cdf_save_setting, label_vatts, filename=file, varname=label_var
            endif
        endif
    endforeach

end
