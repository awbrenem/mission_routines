
    cdf_file = ''
    probe = 'a'
    time_ranges = cdf_read_var('time_range', filename=cdf_file)
    median_srs = cdf_read_var('median_sample_rate', filename=cdf_file)
    mean_srs = cdf_read_var('mean_sample_rate', filename=cdf_file)
    bad_index = where(median_srs ne mean_srs, count)
    if count eq 0 then stop

    foreach sec_id, bad_index do begin
        time_range = time_ranges[sec_id,*]
        msg = strjoin(time_string(time_range,tformat='YYYY-MM-DD/hh:mm:ss.ffffff  '))
        msg += string(median_srs[sec_id],format='(I10)')
        msg += string(mean_srs[sec_id],format='(I10)')
        lprmsg, msg

;        rbsp_efw_read_burst, time_range, probe=probe
;        stop
    endforeach
