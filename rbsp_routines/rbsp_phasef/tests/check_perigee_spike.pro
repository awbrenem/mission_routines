;+
; Check spikes around the day change.
;-

;---Settings.
    time_range = time_double(['2013-02-22','2013-02-24'])
    time_range = time_double(['2014-10-24','2014-10-26'])
    time_range = time_double(['2014-09-05','2014-09-07'])
    probe = 'a'


;---Derived quantities.
    prefix = 'rbsp'+probe+'_'
    common_time_step = 10.
    common_times = make_bins(time_range, common_time_step)
    xyz = constant('xyz')
    uvw = constant('uvw')
    the_time_range = mean(time_range)+[-1,1]*60


;---Load data.
;    rbsp_read_q_uvw2gse, time_range, probe=probe
    rbsp_efw_phasef_read_wobble_free_var, time_range, probe=probe
    timespan, time_range[0], total(time_range*[-1,1]), /seconds
    ;rbsp_load_efw_waveform, probe=probe, datatype='esvy', type='cal', coord='uvw', /noclean, trange=time_range
    ;stop
    rbsp_load_efw_waveform, probe=probe, datatype='esvy', type='raw', coord='uvw', /noclean, trange=time_range
    evar_old = prefix+'efw_esvy'
    options, evar_old, 'colors', constant('rgb')
    sgopen, 0, xsize=8, ysize=5
    tplot, prefix+['efw_esvy','e_mgse'], trange=the_time_range


end