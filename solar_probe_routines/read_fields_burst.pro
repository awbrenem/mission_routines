;Load Solar Probe fields burst CDF files and put into a format that tplot can understand. 
;i.e. split up each burst into a separate tplot file with proper time formatting. 


  ;   X               DOUBLE    Array[80]
  ;   Y               FLOAT     Array[80, 524288]
  ;   V               DOUBLE    Array[80, 524288]

;****
;e.g. 3D waveform data
  ;timespan,'2018-11-03'
  ;psp_fld_load, type = 'dfb_dbm_scm'
  ;psp_fld_load, type = 'dfb_dbm_dvac'

  ;inputvars = 'psp_fld_l2_dfb_dbm_dvac'+['12','34','z']
  ;inputvars = 'psp_fld_l2_dfb_dbm_scmhg'+['u','v','w']

  ;read_fields_burst,inputvars

;****
;e.g. 1D waveform data
  ;timespan,'2018-11-03'
  ;psp_fld_load, type = 'dfb_dbm_scm'
  ;psp_fld_load, type = 'dfb_dbm_dvac'

  ;inputvars = 'psp_fld_l2_dfb_dbm_dvac12'
  ;inputvars = 'psp_fld_l2_dfb_dbm_scmhgu'
  ;read_fields_burst,inputvars


pro read_fields_burst,inputvars



  ;Load 1D data
  if n_elements(inputvars) eq 1 then begin 

    get_data,inputvars,data=d

    ;Loop through all n bursts (you can copy and paste this version into the prompt)
    nbursts = n_elements(d.x)
    for n=0,nbursts-1 do begin $
      ttmp = reform(d.x[n] + d.v[n,*]) & $
      ytmp = reform(d.y[n,*]) & $
      suffixtmp = time_string(d.x[n],prec=6) & $
      store_data,'burst_'+inputvars+'_'+suffixtmp,ttmp,ytmp
    endfor

  endif



  ;Load 3D data
  if n_elements(inputvars) eq 3 then begin 


    ;determine tplot variable pre-suffix (before time is tacked on)
    get_data,inputvars[0],data=d
    tst = strmid(inputvars[0],0,1,/reverse_offset)
    if tst eq '2' then suffix = 'psp_fld_l2_dfb_dbm_dvac'
    if tst eq 'u' then suffix = 'psp_fld_l2_dfb_dbm_scm'


    get_data,inputvars[0],data=d1
    get_data,inputvars[1],data=d2
    get_data,inputvars[2],data=d3


    ;Loop through all n bursts (you can copy and paste this version into the prompt)
    nbursts = n_elements(d1.x)
    for n=0,nbursts-1 do begin $
      ttmp = reform(d1.x[n] + d1.v[n,*]) & $
      ytmp1 = reform(d1.y[n,*]) & $
      ytmp2 = reform(d2.y[n,*]) & $
      ytmp3 = reform(d3.y[n,*]) & $
      ytmp = [[ytmp1],[ytmp2],[ytmp3]] & $
      suffixtmp = time_string(d1.x[n],prec=6) & $
      store_data,'burst_'+suffix+'_'+suffixtmp,ttmp,ytmp
    endfor


  endif

  ;  ;Plot example burst
  ;  burst_to_plot = 0
  ;  ;***copy and paste the below line repeatedly to plot successive bursts
  ;  burst_to_plot ++ & get_data,'burst_ch1_mV'+strtrim(burst_to_plot,2),data=d & timespan,d.x[0],d.x[n_elements(d.x)-1] - d.x[0],/seconds & tplot,['*mV','*nT']+strtrim(burst_to_plot,2)


end
