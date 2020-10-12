;Load Solar Probe fields burst CDF files and put into a format that tplot can understand. 
;i.e. split up each burst into a separate tplot file with proper time formatting. 

;e.g. 
;timespan,'2018-11-03'
;psp_fld_load, type = 'dfb_dbm_scm'
;psp_fld_load, type = 'dfb_dbm_dvac'

;inputvar = 'psp_fld_l2_dfb_dbm_dvac12'
;read_fields_burst,inputvar


pro read_fields_burst,inputvar


get_data,inputvar,data=d

;   X               DOUBLE    Array[80]
;   Y               FLOAT     Array[80, 524288]
;   V               DOUBLE    Array[80, 524288]



;Loop through all n bursts (you can copy and paste this version into the prompt)
nbursts = n_elements(d.x)
for n=0,nbursts-1 do begin $


for n=0,3 do begin $
  ttmp = reform(d.x[n] + d.v[n,*]) & $
  ytmp = reform(d.y[n,*]) & $
  suffixtmp = time_string(d.x[n],prec=6) & $
  store_data,'burst_'+inputvar+'_'+suffixtmp,ttmp,ytmp
endfor


;  ;Plot example burst
;  burst_to_plot = 0
;  ;***copy and paste the below line repeatedly to plot successive bursts
;  burst_to_plot ++ & get_data,'burst_ch1_mV'+strtrim(burst_to_plot,2),data=d & timespan,d.x[0],d.x[n_elements(d.x)-1] - d.x[0],/seconds & tplot,['*mV','*nT']+strtrim(burst_to_plot,2)


end
