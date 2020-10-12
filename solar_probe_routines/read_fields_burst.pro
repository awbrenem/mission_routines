;Load Solar Probe fields burst CDF files from Keith
;e.g. spp_fld_l1_tds_wf_20190405_v00.cdf


rbsp_efw_init


;---------------------------------------------------
;load up file
path = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/solar_probe_routines/'
fn = 'spp_fld_l1_tds_wf_20190405_v00.cdf'
cdf2tplot,path+fn



;-----------------------------------------------------
;Plot locations of bursts along with peak values
options,'SPP_FLD_TDS_WF_Burst_ID','psym',-4
options,'SPP_FLD_TDS_WF_Burst_ID','symsize',1.5
;Plot Peak values in each channel.
options,['SPP_FLD_TDS_WF_Peak_Ch?_mV','SPP_FLD_TDS_WF_Peak_Ch?_nT','SPP_FLD_TDS_WF_Burst_Quality'],'psym',-5
tplot,['SPP_FLD_TDS_WF_Burst_ID','SPP_FLD_TDS_WF_Burst_Quality','SPP_FLD_TDS_WF_Peak_Ch?_mV','SPP_FLD_TDS_WF_Peak_Ch?_nT']
;-----------------------------------------------------


;Put each burst into a format that tplot can handle
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch1_mV',data=e1
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch2_mV',data=e2
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch3_mV',data=e3
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch4_mV',data=e4
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch5_mV',data=e5
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch4_nT',data=m4
get_data,'SPP_FLD_TDS_WF_Burst_Time_Series_Ch5_nT',data=m5


;Loop through all n bursts (you can copy and paste this version into the prompt)
nbursts = n_elements(e1.x)
for n=0,nbursts-1 do begin $
  goo = where(e1.v[n,*] ge 1d9) & $
  store_data,'burst_ch1_mV'+strtrim(n,2),reform(e1.v[n,goo]),reform(e1.y[n,goo]) & $
  goo = where(e2.v[n,*] ge 1d9) & $
  store_data,'burst_ch2_mV'+strtrim(n,2),reform(e2.v[n,goo]),reform(e2.y[n,goo]) & $
  goo = where(e3.v[n,*] ge 1d9) & $
  store_data,'burst_ch3_mV'+strtrim(n,2),reform(e3.v[n,goo]),reform(e3.y[n,goo]) & $
  goo = where(e4.v[n,*] ge 1d9) & $
  store_data,'burst_ch4_mV'+strtrim(n,2),reform(e4.v[n,goo]),reform(e4.y[n,goo]) & $
  goo = where(e5.v[n,*] ge 1d9) & $
  store_data,'burst_ch5_mV'+strtrim(n,2),reform(e5.v[n,goo]),reform(e5.y[n,goo]) & $
  goo = where(m4.v[n,*] ge 1d9) & $
  store_data,'burst_ch4_nT'+strtrim(n,2),reform(m4.v[n,goo]),reform(m4.y[n,goo]) & $
  goo = where(m5.v[n,*] ge 1d9) & $
  store_data,'burst_ch5_nT'+strtrim(n,2),reform(m5.v[n,goo]),reform(m5.y[n,goo])
endfor

  ;Plot example burst
  burst_to_plot = 0
  ;***copy and paste the below line repeatedly to plot successive bursts
  burst_to_plot ++ & get_data,'burst_ch1_mV'+strtrim(burst_to_plot,2),data=d & timespan,d.x[0],d.x[n_elements(d.x)-1] - d.x[0],/seconds & tplot,['*mV','*nT']+strtrim(burst_to_plot,2)


stop
end
