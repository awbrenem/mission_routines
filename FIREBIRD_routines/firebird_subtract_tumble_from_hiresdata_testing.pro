;TEST CODE ONLY: Attempt to subtract FIREBIRD tumble away from hires data so that I can roughly estimate the amplitude of the spikes. 
;For implementation in master_conjunction_list_part3.pro

;Results: 
;--detrend hires flux at 0.75 sec. This picks up isolated uB really well, and isn't so long that it does poorly at removing suddenish increases
;   from boundary crossings or whatever else. 
;--can't detrend survey data due to long cadence (5 sec). This signal is dominated by roll, and there's really not much microburst info here. 
;Really all this gives is a guess at the total precipitation.


rbsp_efw_init
;datetime = '2015-11-18'
;datetime = '2015-12-12'
datetime = '2019-10-04'
;datetime = '2015-06-11'  ;2 sec dettime doesn't work well here. 0.5 works much better

timespan,datetime
sc = '4'
firebird_load_context_data_cdf_file,sc

sctmp = sc
firebird_load_data,sctmp


dettime = 0.75  ;sec

sctmp = sc
cal = firebird_get_calibration_counts2flux(datetime,sctmp)
;chn = strtrim(cal.CHANNEL_USED_FOR_SURVEY_CALIBRATION - 1,2)


split_vec,'fu'+sc+'_fb_col_hires_flux'
tplot,['fu'+sc+'_fb_col_hires_flux_0','flux_context_FU4']

rbsp_detrend,'fu'+sc+'_fb_col_hires_flux_0',dettime
rbsp_detrend,'flux_context_FU'+sc,15
ylim,['fu'+sc+'_fb_col_hires_flux_0_smoothed','fu'+sc+'_fb_col_hires_flux_0'],0,0,0
ylim,['flux_context_FU'+sc,'flux_context_FU'+sc+'_detrend'],0,0,0
store_data,'comb',data=['fu'+sc+'_fb_col_hires_flux_0','fu'+sc+'_fb_col_hires_flux_0_smoothed']
options,'comb','colors',[0,250]
options,'comb','ytitle','Hires flux vs smoothed'

store_data,'comb_hires+survey',data=['flux_context_FU'+sc,'fu'+sc+'_fb_col_hires_flux_0']
options,'comb_hires+survey','colors',[250,0]

tplot,['comb_hires+survey','comb','fu'+sc+'_fb_col_hires_flux_0_detrend','flux_context_FU'+sc,'flux_context_FU'+sc+'_detrend']


stop




;Find where data gaps exist and eliminate detrended data nearby. Otherwise poor fitting caused by the gap can leave large "microburst" counts
;Gaps also tend to coincide with dropouts, and so these are removed too. 

get_data,'fu'+sc+'_fb_col_hires_flux_0',data=d

get_data,'fu'+sc+'_fb_col_hires_flux_0_detrend',data=d_det
tder = deriv(d.x)
store_data,'datagap_test',d.x,tder
options,'datagap_test','ytitle','data gaps'
ylim,'datagap_test',0.01,20,1

tplot,['comb','datagap_test']



gaplim = 2*cal.cadence/1000.
goo = where(tder gt gaplim)
data = d_det.y 
datcadence = sample_rate(d_det.x,/average)

npts_pad = ceil(datcadence * dettime/2.)

;remove data during dropouts +/- the detrend time (2 sec)
for i=0,n_elements(goo)-1 do begin $
    boo = where((d.x ge d.x[goo[i]-npts_pad]) and (d.x le d.x[goo[i]+npts_pad])) & $
    if boo[0] ne -1 then data[boo] = !values.f_nan
endfor

store_data,'fb_corr',d.x,data
store_data,'detcomb',data=['fu'+sc+'_fb_col_hires_flux_0_detrend','fb_corr'] & options,'detcomb','colors',[0,250]
options,'detcomb','ytitle','full detrended data vs!Cgap-reduced version'




tplot,['comb','datagap_test','detcomb']


;Compare to Shumko's microburst list 

ub = load_firebird_microburst_list('4')
;timed = time_double(ub.time)
;goo = where((strmid(ub.time,0,10) eq datetime))

store_data,'ub_shumko',data={x:ub.time,y:ub.flux_ch1}

timebar,ub.time


stop

end



