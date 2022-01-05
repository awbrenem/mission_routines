;TEST CODE ONLY: Attempt to subtract FIREBIRD tumble away from hires data so that I can roughly estimate the amplitude of the spikes. 
;This is mainly used to test Shumko uB identification against actual data. 
;NOTE: I've found that using 2-3 seconds as a detrend time on the FIREBIRD counts data gives generally consistent results with Shumko's uB list calculated using 0.75 sec background subtract
;Using 0.75 sec for both doesn't compare well. 
;**This tells me that my detrend way of identifying uB amplitudes is fickle. Much better to use Shumko's method.
;
;
;
;For implementation in master_conjunction_list_part3.pro
;This will be compared to Shumko's microburst id list. 

;Results: 
;--detrend hires flux at 0.75 sec. This picks up isolated uB really well, and isn't so long that it does poorly at removing suddenish increases
;   from boundary crossings or whatever else. 
;--can't detrend survey data due to long cadence (5 sec). This signal is dominated by roll, and there's really not much microburst info here. 
;Really all this gives is a guess at the total precipitation.


rbsp_efw_init
;datetime = '2015-11-17'    ;not good
;datetime = '2015-02-02'   ;Some questionable uB ids on this day
;datetime = '2015-08-27'
datetime = '2015-08-28'
;datetime = '2016-08-25'
;datetime = '2015-07-04'   ;One event abnormally high flux
;datetime = '2019-10-04'
;datetime = '2015-06-11'  ;2 sec dettime doesn't work well here. 0.5 works much better

timespan,datetime
sc = '4'
firebird_load_context_data_cdf_file,sc

sctmp = sc
firebird_load_data,sctmp

dettime = 2.  ;sec

e_channel = '0'   ;250 keV channel (lowest FB energy bin)


;For getting FIREBIRD data cadence
cal = firebird_get_calibration_counts2flux(datetime,sctmp)




split_vec,'fu'+sc+'_fb_col_hires_flux'
split_vec,'fu'+sc+'_fb_col_hires_counts'
rbsp_detrend,'fu'+sc+'_fb_col_hires_flux_'+e_channel,dettime
rbsp_detrend,'fu'+sc+'_fb_col_hires_counts_'+e_channel,dettime


ylim,['fu'+sc+'_fb_col_hires_flux_'+e_channel+'_smoothed','fu'+sc+'_fb_col_hires_flux_'+e_channel],0,0,0
store_data,'comb',data=['fu'+sc+'_fb_col_hires_flux_'+e_channel,'fu'+sc+'_fb_col_hires_flux_'+e_channel+'_smoothed']
options,'fu'+sc+'_fb_col_hires_flux_'+e_channel,'color',250
options,'comb','ytitle','Hires flux vs smoothed'




;-----------------------------------------------------------
;Find where data gaps exist and eliminate detrended data nearby. Otherwise poor fitting caused by the gap can leave large "microburst" counts
;Gaps also tend to coincide with dropouts, and so these are removed too. 

get_data,'fu'+sc+'_fb_col_hires_counts_'+e_channel,data=d
get_data,'fu'+sc+'_fb_col_hires_counts_'+e_channel+'_detrend',data=dcounts_det
get_data,'fu'+sc+'_fb_col_hires_flux_'+e_channel+'_detrend',data=dflux_det


tder = deriv(d.x)
store_data,'datagap_test',d.x,tder
options,'datagap_test','ytitle','data gaps'
ylim,'datagap_test',0.01,20,1

;tplot,['comb','datagap_test']



gaplim = 2*cal.cadence/1000.
goo = where(tder gt gaplim)
data_counts_det = dcounts_det.y 
data_flux_det = dflux_det.y

datcadence = sample_rate(dcounts_det.x,/average)

npts_pad = ceil(datcadence * dettime/2.)

;remove data during dropouts +/- the detrend time
for i=0,n_elements(goo)-1 do begin $
    boo = where((d.x ge d.x[goo[i]-npts_pad]) and (d.x le d.x[goo[i]+npts_pad])) & $
    if boo[0] ne -1 then data_counts_det[boo] = !values.f_nan
    if boo[0] ne -1 then data_flux_det[boo] = !values.f_nan
endfor


store_data,'fu4_fb_col_hires_counts_'+e_channel+'_detrend_nogaps',d.x,data_counts_det
store_data,'fu4_fb_col_hires_flux_'+e_channel+'_detrend_nogaps',d.x,data_flux_det





;----------------------------------------------------------------------------



;Compare to Shumko's microburst list
ub = firebird_load_shumko_microburst_list('4',filename='FU4_microbursts_bw=0.75sec.csv')
;ub = firebird_load_shumko_microburst_list('4',filename='FU4_microbursts_bw=2sec.csv')




tname = 'fu4_fb_col_hires_counts_'+e_channel+'_smoothed'
tname_timecorrection = 'fu4_fb_count_time_correction'

ub_flux = firebird_convert_shumko_microbursts2flux(sc, ub, tname, tname_timecorrection, float(e_channel))

;If we have Shumko microbursts, then create tplot variable. 
if finite(ub_flux.ub_times[0]) then begin 

 store_data,'ub_flux',time_double(ub_flux.ub_times),ub_flux.ub_flux
 options,'ub_flux','psym',4
 options,'ub_flux','symsize',3
 store_data,'combflux',data=['fu4_fb_col_hires_flux_'+e_channel+'_detrend_nogaps','ub_flux']

endif

tplot,['fu4_fb_count_time_correction','comb','fu4_fb_col_hires_flux_'+e_channel+'_detrend_nogaps','datagap_test','combflux']





stop



end



