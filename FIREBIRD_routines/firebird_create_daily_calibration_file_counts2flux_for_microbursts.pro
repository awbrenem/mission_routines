;Create a text file with a daily conversion factor to fix the simple calibration to Shumko microburst fluxes. 
;Mike's are based on a simple counts -> flux conversion. I now have a much better conversion via firebird_get_calibration_counts2flux.pro
;
;This code needs to only be run once. Once the text file is created, then I can just load it. 
;Currently this is done in the routine load_firebird_microburst_list.pro  
;
;
datetime = '2016-08-25'




timespan,datetime
sc = '4'
firebird_load_context_data_cdf_file,sc

sctmp = sc
firebird_load_data,sctmp


dettime = 0.75  ;sec


split_vec,'fu'+sc+'_fb_col_hires_flux'



sctmp = sc
cal = firebird_get_calibration_counts2flux(datetime,sctmp)
;chn = strtrim(cal.CHANNEL_USED_FOR_SURVEY_CALIBRATION - 1,2)

flux_conv = 1/((cal.cadence/1000.) * (cal.energy_range_collimated[0,1] - cal.energy_range_collimated[0,0]) * cal.g_factor_collimated[0])
; flux = counts/cadence/energy_width/geometric_factor



;Compare the flux data from firebird_load_data.pro to my own calibration here based on firebird_get_calibration_counts2flux.pro

split_vec,'fu4_fb_col_hires_counts'
get_data,'fu4_fb_col_hires_counts_0',data=d
store_data,'fu4_fb_col_hires_flux_0_aaron',data={x:d.x,y:d.y*flux_conv}

tplot,['fu4_fb_col_hires_counts_0','fu4_fb_col_hires_flux_0','fu4_fb_col_hires_flux_0_aaron']


;***tmpp
copy_data,'fu4_fb_col_hires_flux_0','fu4_fb_col_hires_flux_0_backup'
copy_data,'fu4_fb_col_hires_flux_0_aaron','fu4_fb_col_hires_flux_0'


tplot,['fu'+sc+'_fb_col_hires_flux_0','flux_context_FU4']

rbsp_detrend,'fu'+sc+'_fb_col_hires_flux_0',dettime
rbsp_detrend,'fu'+sc+'_fb_col_hires_flux_0_backup',dettime
rbsp_detrend,'flux_context_FU'+sc,15
ylim,['fu'+sc+'_fb_col_hires_flux_0_smoothed','fu'+sc+'_fb_col_hires_flux_0'],0,0,0
ylim,['flux_context_FU'+sc,'flux_context_FU'+sc+'_detrend'],0,0,0
store_data,'comb',data=['fu'+sc+'_fb_col_hires_flux_0','fu'+sc+'_fb_col_hires_flux_0_smoothed']
options,'comb','colors',[0,250]
options,'comb','ytitle','Hires flux vs smoothed'

store_data,'comb_hires+survey',data=['flux_context_FU'+sc,'fu'+sc+'_fb_col_hires_flux_0']
options,'comb_hires+survey','colors',[250,0]


;tplot,['comb_hires+survey','comb','fu'+sc+'_fb_col_hires_flux_0_detrend','flux_context_FU'+sc,'flux_context_FU'+sc+'_detrend']
;stop






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

;remove data during dropouts +/- the detrend time
for i=0,n_elements(goo)-1 do begin $
  boo = where((d.x ge d.x[goo[i]-npts_pad]) and (d.x le d.x[goo[i]+npts_pad])) & $
  if boo[0] ne -1 then data[boo] = !values.f_nan
endfor

store_data,'fb_corr',d.x,data


;Compare to Shumko's microburst list
ub = load_firebird_microburst_list('4',filename='FU4_microbursts_bw=0.75sec.csv')
store_data,'ub_shumko_flux',data={x:time_double(ub.time),y:ub.flux_ch1}
options,'ub_shumko_flux','psym',1
options,'ub_shumko_flux','symsize',3
options,'ub_shumko_flux','color',250

ub2 = load_firebird_microburst_list('4',filename='FU4_microbursts_bw=2sec.csv')
store_data,'ub2_shumko_flux',data={x:time_double(ub2.time),y:ub2.flux_ch1}
options,'ub2_shumko_flux','psym',5
options,'ub2_shumko_flux','symsize',3
options,'ub2_shumko_flux','color',250



;Do a better calibration to Shumko microbursts. Mike used a simple conversion from counts to flux for these. I have a much better conversion
;First uncalibrate
div_data,'fu'+sc+'_fb_col_hires_flux_0_backup','fu'+sc+'_fb_col_hires_counts_0'
get_data,'fu4_fb_col_hires_flux_0_backup/fu4_fb_col_hires_counts_0',data=dd
cal2 = mean(dd.y,/nan)

ubcounts = ub.flux_ch1 / cal2
ub2counts = ub2.flux_ch1 / cal2

ubflux_corr = ubcounts * flux_conv
store_data,'ub_shumko_flux_corr',data={x:time_double(ub.time),y:ubflux_corr}
options,'ub_shumko_flux_corr','psym',1
options,'ub_shumko_flux_corr','symsize',3
options,'ub_shumko_flux_corr','color',50

ub2flux_corr = ub2counts * flux_conv
store_data,'ub2_shumko_flux_corr',data={x:time_double(ub2.time),y:ub2flux_corr}
options,'ub2_shumko_flux_corr','psym',5
options,'ub2_shumko_flux_corr','symsize',3
options,'ub2_shumko_flux_corr','color',50



store_data,'detcomb',data=['fu'+sc+'_fb_col_hires_flux_0_detrend','fu'+sc+'_fb_col_hires_flux_0_backup_detrend','fb_corr','ub_shumko_flux','ub_shumko_flux_corr','ub2_shumko_flux','ub2_shumko_flux_corr']
options,'fu'+sc+'_fb_col_hires_flux_0_backup_detrend','color',210
options,'detcomb','colors',[0,150,250,250]
options,'detcomb','ytitle','full detrended data vs!Cgap-reduced version'
ylim,'detcomb',0,0,0



tplot,['comb','datagap_test','detcomb','fu4_fb_col_hires_counts']


;tplot,['detcomb','detcomb_counts']




stop

end
