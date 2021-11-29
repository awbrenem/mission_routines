;Test 19 dB correction to EMFISIS data. See Malaspina email on Feb 2, 2018.

;RBSP-B: on 2017-02-16, a gain level of 19dB applied by EMFISIS to their
;SCM waveform data was turned off. After this date, EFW needs to boost the SCM
;signal by 19dB.

;2017-02-17

;date = '2015-02-17'
;t0 = time_double(date + '21:32')
;t1 = time_double(date + '21:33')
;bt = '2'
;sc = 'b'

;ATTENUATOR END OF DAY FLIP (20:51:45)
date = '2015-09-22'
bt = '2'
sc = 'a'
;t0 = time_double(date + '/05:45:10')
;t1 = time_double(date + '/05:46:00')
;t0z = time_double(date + '/05:45:10')
;t1z = time_double(date + '/05:46:00')


;ATTENUATOR ON
;date = '2015-02-10'
;bt = '1'
;sc = 'b'
;t0 = time_double(date + '/05:45:10')
;t1 = time_double(date + '/05:46:00')
;t0z = time_double(date + '/05:45:10')
;t1z = time_double(date + '/05:46:00')

;;ATTENUATOR OFF
;date = '2015-03-15'
;bt = '1'
;sc = 'b'
;t0 = time_double(date + '/03:09:01.400')
;t1 = time_double(date + '/03:09:02.600')
;t0z = time_double(date + '/03:09:01.400')
;t1z = time_double(date + '/03:09:02.600')

date = '2015-02-16'  ;RBSPa attenuator ON halfway through day
bt = '2'
sc = 'a'


;date = '2015-11-08'  ;RBSPa attenuator ON
;t0 = time_double(date + '/00:00')
;t1 = time_double(date + '/00:10')
;t0z = time_double(date + '/00:04:21')
;t1z = time_double(date + '/00:04:24')
;bt = '2'
;sc = 'a'

;Load burst data
timespan,date
rbsp_load_efw_waveform,probe=sc,type='calibrated',$
  datatype=['mscb'+bt];,trange=[t0,t1]
yv = tsample('rbsp'+sc+'_efw_mscb'+bt,[t0z,t1z],times=tt)
;burstmag = sqrt(yv[*,0]^2 + yv[*,1]^2 + yv[*,2]^2)
burstmag = yv[*,2]
power_burst = fft_power_calc(tt,burstmag,/read_win,samp_ra=16384.)   ;units^2/Hz
help,power_burst
print,max(power_burst.freq)
;powerburst = power_burst.power_a/10^(19./10.)  ;reduce power that has been incorrectly artificially boosted


;Load EMFISIS DC field
rbsp_load_emfisis,probe=sc,type='calibrated',coord='uvw',cadence='hires',level='l2'
split_vec,'rbsp'+sc+'_emfisis_l2_uvw_Mag',suffix='_'+['u','v','w']
yv = tsample('rbsp'+sc+'_emfisis_l2_uvw_Mag_w',[t0z,t1z],times=tt)
store_data,'BoDC_mag',tt,yv
power_bodc = fft_power_calc(tt,yv,/read_win,samp_ra=64.)  ;units^2/Hz
help,power_bodc
print,max(power_bodc.freq)

;Load EMFISIS spectral data
;path = '~/Downloads/rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_20151108_v1.6.3.cdf'
;path = '~/Downloads/rbsp-b_WFR-spectral-matrix-diagonal_emfisis-L2_20150210_v1.4.3.cdf'
path = '~/Downloads/rbsp-b_WFR-spectral-matrix-diagonal_emfisis-L2_20150315_v1.4.6.cdf'

cdf2tplot,path
get_data,'BwBw',data=dspec
store_data,'BwBw',data={x:dspec.x,y:dspec.y,v:reform(dspec.v)}
powspec = tsample('BwBw',[t0-0.5*60.,t1+0.5*60.],times=tt)


;Plot spectra of all quantities for comparison.
plot,reform(dspec.v),powspec[0,*],xrange=[0.01,10000],/xlog,xstyle=1,/ylog,yrange=[1d-14,1d-2],/nodata
oplot,power_bodc.freq,power_bodc.power_a
oplot,power_burst.freq,powerburst,color=250
for i=0,n_elements(powspec[*,0])-1 do oplot,reform(dspec.v),powspec[i,*],color=50
