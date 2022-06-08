;Fit a best fit as well as max/min slopes (via lineslope_minmax.pro) to a 
;dispersive microburst with associated errors in both time and energy. 


;There are two uses for this program:
;option = '1': Fit a line and error slopes to an actual microburst on FIREBIRD. The data is saved
;     at the end of the program so that it can be inputted into microburst_simulator.pro to mimic a real life microburst.
;     This can be useful if you want to then rerun this program (option = '2') to see how a hypothetical detector with higher energy
;     resolution would better detect dispersion
;option = '2': Fit a line and error slopes to a simulated microburst from microburst_simulator.pro


;Adjustable parameters;
;nchannels = number of energy channels in detector. If option='1' is chosen this is set to the FIREBIRD channels. 
;terror = error in time for each datapoint. Calculated based on uncertainty of 1/SNR.    
;eerror = error in energy for each datapoint. Data point is at the center of each energy bin and the error bars extend 
;         to the top and bottom of each bin
;noiselevel = the level of "noise" in the signal. When the microburst amplitude in a given channel is at or less than 
;           this value then its associated terror/eerror values are set to be very high. Note that I've compared this approach to 
;           the approach of ignoring the higher energy channels when their counts are very low and the results are 
;           very comparable. 
;           
            ;****NOTE THAT THERE ARE MULTIPLE SOURCES OF "NOISE"
            ;1) INSTRUMENT NOISE. VERY SMALL AND GENERALLY ONLY IMPORTANT FOR HIGHER ENERGIES WHERE FLUX IS LOW
            ;2) NOISE FROM ADJACENT MICROBURSTS AND VARIATIONS IN SMALL-SCALE PRECIPITATION. 

;NOTE: it's better to remove datapoints where a clear microburst peak can't be identified rather than 
;keep it but assign it a large time errorbar.

;TESTING: 
; test_plot1 - shows a comparison b/t a real microburst (left) and a simulated one (right) made to go through a detector 
;              identical to the real FIREBIRD detector. Results are very close, indicating that whole process this is working. 


;---------------------------------------------------------------------------------------------------------


;option = '1' ;load FIREBIRD Data (SAVES DATA AT END THAT CAN BE LOADED WITH microburst_simulator.pro)
             ;This is how I'll choose the test parameters in microburst_simulator.pro that mimic a real life event.
option = '2'  ;load simulated FIREBIRD microburst (from microburst_simulator.pro)



;OPTION 2 ONLY: Choose the spectrum tplot variable 
spec = 'ub_spec_after_detection'  ;after artificial uB goes through fake detector
;spec = 'ub_spec_wnoise'   ;full [1000,1000] uB
;spec = 'ub_spec_nonoise'    ;full [1000,1000] uB


;Select filename (don't include .tplot)
;filename = 'fb_ub_63channel_cadence=20msec_Emin=132_Emax=1007_recreation_of_20160830_2047-2048'
;filename = 'fb_ub_63channel_cadence=20msec_Emin=82_Emax=1007_recreation_of_20160830_2047-2048'
;filename = 'fb_ub_32channel_cadence=20msec_Emin=139_Emax=1014_recreation_of_20160830_2047-2048'
;filename = 'fb_ub_16channel_cadence=20msec_Emin=154_Emax=1029_recreation_of_20160830_2047-2048'
;filename = 'fb_ub_6channel_cadence=20msec_Emin=212_Emax=1087_recreation_of_20160830_2047-2048'
filename = 'fb_ub_6channel_cadence=20msec_Emin=136_Emax=1096_recreation_of_20160830_2047-2048'

;filename = 'fb_ub_5channel_cadence=20msec_Emin=251_Emax=853_recreation_of_20160830_2047-2048'
;filename = 'fb_ub_63channel_cadence=20msec_Emin=256_Emax=1006_recreation_of_20160830_2047-2048'


;filename = 'fb_ub_63channel_cadence=20msec_Emin=82_Emax=1007_recreation_of_20160830_2047-2048_nonoise'
;filename = 'fb_ub_63channel_cadence=20msec_Emin=82_Emax=1007_recreation_of_20160830_2047-2048_onlyinstrumentNoise'
;if nchannels eq 5 or nchannels eq 10 or nchannels eq 20 or nchannels eq 40 then filename = 'fb_ub_'+strtrim(nchannels,2)+'channel_recreation_of_20160830_2047-2048'
;if nchannels eq 1000 then filename = 'fb_ub_10channel_recreation_of_20160830_2047-2048'
;fb_ub_5channel_nonoise_recreation_of_20160830_2047-2048.



;;MAKE SURE THE NUMBER OF ENERGY CHANNELS IS SET CORRECTLY 
;sstr = strsplit(filename,'_',/extract)
;sstrtmp = sstr[2]
;tmp = strpos(sstrtmp,'channel')
;nchannels = float(strmid(sstrtmp,0,tmp))







rbsp_efw_init
device,decomposed=0
loadct,39

charsz_plot = 0.8  ;character size for plots
charsz_win = 1.2
!p.charsize = charsz_win
tplot_options,'xmargin',[20.,15.]
tplot_options,'ymargin',[3,6]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',2
tplot_options,'ythick',2
tplot_options,'labflag',-1



;--------------------------------------------------
;Option 1: load FIREBIRD microbursts 
;--------------------------------------------------

if option eq '1' then begin 

    nchannels = floor(5)  ;override previously set value
    fblow = [220.,283.,384.,520.,721.]
    fbhig = [283.,384.,520.,721.,985.]
    ecenter = (fblow + fbhig)/2.

    timespan,'2016-08-30'     
    firebird_load_data,'3'    


    rbsp_detrend,'fu3_fb_col_hires_flux',0.1

    ;Select rough start and stop times of microburst. A good approach is to select times where only a single 
    ;microburst is present. 
    t0 = time_double('2016-08-30/20:47:28.400')
    t1 = time_double('2016-08-30/20:47:29.600')
    device,decomposed=0
    loadct,39
    tplot,['fu3_fb_col_hires_flux','fu3_fb_col_hires_flux_smoothed']
    tlimit,t0,t1
    stop

    ytmp = tsample('fu3_fb_col_hires_flux_smoothed',[t0,t1],times=tms)

    ;Turn this data into a spectrogram
    ;First adjust the times so that they start at 2014-01-01/00:00, which matches the simulated microbursts
    ;from "option 2"
    tshift = tms[0] - time_double('2014-01-01')

    store_data,'ub_spec_after_detection',data={x:tms-tshift,y:ytmp,v:ecenter}
    options,'ub_spec_after_detection','spec',1 
    ylim,'ub_spec_after_detection',200,1000
    zlim,'ub_spec_after_detection',0,5


    timespan,'2014-01-01',1,/sec 
    ylim,'ub_spec_after_detection',250,850
    tplot,'ub_spec_after_detection'
    

    get_data,'ub_spec_after_detection',data=dd
    tms = dd.x
    ytmp = dd.y

    cadence = dd.x[1] - dd.x[0] ;sec

endif


;--------------------------------------------------
;Option 2: load simulated FIREBIRD data
;--------------------------------------------------


if option eq '2' then begin 

    path = '/Users/abrenema/Desktop/code/Aaron/github/mission_routines/IMPAX/'
    ;SELECT DATA TO LOAD
    tplot_restore,filename=path + filename + '.tplot'
    restore,path + filename + '.sav'


    if nchannels eq 1000 then begin 
        fblow = dd.v
        fbhig = shift(dd.v,-1) & fbhig[999] = 1001.
    endif   

    ecenter = (fblow + fbhig)/2.




    timespan,'2014-01-01/00:00',1,/sec

 
    copy_data,spec,'ub_spec_after_detection'
    get_data,'ub_spec_after_detection',data=dd
    tms = dd.x
    ytmp = dd.y
    cadence = dd.x[1] - dd.x[0] ;sec


endif



fluxpeak = fltarr(nchannels)
fluxadj_high = fltarr(nchannels)
fluxadj_low = fltarr(nchannels)
tpeak = dblarr(nchannels)





tpeak = dblarr(nchannels)
;extract time of max value in each channel 
for i=0,nchannels-1 do begin
    fluxpeak[i] = max(ytmp[*,i],wh)

    tpeak[i] = tms[wh]
    if wh le n_elements(tms)-2 then fluxadj_high[i] = ytmp[wh+1,i] else fluxadj_high[i] = ytmp[wh,i] 
    if wh gt 0 then fluxadj_low[i] = ytmp[wh-1,i] else fluxadj_low[i] = ytmp[wh,i] 

endfor

;(CHECK) Plot the peak flux value as well as the adjacent values
for i=0,nchannels-1 do print,fluxadj_low[i],' ',fluxpeak[i],' ',fluxadj_high[i]





;------------------------------------------------------------
;Determine error bars in TIME as well as time values
;------------------------------------------------------------


;(1)---------------------------------------------
;Define time error based on sampling rate - This is used as a backup value if 
;the determination based on SNR doesn't work. 


;Error bars in time will extend to the halfway point b/t adjacent values
;*****TEMPORARY - COME UP WITH BETTER WAY TO DECIDE TIME ERROR BARS
tmin = (tpeak - cadence/2.)
tmax = (tpeak + cadence/2.)


;For each peak, if either adjacent data point has the same value then extend the error bar
;(This happens frequently)
for i=0,4 do begin $
    if fluxpeak[i] eq fluxadj_low[i] then tmin[i] -= cadence/2. & $
    if fluxpeak[i] eq fluxadj_high[i] then tmax[i] += cadence/2.
endfor


;Reshift peak values if either tmin or tmax has been extended due to 
;neighboring data points with the same value 
goohigh = where((tmax - tpeak) gt cadence)
goolow = where((tpeak - tmin) gt cadence)
if goohigh[0] ne -1 then for i=0,n_elements(goohigh)-1 do tpeak[goohigh[i]] += cadence/2.
if goolow[0] ne -1 then for i=0,n_elements(goolow)-1 do tpeak[goolow[i]] -= cadence/2.

terrorTMP = (tmax-tmin)/2.


;(2)---------------------------------------------
;Define time error based on SNR 
;****THIS WORKS VERY POORLY FOR BELL-CURVES B/C THE STD IS VERY HIGH. 
;****MAY BE BETTER TO TAKE THE STD OF THE HIGH-PASSED DATA. 
;*****;*****;*****;*****;*****;*****;*****;*****
;*****;*****;*****;*****;*****;*****;*****;*****


;stderr = fltarr(n_elements(dd.v))
;for i=0,nchannels-1 do stderr[i] = stddev(dd.y[*,i])/sqrt(n_elements(dd.x))

;SNR (ratio of mean to stdev)
snr = fltarr(n_elements(dd.v))
for i=0,nchannels-1 do snr[i] = mean(dd.y[*,i])/stddev(dd.y[*,i])

snr = snr^2.
uncertainty = 1/snr  ;percentage uncertainty




tmin = (tpeak - cadence/2.)
tmax = (tpeak + cadence/2.)
terror = fltarr(n_elements(dd.v))
;testvalue = fltarr(n_elements(dd.v))

;Now determine how far error bars must extend to left/right of peak based on the stdev
for i=0,nchannels-1 do begin
  ;testvalue = max(dd.y[*,i],whmax) - stderr[i]
  testvalue = max(dd.y[*,i],whmax) * (1 - uncertainty[i])
  print,max(dd.y[*,i],whmax)
  print,testvalue
  print,'****'
  goo = where(dd.y[*,i] ge testvalue)
  if goo[0] ne -1 then begin 
    elocL = goo[0]
    elocR = goo[-1]
    terror[i] = dd.x[elocR] - dd.x[elocL]
    if terror[i] eq 0. then terror[i] = (tmax[goo] - tmin[goo])/2.
  endif else terror[i] = 2.*terrorTMP[i]

;  loadct,39
;  plot,dd.x,dd.y[*,i]
;  oplot,[time_double('2014-01-01'),time_double('2014-01-01/00:00:01')],[testvalue,testvalue]
;  oplot,[dd.x[whmax]-terror[i],dd.x[whmax]-terror[i]],[0,1000]
;  oplot,[dd.x[whmax]+terror[i],dd.x[whmax]+terror[i]],[0,1000]
;  stop
endfor

;goo = where(terror eq 0.)
;if goo[0] ne -1 then terror[goo] = (tmax[goo] - tmin[goo])/2.  




ecenter = (fbhig + fblow)/2.
loadct,39


;Compare the error from using only the sampling rate to the error determined from the SNR.
plot,ecenter,terror
oplot,ecenter,terrorTMP,color=250
stop





;;reference times to zero for plotting and fitting
times = tpeak - tpeak[0]






;stop

;------------------------------------------------------------
;Determine error bars in ENERGY as well as time values
;------------------------------------------------------------

;Energy error bars come from the finite energy bin size


ecenter = (fblow + fbhig)/2.
eerror = (fbhig - fblow)/2.

;****TEMPORARY 
;eerror[4] = 10000.


;;Remove data points that have very low signal to noise 
;Test highest bin to see what noise level should be. 


;*************
;Set noise level (higher values mean a datapoint is less likely to exceed noise)
if option eq '1' then noiselevel = 0.02
if option eq '2' then noiselevel = 0.02*noise_med  ;set from microburst_simulator.pro
;noiselevel = 0.18
good = where(fluxpeak ge noiselevel)
if good[0] ne -1 then begin
  times = times[good]
  ecenter = ecenter[good]
  eerror = eerror[good]
  terror = terror[good]
endif
  



fit = lineslope_minmax(times,ecenter,eerror,xerr=terror)




;Best fit line
fitline = fit.fitline
;Max slope line 
fitlinemax = fit.fitlinemax
;Min slope line 
fitlinemin = fit.fitlinemin


;*************FIX THIS TO BE MORE GENERAL***************
;High time resolution versions of the fitlines
maxsec = 0.5 
timesHR = maxsec*indgen(1000)/999. - 0.2
;Best fit line
fitlineHR = fit.coeff[1]*timesHR + fit.coeff[0]
;Max slope line 
fitlinemaxHR = (fit.coeff[1]+fit.sigma[1])*timesHR + (fit.coeff[0]-fit.sigma[0])
;Min slope line 
fitlineminHR = (fit.coeff[1]-fit.sigma[1])*timesHR + (fit.coeff[0]+fit.sigma[0])


stop




;-------------------------------------------------------
;Overplot this onto the spectrogram tplot variable. 
store_data,'epeak',tpeak,ecenter
store_data,'fitline',tpeak,fitline
store_data,'fitlinemin',tpeak,fitlinemin
store_data,'fitlinemax',tpeak,fitlinemax

;;Hires version
store_data,'fitlineHR',tpeak[0]+timesHR,fitlineHR
store_data,'fitlineminHR',tpeak[0]+timesHR,fitlineminHR
store_data,'fitlinemaxHR',tpeak[0]+timesHR,fitlinemaxHR


options,['fitlineHR','fitline'],'thick',2 & options,['fitlineHR','fitline'],'color',0
options,['fitlineminHR','fitlinemaxHR'],'thick',2 & options,['fitlineminHR','fitlinemaxHR'],'color',250
options,['fitlinemin','fitlinemax'],'thick',2 & options,['fitlinemin','fitlinemax'],'color',250



options,'epeak','psym',4
options,'epeak','thick',3
options,'epeak','color',99

times_left = tpeak - terror
times_right = tpeak + terror
store_data,'epeak_left',times_left,ecenter
store_data,'epeak_right',times_right,ecenter
options,'epeak_left','psym',2
options,'epeak_right','psym',2
options,'epeak_left','color',220
options,'epeak_right','color',220


copy_data,'ub_spec_after_detection','ub_spec_after_detection_line'
options,'ub_spec_after_detection_line','spec',0
get_data,'ub_spec_after_detection_line',data=dtmp
if is_struct(dtmp) then begin
  store_data,'ub_spec_after_detection_line_0',dtmp.x,dtmp.y[*,0]
  store_data,'ub_spec_after_detection_line_1',dtmp.x,dtmp.y[*,1]
  store_data,'ub_spec_after_detection_line_2',dtmp.x,dtmp.y[*,2]
  store_data,'ub_spec_after_detection_line_3',dtmp.x,dtmp.y[*,3]
  store_data,'ub_spec_after_detection_line_4',dtmp.x,dtmp.y[*,4]
endif

;split_vec,'ub_spec_after_detection_line' 
ylim,'ub_spec_after_detection_line_?',0,0,0
if is_struct(dtmp) then ylim,'ub_spec_after_detection_line',0.01,max(dtmp.y,/nan),1



store_data,'speccomb_real',data=['ub_spec_after_detection','epeak','fitlineHR','fitlineminHR','fitlinemaxHR','epeak_left','epeak_right']
options,'speccomb_real','ytitle','Real uB!CkeV'
ylim,'speccomb',200,1000,0
ylim,'speccomb_real',200,1000,0
zlim,'ub_spec_after_detection',0.01,100,1
zlim,'ub_spec_after_detection',0.01,100,1


;store_data,'fitlinecomb',data=['epeak','fitline','fitlinemin','fitlinemax','epeak_left','epeak_right']
;ylim,'fitlinecomb',200,1000

options,'speccomb_real','panel_size',6
options,'ub_spec_after_detection_line_?','panel_size',0.7



;Plot for comparison b/t detectors of different energy channels
if option eq '1' then options,'speccomb_real','ytitle','real uB!CkeV'
if option eq '2' then options,'speccomb_real','ytitle','simulated uB!CkeV'
tplot_options,'title','Best fit for nchannels='+strtrim(nchannels,2)
ylim,'speccomb_real',0,1000,0
loadct,39
tplot,'speccomb_real'


if nchannels eq 5 then begin
  loadct,39
  tplot,['ub_spec_after_detection_line_4',$
  'ub_spec_after_detection_line_3',$
  'ub_spec_after_detection_line_2',$
  'ub_spec_after_detection_line_1',$
  'ub_spec_after_detection_line_0',$
  'speccomb_real']
endif



;--------------------------------------------------------------------
;Find the delta-times for each of the fit lines for the ray tracing. 
;nelem = n_elements(times)
;dt_times = times[nelem-1] - times[0]
;dt_timesFIT = fitline[nelem-1] - fitline[0]
;dt_timesFITmin = fitlinemax[nelem-1] - fitlinemax[0]
;dt_timesFITmax = fitlinemin[nelem-1] - fitlinemin[0]

emin = 220. 
emax = 721.

goo1 = where(fitlineHR ge emin)
goo2 = where(fitlineHR ge emax)
dt_timesFIT = timesHR[goo2[0]] - timesHR[goo1[0]]
goo1 = where(fitlineminHR ge emin)
goo2 = where(fitlineminHR ge emax)
dt_timesFITmin = timesHR[goo2[0]] - timesHR[goo1[0]]
goo1 = where(fitlinemaxHR ge emin)
goo2 = where(fitlinemaxHR ge emax)
dt_timesFITmax = timesHR[goo2[0]] - timesHR[goo1[0]]


;--------------------------------------------------
;Print results of the linefits for use in microburst_simulator.pro
;--------------------------------------------------

print,'*******************************'
print,filename
print,'Slope analysis (arrival time difference b/t lowest and highest energy (emin/emax))'
print,'Probability of good fit ' + string(fit.prob_of_good_fit)
print,'Emin = ' + strtrim(Emin) + ' ; Emax = ' + strtrim(Emax)
print,'delta-time bestfit',' ',string(1000.*dt_timesFIT) + ' msec'
print,'delta-time min',' ',string(1000.*dt_timesFITmin) + ' msec'
print,'delta-time max',' ',string(1000.*dt_timesFITmax) + ' msec'
;-----------
print,'Best fit line and error slope equations (y=mx+b)
print,'BEST FIT: y = '+strtrim(fit.coeff[1],2)+'x + ' + strtrim(fit.coeff[0],2)
print,'MAX SLOPE: y = '+strtrim((fit.coeff[1]+fit.sigma[1]),2)+'x + ' + strtrim((fit.coeff[0]-fit.sigma[0]),2)
print,'MIN SLOPE: y = '+strtrim((fit.coeff[1]-fit.sigma[1]),2)+'x + ' + strtrim((fit.coeff[0]+fit.sigma[0]),2)





;****TEST RESULTS WITH NO NOISE
;***ARRIVAL TIME DIFFERENCE OF UB B/T 220 AND 721 KEV
;5 channel --> [45.55,53.55,65.57]  (220-721 keV)
;10 channel--> [55.05,58.56,62.56]  (220-721 keV)
;20 channel -->[57.56,59.06,61.06]  (220-721 keV)
;40 channel -->[58.56,59.06,60.06]  (220-721 keV)
;1000 channel->[60.06] (CORRECT ANSWER FOR SIMULATED DATA)  (220-721 keV)
;REAL DATA --> [58.06,68.07,83.08]  (220-721 keV)
stop
stop
stop


;For option 1 only, save the line fits for the realistic data
;in order to best modify the simulated data to match the real data.
if option eq '1' then begin 
  save,/variables,filename='~/Desktop/realistic_uB_linefits_20160830_2047-2048'
  ;tplot_save,'ub_spec_after_detection',filename='~/Desktop/realistic_uB_spec_20160830_2047-2048'
  tplot_save,'*',filename='~/Desktop/realistic_uB_spec_20160830_2047-2048'
endif


end  




;***OBSOLETE
;-------------------------------------------------------
;Make a plot that shows the fit lines and error bars
;-------------------------------------------------------
;


;;Flip axes so that time is on x-axis
;nelem = n_elements(fitline)
;;xr=[-0.02,1.4*max(times)]
;
;pt = errorplot(times,ecenter,terror,eerror,xtitle='times (relative, sec)',ytitle='Energy (keV)')
;pt.xrange = [-0.25,0.25]
;pt.yrange = yr
;;pt2 = errorplot(timesHR,fitlineHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
;;pt3 = errorplot(timesHR,fitlineminHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
;;pt4 = errorplot(timesHR,fitlinemaxHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
;pt5 = errorplot(times,fitline,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[200,0,0])
;pt6 = errorplot(times,fitlinemin,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[0,200,0])
;pt7 = errorplot(times,fitlinemax,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[0,200,0])

