;Fit a best fit as well as max/min slopes (via lineslope_minmax.pro) to a 
;dispersive microburst with associated errors. 
;****USE IN CONCERT WITH microburst_simulator.pro



;************IMPROVEMENTS***************
;1) USE FITEXY.PRO INSTEAD OF LNFIT -- THIS ALLOWS ERRORS BOTH IN TIME AND ENERGY 
;2) BETTER DETERMINE THE UNCERTAINTY IN TIME BASED ON SIGNAL/NOISE RATIO. 
;************IMPROVEMENTS***************



;option = '1' ;load FIREBIRD Data (SAVES DATA AT END THAT CAN BE LOADED WITH microburst_simulator.pro)
             ;This is how I'll choose the test parameters in microburst_simulator.pro that mimic a real life event.
option = '2'  ;load simulated FIREBIRD microburst


;OPTION 2 ONLY: Choose the spectrum tplot variable 
spec = 'ub_spec_after_detection'
;spec = 'ub_spec_wnoise'   ;full [1000,1000] uB
;spec = 'ub_spec_nonoise'    ;full [1000,1000] uB

;MAKE SURE THE NUMBER OF ENERGY CHANNELS IS SET CORRECTLY 
;nchannels = floor(1000)  ;number of detector ENERGY channels in "spec"
;nchannels = floor(10)  ;number of detector ENERGY channels in "spec"
nchannels = floor(40)  ;number of detector ENERGY channels in "spec"


if nchannels eq 5 or nchannels eq 10 or nchannels eq 20 or nchannels eq 40 then filename = 'fb_ub_'+strtrim(nchannels,2)+'channel_nonoise_recreation_of_20160830_2047-2048.tplot'
if nchannels eq 1000 then filename = 'fb_ub_10channel_nonoise_recreation_of_20160830_2047-2048.tplot'
;fb_ub_5channel_nonoise_recreation_of_20160830_2047-2048.tplot



fluxpeak = fltarr(nchannels)
fluxadj_high = fltarr(nchannels) 
fluxadj_low = fltarr(nchannels) 
tpeak = dblarr(nchannels)

rbsp_efw_init

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

    cadence = 0.05  ;sec

    rbsp_detrend,'fu3_fb_col_hires_flux',0.1

    t0 = time_double('2016-08-30/20:47:27.400')
    t1 = time_double('2016-08-30/20:47:28.600')
    tlimit,t0,t1
    tplot,['fu3_fb_col_hires_flux','fu3_fb_col_hires_flux_smoothed']


    ytmp = tsample('fu3_fb_col_hires_flux_smoothed',[t0,t1],times=tms)

    ;Turn this data into a spectrogram
    ;First adjust the times so that they start at 2014-01-01/00:00, which matches the simulated microbursts
    ;from "option 2"
    tshift = tms[0] - time_double('2014-01-01')

    store_data,'ub_spec_after_detection_realuB',data={x:tms-tshift,y:ytmp,v:ecenter}
    options,'ub_spec_after_detection_realuB','spec',1 
    ylim,'ub_spec_after_detection_realuB',200,1000

    timespan,'2014-01-01',1,/sec 
    tplot,'ub_spec_after_detection_realuB'

stop
endif


;--------------------------------------------------
;Option 2: load simulated FIREBIRD data
;--------------------------------------------------


if option eq '2' then begin 

    path = '/Users/abrenema/Desktop/code/Aaron/github/mission_routines/IMPAX/'
    ;SELECT DATA TO LOAD
    tplot_restore,filename=path + filename



    get_data,spec,data=dd,dlim=dlim
    cadence = dd.x[1] - dd.x[0] ;sec

    if nchannels eq 5 then begin 
        ;Define FIREBIRD detector channels
        ;Values from Crew16 for the collimated detector on FU4
        fblow = [220.,283.,384.,520.,721.]
        fbhig = [283.,384.,520.,721.,985.]
    endif   

    if nchannels eq 10 then begin 
        fblow = [220.,251.,283.,333.,384.,452.,520.,620.,721.,853.]
        fbhig = [251.,283.,333.,384.,452.,520.,620.,721.,853.,985.]
    endif   

    if nchannels eq 20 then begin 
        fblow = [220.,262.,304.,346.,388.,430.,472.,514.,556.,598.,641.,683.,725.,767.,809.,851.,893.,935.,977.,1020.]
        fbhig = [262.,304.,346.,388.,430.,472.,514.,556.,598.,641.,683.,725.,767.,809.,851.,893.,935.,977.,1020.,1040.]
    endif   
    if nchannels eq 40 then begin 
       fblow = [220.00000,240.51282,261.02563,281.53845,302.05127,322.56409,343.07690,363.58972,384.10257,$
        404.61539,425.12820,445.64102,466.15384,486.66666,507.17947,527.69232,548.20514,568.71796,$
        589.23077,609.74359,630.25641,650.76923,671.28204,691.79486,712.30768,732.82050,753.33331,$
        773.84613,794.35895,814.87177,835.38464,855.89746,876.41028,896.92310,917.43591,937.94873,$
        958.46155,978.97437,999.48718,1020.0000]
        fbhig = [240.51282,261.02563,281.53845,302.05127,322.56409,343.07690,363.58972,384.10257,404.61539,$
        425.12820,445.64102,466.15384,486.66666,507.17947,527.69232,548.20514,568.71796,589.23077,$
        609.74359,630.25641,650.76923,671.28204,691.79486,712.30768,732.82050,753.33331,773.84613,$
        794.35895,814.87177,835.38464,855.89746,876.41028,896.92310,917.43591,937.94873,958.46155,$
        978.97437,999.48718,1020.0000,1040.0000]
    endif

    if nchannels eq 1000 then begin 
        fblow = dd.v
        fbhig = shift(dd.v,-1) & fbhig[999] = 1001.
    endif   






    timespan,'2014-01-01/00:00',1,/sec
    ylim,spec,200,1000
    tplot,spec

    tms = dd.x
    ytmp = dd.y

;    ;smooth the data
;    ytmp2 = smooth(ytmp,2)
;    store_data,'tst',data={x:dd.x,y:ytmp2,v:dd.v}
;    options,'tst','spec',1
;    ylim,'tst',200,1000 & zlim,'tst',0,150
;    ytmp = ytmp2


endif



tpeak = dblarr(nchannels)
;extract time of max value in each channel 
for i=0,nchannels-1 do begin
    fluxpeak[i] = max(ytmp[*,i],wh)
    tpeak[i] = tms[wh]
    fluxadj_high[i] = ytmp[wh+1,i]
    fluxadj_low[i] = ytmp[wh-1,i]
;plot,ytmp[*,i]
;stop
endfor

;(CHECK) Plot the peak flux value as well as the adjacent values
for i=0,nchannels-1 do print,fluxadj_low[i],' ',fluxpeak[i],' ',fluxadj_high[i]









;------------------------------------------------------------
;Determine error bars in TIME as well as time values
;------------------------------------------------------------

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



;;reference times to zero for plotting
times = tpeak - tpeak[0]

;Time uncertainty
terror = tmax-tmin



;------------------------------------------------------------
;Determine error bars in ENERGY as well as time values
;------------------------------------------------------------

;Energy error bars come from the finite energy bin size


ecenter = (fblow + fbhig)/2.
eerror = (fbhig - fblow)/2.

;****TEMPORARY 
;eerror[4] = 10000.


;Make error bars huge for lines that have very low signal to noise 
;noiselevel = 20.  ;***TEMPORARY WAY OF IGNORING NOISY DATA. 
noiselevel = 0.  ;***TEMPORARY WAY OF IGNORING NOISY DATA. 
for i=0,nchannels-1 do if fluxpeak[i] le 2*noiselevel then terror[i] = 100.
for i=0,nchannels-1 do if fluxpeak[i] le 2*noiselevel then eerror[i] = 10000.



;p = errorplot(times,ecenter,terror,eerror,xrange=[-2*cadence,20*cadence],$
;xtitle='time (sec, relative)',ytitle='Energy (keV)',title='Microburst')



fit = lineslope_minmax(times,ecenter,eerror,xerr=terror)



;Best fit line
fitline = fit.fitline
;Max slope line 
fitlinemax = fit.fitlinemax
;Min slope line 
fitlinemin = fit.fitlinemin


;High time resolution versions of the fitlines
maxsec = 0.5 
timesHR = maxsec*indgen(1000)/999. - 0.1
;Best fit line
fitlineHR = fit.coeff[1]*timesHR + fit.coeff[0]
;Max slope line 
fitlinemaxHR = (fit.coeff[1]+fit.sigma[1])*timesHR + (fit.coeff[0]-fit.sigma[0])
;Min slope line 
fitlineminHR = (fit.coeff[1]-fit.sigma[1])*timesHR + (fit.coeff[0]+fit.sigma[0])



;Flip axes so that time is on x-axis
nelem = n_elements(fitline)
;xr=[-0.02,1.4*max(times)]
xr = [-0.02,0.1]
yr=[0,1000]

pt = errorplot(times,ecenter,terror,eerror,xtitle='times (relative, sec)',ytitle='Energy (keV)',xrange=xr,yrange=yr)
;pt2 = errorplot(timesHR,fitlineHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
;pt3 = errorplot(timesHR,fitlineminHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
;pt4 = errorplot(timesHR,fitlinemaxHR,replicate(0.,1000),replicate(0.,1000),/overplot,color=[0,200,0])
pt5 = errorplot(times,fitline,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[200,0,0])
pt6 = errorplot(times,fitlinemin,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[0,200,0])
pt7 = errorplot(times,fitlinemax,replicate(0.,nelem),replicate(0.,nelem),/overplot,color=[0,200,0])


;Overplot this onto the spectrogram tplot variable. 

store_data,'epeak',tpeak,ecenter
store_data,'fitline',tpeak,fitline
store_data,'fitlinemin',tpeak,fitlinemin
store_data,'fitlinemax',tpeak,fitlinemax
options,'epeak','psym',4

store_data,'speccomb',data=['ub_spec_after_detection','epeak']
store_data,'speccomb_real',data=['ub_spec_after_detection_realuB','epeak']
options,'speccomb','ytitle','Simulated uB!CkeV'
options,'speccomb_real','ytitle','Real uB!CkeV'
ylim,'speccomb',220,850,0
ylim,'speccomb_real',220,850,0
zlim,'ub_spec_after_detection_realuB',0.01,100,1
zlim,'ub_spec_after_detection',0.01,100,1

tplot,['speccomb_real','speccomb']

store_data,'fitlinecomb',data=['epeak','fitline','fitlinemin','fitlinemax']
ylim,'fitlinecomb',220,1000
tplot,'fitlinecomb'


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



;**********************
print,'Slope analysis (arrival time difference b/t lowest and highest energy (emin/emax))'
;print,'delta-time data',' ',dt_times + ' msec'
print,'delta-time bestfit',' ',1000.*dt_timesFIT + ' msec'
print,'delta-time min',' ',1000.*dt_timesFITmin + ' msec'
print,'delta-time max',' ',1000.*dt_timesFITmax + ' msec'


;****TEST RESULTS WITH NO NOISE
;***ARRIVAL TIME DIFFERENCE OF UB B/T 220 AND 721 KEV
;5 channel --> [45.55,53.55,65.57]  (220-721 keV)
;10 channel--> [55.05,58.56,62.56]  (220-721 keV)
;20 channel -->[57.56,59.06,61.06]  (220-721 keV)
;40 channel -->[58.56,59.06,60.06]  (220-721 keV)
;1000 channel->[60.06] (CORRECT ANSWER FOR SIMULATED DATA)  (220-721 keV)
;REAL DATA --> [58.06,68.07,83.08]  (220-721 keV)
stop

;Save the line fits for the realistic data 
;in order to best modify the simulated data to match the real data. 
save,/variables,filename='~/Desktop/realistic_uB_linefits_20160830_2047-2048'
tplot_save,'ub_spec_after_detection_realuB',filename='~/Desktop/realistic_uB_spec_20160830_2047-2048'

end  


