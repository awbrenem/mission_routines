;Create a microburst and simulate its detection at the location of a satellite (e.g. FIREBIRD) as a function of energy and time.
;***To get microburst best fit and error slopes first run microburst_fit_slope.pro 

;*********************************************
;TODO: incorporate 10-50msec integration time that FIREBIRD uses for uB simulation.
;*********************************************



;f(E,t) = A(t)*exp(-1*(E - Eo)/deltaE(t))

;A(t) = For each time step this defines the peak amplitude value and the energy it occurs at.
;dE = For each time step this defines the width of Gaussian curve
;Eo = time-offset of curve. Currently this is time-independent


;We want to define the functional forms of A, Eo, and dE such that
;the total flux (integrated over all energies) mimicks what probably happens in real life. Namely,
;it should start low, ramp up, then ramp down.

;Steps:
;1) define f(E,t) as time-varying exponential f(E) = fo*exp(-E/Eo)
;2) For each individual channel multiply f(E) by energy response (Gaussian).
;   Fch_integrated = integral(f(E)*fch(E)) from 0 to infinity
;   This is the instantaneous integrated energy response. If FB detectors respond very quickly (which they do)
;   than we may not need to worry about a ramp-up time.
;3) Do this for all the channels and combine the signals.

;-------
;Gets more complicated if FB has a slow ramp-up response time.
;for ex: f(t) = int(F(E,t)) from 250-400



;******************************************
;Test to see how well we'll be able to resolve the microbursts 
;given an instrument noise level plus a background level of 
;"microburst noise" that I think occurs due to edge-detection of 
;other microbursts. 
;*************************************
;Highest bin with visible microburst --> 2016-08-30/20:47:28 dispersive smallish amplitude uB
;f0tmp = 120. ;Flux at 0 keV
;E0tmp = 100.
;epow = 1.2
  ;5bin : 721-985
  ;10bin: 721-853
  ;20bin: 767-809
  ;40bin: 753-773
  ;65bin: 720-732
;**RESULT: EVEN WITH 65 BINS WE CAN RESOLVE THIS SMALLISH AMPLITUDE MICROBURST
;*************************************

;Highest bin with visible microburst () --> typical FB microburst using Arlo's statistics
;More typical uB with epow=1
;f0tmp = 1000. ;Flux at 0 keV
;E0tmp = 70.
;epow = 1.
  ;5bin :
  ;10bin:
  ;20bin: 
  ;40bin: 
  ;50bin: 
  ;60bin: 
;*************************************



rbsp_efw_init
!p.charsize = 1.5

;***********************************
;***********************************
;Harlan says 12 keV channels doable
;FWHM is 12 keV for noise. (typically set threshold 3x noise)
;Can go even lower in energy. 
;--mention that an even lower energy threshold would make dispersion more obvious (Saito)

;Max noise at 220 keV (65 bins) is about 5 in flux units. 
;Instrument noise is constant across bins at 0.05 in flux units and is only apparent in highest channels

;Estimate of ~31 counts/flux_unit 

;***********************************


;Load up spec of real microburst (from microburst_fit_slope.pro) so that I can set the 
;simulated uB parameters
datetime_real = '20160830_2047-2048'
tplot_restore,filename='/Users/abrenema/Desktop/code/Aaron/github/mission_routines/IMPAX/realistic_uB_spec_'+datetime_real+'.tplot'
;Note that times for the realistic spec have been shifted so that they start at 2014-01-01, just like
;the simulated microburst here. 
loadct,39
zlim,'ub_spec_after_detection_realuB',0,0,0
tplot,['ub_spec_after_detection_realuB']

;Determine max flux of 220 keV bin 
get_data,'ub_spec_after_detection_realuB',data=dd

;Flux max for 220 keV bin
print,max(dd.y[*,0])
;       16.475281

stop

;------------------------------------------------------------------------
;Define microburst parameters
;Can define flux at 220 keV from either
;1) direct observations of flux at 220 keV from actual microburst
;2) flux at a lower energy extrapolated to 220 keV based on Arlo's fits.
;------------------------------------------------------------------------

;---
;Flux max for the ~220 keV FIREBIRD bin (from observations)
f0 = max(dd.y[*,0])
print,f0
;f0 = 19.   ;flux
;f0 = 1000.  ;counts   


;Determine f0 (flux in 1/(cm2-s-sr-keV) at 200 keV) by extrapolating Arlo's flux results (for 0 keV) to 220 keV. 

;;Typical uB
;f0tmp = 1000. ;Flux at 0 keV
;E0tmp = 70.

;;uB at 2016-08-30/20:47:28
;f0tmp = 140. ;Flux at E0tmp (keV)
;E0tmp = 100.
;
;f0 = f0tmp*exp(-1*220./E0tmp)  ;typical uB flux at 220 keV
;---


;---
;Set detector cadence
cadence_newdetector = 0.05 ;sec 
;---


;---
;multiplicative factor for the noise.
noise_scalefac = 0.01
;---


;---
;Add in contaminating microburst "noise" at sample rate cadence
;These are the sample rate bumps you see in real data that may correspond to 
;edge-detected microbursts?
;This "noise" is channel dependent. I'll set it at a fraction of the peak value 
;detected in each channel (noise_from_uB_fraction). NOTE: this value will need to be 
;abnormally high b/c the detectors integrate over many adjacent (1000x) energy bins and 
;the noise tends to disappear to unrealistically low values. For example, if you want the 
;noise in the realistic uB (say energy channel 500) to be 10% of the max value in that channel, 
;you would set noise_from_uB_fraction >> 0.1. 
;Just adjust so that things look real. 
;noise_from_uB_fraction = 3.  ;--> determined by comparison to 2016-08-30/20:47:28 dispersive event
noise_from_uB_fraction = 1.  ;--> determined by comparison to 2016-08-30/20:47:28 dispersive event
;---







;---
;energy profile ("infinite" resolution uB) --> exp(-(E-Eo)^epow/dE^epow)
nenergies = 1000. ;number of microburst energy steps
emin_keV = 0. 
emax_keV = 1000.
energies_uB = (emax_keV - emin_keV)*indgen(nenergies)/(nenergies-1) + emin_keV   ;creates a 1000 element array
Eo = 0. ;Center energy for Gaussian
dE = 75. ;Width of energy Gaussian curve (keV)
epow = 1.  ;typically=1     flux[*,i] *= exp(-1*(energies_uB - Eo)^epow/dE^epow)
;---


;---
;time profile  (exp(-(t-to)^2/dt^2))
ntsteps = 1000. ;number of microburst time steps for "infinite" resolution uB
tend_sec = 1.15  ;total time duration considered in seconds (microburst contained within this)
;dt = 0.2  ;Width (sec) of the exponential falloff
dt = 0.1  ;Width (sec) of the exponential falloff
to = 0.56  ;time offset for the zeroth bin (sec). Choose to align with actual microburst from microburst_fit_slope.pro
t_dispersion = 0.134 ;sec  - delta time b/t lowest and highest energy channels from best fit line (emin_keV; emax_keV)
;---




;-----------------------------------------------------------------
;;Define FIREBIRD detector channels
;-----------------------------------------------------------------



;;Values from Crew16 for the collimated detector on FU4
fblow = [220.,283.,384.,520.,721.]
fbhig = [283.,384.,520.,721.,985.]
fb_comparison = 1.


;;Stupidly fat channels
;fblow = [220.,520.]
;fbhig = [520.,721.]
;fb_comparison = 0.


;;Hypothetical 10 channels
;fblow = [220.,251.,283.,333.,384.,452.,520.,620.,721.,853.]
;fbhig = [251.,283.,333.,384.,452.,520.,620.,721.,853.,985.]
;;plot,indgen(n_elements(fblow)),fblow,psym=-2
;;oplot,indgen(n_elements(fblow2))/2.,fblow2,color=250,psym=-4
;fblow = fblow2 
;fbhig = fbhig2
;fb_comparison = 0.


;;Extremely hires channels
;nch = 65.
;fblow = (1000 - 200.)*indgen(nch)/(nch-1) + 220.
;fbhig = shift(fblow,-1)
;fbhig[n_elements(fbhig)-1] = 1040.
;fb_comparison = 0.

;nch = 5.
;fblow = (1000 - 200.)*indgen(nch)/(nch-1) + 220.
;fbhig = shift(fblow,-1)
;fbhig[n_elements(fbhig)-1] = 1040.
;fb_comparison = 0.


nchannels = n_elements(fblow)
binwidth = fbhig - fblow
Ecenter = (fblow + fbhig)/2.

;String channel names
chnameslow = strtrim((floor(fblow)),2)
chnameshig = strtrim((floor(fbhig)),2)


;----------------------------------------------------------------------------
;Define time variation of flux --> exp(-t^2/dt^2) for first energy bin. 
;----------------------------------------------------------------------------


tstep_sec = tend_sec/ntsteps   ;number of seconds per time tick
times_sec = indgen(ntsteps)*tstep_sec ;times (sec) b/t 0 and tend_sec

flux_tprofile_zeroenergy = exp(-1.*times_sec^2/dt^2)
flux_tprofile_zeroenergy = [reverse(flux_tprofile_zeroenergy),flux_tprofile_zeroenergy] ;full Gaussian

!p.multi = [0,0,1]
ttmp = [-1*reverse(times_sec),times_sec]
plot,ttmp,flux_tprofile_zeroenergy,xtitle='time(sec)',ytitle='Flux profile for zeroth energy'


;------------------------------------------------------------------------------
;Now add in dispersion for the higher energy bins
;------------------------------------------------------------------------------

;First calculate the amount of dispersion from one energy channel to the next. 
;This [e,t] array consists of only normalized Gaussians^2 for  

dt_sec_singlechannel_jump = t_dispersion/nenergies   ;dispersive time shift (sec) from one energy channel to the next

;sec_per_tbin = dt_sec_singlechannel_jump*ntsteps
sec_per_tbin = tend_sec/ntsteps  ;how many seconds in each time grid chunk?


dt_bin_singlechannel_jump = dt_sec_singlechannel_jump/sec_per_tbin   ;number of grid time steps to shift by to have every successive energy peak arrive dt (sec) later 

offset = to/sec_per_tbin

flux_tprofile = dblarr(nenergies,ntsteps)

for i=0,nenergies-1 do begin $
  tmp = shift(flux_tprofile_zeroenergy,dt_bin_singlechannel_jump*i + offset) & $
  flux_tprofile[i,*] = tmp[ntsteps:(ntsteps*2)-1]
endfor


;!p.multi = [0,0,1]
;plot,times_sec,flux_tprofile[40.,*],yrange=[0,1],xtitle='time(sec)',ytitle='flux for specific energies',title='Dispersion test'
;for e=0,nenergies-1 do oplot,times_sec,flux_tprofile[e,*],color=e/20

t0 = '2014-01-01'  ;dummy time
tms = time_double(t0) + double(times_sec)

store_data,'ub_timevariation',data={x:tms,y:transpose(flux_tprofile),v:energies_uB}
options,'ub_timevariation','spec',1
options,'ub_timevariation','ytitle','Energy (keV)'
options,'ub_timevariation','xtitle','time (unitless)'
options,'ub_timevariation','title','Simulated Microburst time variation (normalized)'

options,'fitlineHR','colors',0
store_data,'tmpcomb',data=['ub_timevariation','fitlineHR']
store_data,'tmpcomb2',data=['ub_spec_after_detection_realuB','fitlineHR']

ylim,['tmpcomb','tmpcomb2'],0,1000
zlim,'ub_timevariation',0,1,0
zlim,'ub_spec_after_detection_realuB',0.01,20,1
loadct,39
;Plot the normalized dispersive time signature. NOTE: don't expect the magnitude to match real microburst at this point.
;***check to see that the dispersion is similar to real microburst 
tplot,['tmpcomb','tmpcomb2']



stop






;Construct flux array f(E,t). It'll have size [nenergies, ntsteps]
flux = flux_tprofile 
for i=0,ntsteps-1 do flux[*,i] *= exp(-1*((energies_uB - Eo)^epow)/(dE^epow))

;stop

;change from normalized flux to real values (from FIREBIRD)
;Since f0 is the max flux at 220 keV, we'll need to determine the max flux at 0 keV
f0_0keV = f0/exp(-1*((220.-Eo)^epow)/(dE^epow))


flux *= f0_0keV

store_data,'ub_spec_nonoise',data={x:tms,y:transpose(flux),v:energies_uB}
options,'ub_spec_nonoise','spec',1
options,'ub_spec_nonoise','ytitle','Energy (keV)'
options,'ub_spec_nonoise','xtitle','time (unitless)'
options,'ub_spec_nonoise','title','Simulated Microburst flux (normalized)!Cno noise'

store_data,'tmpcomb3',data=['ub_spec_nonoise','fitlineHR']


ylim,'tmpcomb3',0,1000
zlim,'ub_spec_nonoise',0.1,f0,1
loadct,39
;tplot,['ub_timevariation','ub_spec_nonoise']
tplot,['tmpcomb2','tmpcomb','tmpcomb3']

stop

;fluxN = flux
;for ee=0.,nenergies-1 do begin $
;  etmp = ee & $
;  noise = noise_scalefac*(randomu(etmp,ntsteps) - 1.) & $
;  noiseS = noise * sqrt(max(flux[*,*]) - flux[ee,*]) & $
;  fluxN[ee,*] += noiseS
;endfor
;flux = fluxN

;Now add in instrumental noise spectrum to this. This will only effect channels with low counts.  
noisetimesteps = tend_sec/cadence_newdetector
noisetimes = tend_sec*indgen(noisetimesteps)/(noisetimesteps-1)
noiseN = fltarr(nenergies,noisetimesteps)

fluxN = flux
for ee=0.,nenergies-1 do begin $
  etmp = ee & $
  noise = noise_scalefac*(randomu(etmp,noisetimesteps)) & $
  noiseS = noise * sqrt(max(flux[*,*]) - flux[ee,*]) & $
  noiseN[ee,*] = noiseS
endfor
;Increase the cadence of noise array to the full time cadence. 
noiseF = fltarr(nenergies,ntsteps)
for i=0,nenergies-1 do noiseF[i,*] = interpol(reform(noiseN[i,*]),noisetimes,times_sec)
fluxN = flux + noiseF

;stop
;RBSP_EFW> print,max(noisef)
;     0.109538
;RBSP_EFW> print,median(noisef)
;    0.0548391

flux = fluxN




;Add in contaminating microburst "noise" at sample rate cadence
;These are the sample rate bumps you see in real data that may correspond to 
;edge-detected microbursts?
;This "noise" is channel dependent. I'll set it at a fraction of the peak value 
;detected in each channel (noise_from_uB_fraction). NOTE: this value will need to be 
;abnormally high b/c the detectors integrate over many adjacent (1000x) energy bins and 
;the noise tends to disappear to unrealistically low values. For example, if you want the 
;noise in the realistic uB (say energy channel 500) to be 10% of the max value in that channel, 
;you would set noise_from_uB_fraction >> 0.1. 
;Just adjust so that things look real. 


noisetimesteps = tend_sec/cadence_newdetector
noisetimes = tend_sec*indgen(noisetimesteps)/(noisetimesteps-1)
noiseN = fltarr(nenergies,noisetimesteps)

for ee=0.,nenergies-1 do begin $
  etmp = ee & $
  noiseamp = noise_from_uB_fraction * max(flux[ee,*],/nan) & $
  tmp = randomu(etmp,noisetimesteps) & $ 
  noise = noiseamp*tmp & $
  noiseN[ee,*] = noise
endfor
;Increase the cadence of noise array to the full time cadence. 
noiseF = fltarr(nenergies,ntsteps)
for i=0,nenergies-1 do noiseF[i,*] = interpol(reform(noiseN[i,*]),noisetimes,times_sec)
;!p.multi = [0,0,2]
;plot,noiseN[20,*]
;plot,noiseF[20,*]
fluxN = flux + noiseF
;plot,fluxN[300,*]
;stop




;;***This is my check to see how the integrated flux profile changes.
;Eint_t = fltarr(ntsteps)
;for t=0,ntsteps-1 do Eint_t[t] = int_tabulated(energies_uB,flux[*,t])
;!p.multi = [0,0,1]



;;Plot flux spectra vs energy for all time steps
;nrows = ceil(sqrt(ntsteps))
;!p.multi = [0,nrows,nrows]
;ytitles = 'flux (t=' + strtrim(indgen(ntsteps),2)+')'
;for i=0,ntsteps-1 do plot,energies_uB,fluxN[*,i],yrange=[0,1],xtitle='Energy (keV)',ytitle=ytitles[i]

;stop

;fluxorig = flux

;noisesize = 0.1*fluxmax
;;
;;Add random noise to flux
;for t=0,ntsteps-1 do begin $
;  ttmp = t & $
;  noise = 0.2*(randomu(ttmp,nenergies)) & $
;  noise_mod = noise / flux[*,t] & $
;  flux[*,t] += noise_mod
;endfor
;
;!p.multi = [0,0,3]
;plot,fluxorig[*,100]
;plot,flux[*,100]
;plot,noise_mod
;stop

;;Make a plot of all times stacked on top of each other
;;--define colors
;deltac = 256/ntsteps - 1
;colorv = deltac*indgen(ntsteps)
;!p.multi = [0,0,1]
;plot,energies_uB,flux[*,0],yrange=[0,1],xtitle='Energy (keV)',ytitle='flux'
;for i=0,ntsteps-1 do oplot,energies_uB,flux[*,i],color=colorv[i]



;Plot a microburst spectrogram (Energy vs time with color indicating flux)
;****NOTE: NEED TO COMPARE THIS TO ACTUAL FIREBIRD MICROBURST FLUX.



store_data,'ub_spec_wnoise',data={x:tms,y:transpose(fluxN),v:energies_uB}
options,'ub_spec_wnoise','spec',1
options,'ub_spec_wnoise','ytitle','Energy (keV)'
options,'ub_spec_wnoise','xtitle','time (unitless)'
options,'ub_spec_wnoise','title','Simulated Microburst flux (normalized)'

store_data,'tmpcomb4',data=['ub_spec_wnoise','fitlineHR']


ylim,'tmpcomb4',0,1000
zlim,'ub_spec_wnoise',0.1,f0,1
loadct,39
tplot,['tmpcomb2','tmpcomb4']
;tplot,['ub_timevariation','ub_spec_nonoise','ub_spec_wnoise']

stop


;---------------------------------------------------------------------
;Now let's push this microburst into a simulated FIREBIRD detector
;---------------------------------------------------------------------


fb_channel_resp = dblarr(nenergies,nchannels)


;;Assume a Gaussian response for each FB channel. (This probably isn't correct)
;;dE = fbhig - fblow
;dE = replicate(50.,nchannels)
;for i=0,4 do fb_channel_resp[*,i] = exp(-1*(energies_uB - Ecenter[i])^2/dE[i]^2)


;Assume a square response for each FB channel - may be more accurate than a Gaussian
for i=0,nchannels-1 do begin $
  goo = where((energies_uB ge fblow[i]) and (energies_uB le fbhig[i])) & $
  fb_channel_resp[goo,i] = 1.
endfor
;;Plot imagined energy response for all FB channels
;!p.multi = [0,0,1]
;plot,energies_uB,fb_channel_resp[*,0],yrange=[0,1.2]
;for i=0,nchannels-1 do oplot,energies_uB,fb_channel_resp[*,i],color=50*i
;stop


;Multiply the FB channel response by the incident uB flux for each time.
Fch_integrated = dblarr(nchannels,ntsteps)

;Choose which array to use 
;fluxfin = flux_tprofile  ;normalized time profile 
;fluxfin = flux  ;flux with NO noise
fluxfin = fluxN  ;flux with noise added

;for t=0,ntsteps-1 do begin $
;  fluxtmp = reform(fluxfin[*,t]) & $
;  flux_after_fb = fltarr(nenergies,nchannels) & $
;  for i=0,nchannels-1 do flux_after_fb[*,i] = fluxtmp*fb_channel_resp[*,i] & $
;  ;Now integrate the FB fluxes for each channel over all energies.
;  ;Also divide out bin width
;  for i=0,nchannels-1 do Fch_integrated[i,t] = int_tabulated(energies_uB,flux_after_fb[*,i])/binwidth[i]
;endfor ;each time


;Final output [nchannels, ntsteps]
for t=0,ntsteps-1 do begin $
  fluxtmp = reform(fluxfin[*,t]) & $
  for n=0,nchannels-1 do begin & $
    flux_after_fb_tmp = fluxtmp * reform(fb_channel_resp[*,n]) & $
    Fch_integrated[n,t] = int_tabulated(energies_uB,flux_after_fb_tmp)/binwidth[n]
  endfor
endfor 


;;  ;plot output of FB detector at each time
;  !p.multi = [0,0,3]
;  plot,energies_uB,fluxtmp  ;total flux vs energy at current time
;  plot,energies_uB,fb_channel_resp[*,0],yrange=[0,2]  ;FB gain channels
;  for i=0,4 do oplot,energies_uB,fb_channel_resp[*,i],color=50*i
;  plot,energies_uB,flux_after_fb[*,0]
;  for i=0,4 do oplot,energies_uB,flux_after_fb[*,i],color=50*i
;stop


store_data,'uB_after_detection',tms,transpose(Fch_integrated)
store_data,'ub_spec_after_detection',data={x:tms,y:transpose(Fch_integrated),v:ecenter} ;spectral version
options,'ub_spec_after_detection','spec',1


;******************************
;Finally, downsample the ideal microburst to the cadence of the new detector. 

;cadence_newdetector = 0.05  ;sec (what do you want cadence of new detector to be?)
t0 = tms[0]
t1 = tms[n_elements(tms)-1]
nelem = (t1 - t0)/cadence_newdetector + 2
newtimes = dindgen(nelem)*cadence_newdetector + t0
;print,time_string(newtimes,prec=3)

tinterpol_mxn,'uB_after_detection',newtimes,/overwrite
tinterpol_mxn,'ub_spec_after_detection',newtimes,/overwrite


;options,'ub_spec_after_detection_interp','spec',1
;zlim,['ub_spec_after_detection','ub_spec_after_detection_interp'],0,20,0
;
;loadct,39
;tplot,['ub_spec_after_detection','ub_spec_after_detection_interp']



;*******FIX THIS???
;tshift = newtimes[0] - time_double('2014-01-01')





;Plot resultant uB after FIREBIRD detect

split_vec,'uB_after_detection',suffix='_'+chnameslow
chnames = 'uB_after_detection_'+chnameslow



;;ylim,'uB_after_detection_'+['220','283','384','520','721'],0,max(Fch_integrated)
;ylim,'uB_after_detection_220',0,200
;ylim,'uB_after_detection_283',0,150
;ylim,'uB_after_detection_384',0,80
;ylim,'uB_after_detection_520',0,25
;ylim,'uB_after_detection_721',0,4

;options,'uB_after_detection_220','colors',0
;options,'uB_after_detection_283','colors',50
;options,'uB_after_detection_384','colors',100
;options,'uB_after_detection_520','colors',150
;options,'uB_after_detection_721','colors',200

;
;RBSP_EFW> print,fblow
;      220.000      283.000      384.000      520.000      721.000
;RBSP_EFW> print,ecenter
;      251.500      333.500      452.000      620.500      853.000
;RBSP_EFW> print,fbhig
;      283.000      384.000      520.000      721.000      985.000

;tplot,'ub_spec_after_detection'


store_data,'tmpcomb5',data=['ub_spec_after_detection','fitlineHR']


zlim,['ub_spec_nonoise','ub_spec_wnoise','ub_spec_after_detection'],0.001,150,1
ylim,'tmpcomb5',0,1000,0
zlim,'tmpcomb5',0.01,20,1
;ylim,['ub_spec_nonoise','ub_spec_wnoise','ub_spec_after_detection'],0,1000,0
;tplot,['ub_spec_nonoise','ub_spec_wnoise','uB_after_detection_721','uB_after_detection_520','uB_after_detection_384','uB_after_detection_283','uB_after_detection_220']
loadct,39
tplot,['tmpcomb2','tmpcomb4','tmpcomb5']
stop



;*********************************************************
;When running the simulated microburst through the actual FIREBIRD channels:
;For each energy bin, compare the line plot from simulated to real microbursts. 
; *********************************************************
if fb_comparison eq 1 then begin 

  options,'uB_after_detection_???','color',250
  store_data,'tmp220',data=['ub_spec_after_detection_realuB_line_0','uB_after_detection_220']
  store_data,'tmp283',data=['ub_spec_after_detection_realuB_line_1','uB_after_detection_283']
  store_data,'tmp384',data=['ub_spec_after_detection_realuB_line_2','uB_after_detection_384']
  store_data,'tmp520',data=['ub_spec_after_detection_realuB_line_3','uB_after_detection_520']
  store_data,'tmp721',data=['ub_spec_after_detection_realuB_line_4','uB_after_detection_721']
  
  loadct,39
  ;tplot,[reverse(chnames)]
  tplot,['tmp721','tmp520','tmp384','tmp283','tmp220']
 stop

endif


;tplot,['ub_spec_wnoise','ub_spec_after_detection',reverse(chnames)]

;!p.charsize = 2
;!p.multi = [0,0,5]
;for i=0,4 do plot,Fch_integrated[4-i,*],xtitle='time',ytitle='FB channel = '+strtrim(fblow[4-i],2)



;;-----------------------------------------------------------------------------
;;Compare to the linefits of the realistic uB (from microburst_fit_slope.pro)
;;-----------------------------------------------------------------------------

;restore,'~/Desktop/code/Aaron/github/mission_routines/IMPAX/realistic_uB_linefits_20160830_2047-2048'
;tplot_restore,filename='~/Desktop/code/Aaron/github/mission_routines/realistic_uB_spec_20160830_2047-2048.tplot'
;;Note that times for the realistic spec have been shifted so that they start at 2014-01-01, just like
;;the simulated microburst here. 
;
;get_data,'ub_spec_after_detection_realuB',data=dd
;store_data,'ubREAL_after_detection_220',dd.x,dd.y[*,0]
;store_data,'ubREAL_after_detection_283',dd.x,dd.y[*,1]
;store_data,'ubREAL_after_detection_384',dd.x,dd.y[*,2]
;store_data,'ubREAL_after_detection_520',dd.x,dd.y[*,3]
;store_data,'ubREAL_after_detection_721',dd.x,dd.y[*,4]
;
;store_data,'ch1comb',data=['ubREAL_after_detection_220','uB_after_detection_220'] & options,'ch1comb','colors',[0,250]
;store_data,'ch2comb',data=['ubREAL_after_detection_283','uB_after_detection_283'] & options,'ch2comb','colors',[0,250]
;store_data,'ch3comb',data=['ubREAL_after_detection_384','uB_after_detection_384'] & options,'ch3comb','colors',[0,250]
;store_data,'ch4comb',data=['ubREAL_after_detection_520','uB_after_detection_520'] & options,'ch4comb','colors',[0,250]
;store_data,'ch5comb',data=['ubREAL_after_detection_721','uB_after_detection_721'] & options,'ch5comb','colors',[0,250]

;options,'ch1comb','ytitle','220-283!CkeV'
;options,'ch2comb','ytitle','283-384!CkeV'
;options,'ch3comb','ytitle','384-520!CkeV'
;options,'ch4comb','ytitle','520-721!CkeV'
;options,'ch5comb','ytitle','721-985!CkeV'

;zlim,['ub_spec_after_detection','ub_spec_after_detection_realuB'],0.01,10,1
;ylim,['ub_spec_after_detection','ub_spec_after_detection_realuB'],220,850,0
;options,['ub_spec_after_detection','ub_spec_after_detection_realuB'],'panel_size',2
;options,'ch?comb','panel_size',0.5
;;title = 'Simulated vs real uB-Real uB at 20160830_2047-2048-Simulated has J0='
;tplot,['ub_spec_after_detection_realuB','ub_spec_after_detection','ch5comb','ch4comb','ch3comb','ch2comb','ch1comb']

;-----------------------------------------------------------------------------



stop


;-------------------------------------------------------
;Save the data so I can load it up with microburst_fit_slope.pro 
tplot_save,'*',filename='~/Desktop/fb_ub_'+strtrim(nchannels,2)+'channel_cadence='+strtrim(cadence_newdetector,2)+'_recreation_of_'+datetime_real



;tplot_save,'*',filename='~/Desktop/fb_ub_5channel_cadence='+strtrim(cadence_newdetector,2)
;tplot_save,'*',filename='~/Desktop/fb_ub_10channel_cadence='+strtrim(cadence_newdetector,2)
;tplot_save,'*',filename='~/Desktop/fb_ub_5channel_cadence='+strtrim(cadence_newdetector,2)+'_recreation_of_20160830_2047-2048'
;tplot_save,'*',filename='~/Desktop/fb_ub_10channel_cadence='+strtrim(cadence_newdetector,2)+'_recreation_of_20160830_2047-2048'
;tplot_save,'*',filename='~/Desktop/fb_ub_20channel_cadence='+strtrim(cadence_newdetector,2)+'_recreation_of_20160830_2047-2048'
;tplot_save,'*',filename='~/Desktop/fb_ub_40channel_cadence='+strtrim(cadence_newdetector,2)+'_recreation_of_20160830_2047-2048'

;-------------------------------------------------------



;The more channels I add the more the output should converge on the highres value 

;keV = 220
;plot,indgen(ntsteps),fluxfin[keV,*]
;oplot,indgen(ntsteps),fch_integrated[0,*],color=250



end




;;Define deltaE (width) profile over time for each time step
;deltaE = 2*indgen(ntsteps) + 300.
;;***plot1
;plot,indgen(ntsteps),deltaE,xtitle='time',ytitle='deltaE',title='deltaE(t) functional form'
