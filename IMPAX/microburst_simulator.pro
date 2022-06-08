;Create a microburst (with ~infinite resolution) and simulate its detection in a hypothetical detector as a function of energy and time.
;NOTE: one way to simulate a microburst is to imitate one detected on FIREBIRD. To do this first run microburst_fit_slope.pro to 
;get the best-fit line for a selected microburst. Use this to create the ideal microburst with the correct dispersion and spectrum.  

;After you've created the ideal microburst and run it through a hypothetical detector using this code, you can then 
;run microburst_fit_slope.pro with "option 2" to find a best-fit slope and error slopes. These two programs are thus a way to 
;test how good a hypothetical detector is at detecting microburst dispersion. 

;Finally, the test fit slope and error slopes can be used in my ray tracing routines to identify the region where the microbursts are created. 
;The smaller the error in slopes (by having a better detector) the more restricted the possible source region would be. 
;SEE microburst_fit_slope.pro and crib_raytrace_IMPAX.pro



;Microburst profile
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


;Noise: Two types of noise are added
;1) edge-detected microburst "noise". This is meant to imitate the existence of other random microburst. 
;The noise values scale with the peak flux (all times) of EACH energy channel. This noise is added to the 
;"infinite resolution" microburst (e.g. 1000 energies by 1000 times)
;2) Instrument noise. Energy channel independent. A set value with fluctuations. This is added after the simulated 
;microburst is run through the hypothetical detector. 




;***********************************
;***********************************
;Harlan email on IMPAX instrument AFIRE properties (from April 22, 2022):
;Emin = 125 keV
;nchannels = 64
;dE=14 keV --> dE/E of ~10% at the lowest bin and <1.5% at the top differential channel.


;Harlan email from 2021:
;NOTE: I'm not sure how Harlan is defining "flux units"
;FWHM is 12 keV for noise. (typically set threshold 3x noise)
;Can go even lower in energy.
;--mention that an even lower energy threshold would make dispersion more obvious (Saito)

;Max instrument noise at 220 keV (65 bins) is about 5 in flux units.
;Instrument noise is constant across bins at 0.05 in flux units and is only apparent in highest channels

;Estimate of ~31 counts/flux_unit

;***********************************



;*********************************************
;TODO: incorporate 10-50msec integration time that FIREBIRD uses for uB simulation.
;*********************************************




rbsp_efw_init
!p.charsize = 1.5



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
;Set detector cadence of the idealized detector (Harlan email)
cadence_newdetector = 20 ;msec
;---



;----------------------------------------------------------
;FLUX-DEPENDENT NOISE CAUSED BY EDGE-DETECTED MICROBURSTS
;----------------------------------------------------------

;These are the sample rate bumps you see in real data that may correspond to 
;edge-detected microbursts?
;This "noise" is channel dependent. I'll set it at a fraction of the peak value 
;detected in each channel (noise_from_uB_fraction). 
;e.g. noiseamp = noise_from_uB_fraction * max(flux[ee,*],/nan)
;
noise_from_uB_fraction = 2. 
;noise_from_uB_fraction = 0.

;NOTE: YOU MAY WANT TO SET noise_from_uB_fraction ARTIFICIALLY HIGH (i.e. > 1) FOR THE FOLLOWING REASON: The hypothetical detector 
;will have less than 1000 energy channels, and thus represents an integration over a certain number of energy channels of the
;idealized microburst. This will increase the signal/noise ratio. Set this value so that the end result (after detection in 
;hypothetical detector) looks similar to what you think a real microburst would look like. 

;Set max and min values (from 0-1) for the variation in random noise of contaminating microbursts RELATIVE to the 
;peak value in each energy bin (for all times) 
;For example, minn=0.6 and maxx = 0.8 means that the contaminating noise microbursts can range from 0.6 to 0.8
;times the peak microburst value in each energy bin.
maxx = 1.2
minn = 0.4
;---






;-----------------------------------------
;Set idealized microburst energy profile
;-----------------------------------------

;energy profile ("infinite" resolution uB) --> exp(-(E-Eo)^epow/dE^epow)
nenergies = 1000. ;number of microburst energy steps
emin_keV = 0. 
emax_keV = 1000.
energies_uB = (emax_keV - emin_keV)*indgen(nenergies)/(nenergies-1) + emin_keV   ;creates a 1000 element array
Eo = 0. ;Center energy for Gaussian
dE = 75. ;Width of energy Gaussian curve (keV)
epow = 1.  ;typically=1     flux[*,i] *= exp(-1*(energies_uB - Eo)^epow/dE^epow)



;-----------------------------------------
;Set idealized microburst time profile
;-----------------------------------------

;time profile  (exp(-(t-to)^2/dt^2))
ntsteps = 1000. ;number of microburst time steps for "infinite" resolution uB
tend_sec = 1.15  ;total time duration considered in seconds (microburst contained within this)
;dt = 0.2  ;Width (sec) of the exponential falloff
dt = 0.1  ;Width (sec) of the exponential falloff
to = 0.56  ;time offset for the zeroth bin (sec). Choose to align with actual microburst from microburst_fit_slope.pro
t_dispersion = 0.134 ;sec  - delta time b/t lowest and highest energy channels from best fit line (emin_keV; emax_keV)





fb_comparison = 0.



;-----------------------------------------------------------------
;;Define FIREBIRD detector channels
;-----------------------------------------------------------------



;;Values from Crew16 for the collimated detector on FU4
;fblow = [220.,283.,384.,520.,721.]
;fbhig = [283.,384.,520.,721.,985.]
;fb_comparison = 1.


;;Stupidly fat channels
;fblow = [220.,520.]
;fbhig = [520.,721.]


;;Hypothetical 10 channels
;fblow = [220.,251.,283.,333.,384.,452.,520.,620.,721.,853.]
;fbhig = [251.,283.,333.,384.,452.,520.,620.,721.,853.,985.]



;;Harlan channels from Apr 22, 2022 email
emaxtmp = 1000.  ;keV
emintmp = 40.    ;keV
nch = 6.
fblow = (emaxtmp - emintmp)*indgen(nch)/(nch-1) + emintmp
fbhig = shift(fblow,-1)
deltae = fbhig[0] - fblow[0]
fbhig[n_elements(fbhig)-1] = emaxtmp + deltae



;nch = 5.
;fblow = (1000 - 200.)*indgen(nch)/(nch-1) + 220.
;fbhig = shift(fblow,-1)
;fbhig[n_elements(fbhig)-1] = 1040.



;----------------------------------------------------







nchannels = n_elements(fblow)
binwidth = fbhig - fblow
Ecenter = (fblow + fbhig)/2.

;String channel names
chnameslow = strtrim((floor(fblow)),2)
chnameshig = strtrim((floor(fbhig)),2)

deltae = fbhig[0] - fblow[0]


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
;options,'ub_timevariation','title','Simulated Microburst time variation (normalized)'

options,'fitlineHR','colors',0
store_data,'tmpcomb',data=['ub_timevariation','fitlineHR']
store_data,'tmpcomb2',data=['ub_spec_after_detection_realuB','fitlineHR']

options,'tmpcomb','ytitle','Energy(keV)!Csimulated uB!Cdispersion only'
options,'tmpcomb2','ytitle','Energy(keV)!Creal uB'


ylim,['tmpcomb','tmpcomb2'],0,1000
zlim,'ub_timevariation',0,1,0
zlim,'ub_spec_after_detection_realuB',0.01,20,1
loadct,39
;Plot the normalized dispersive time signature. NOTE: don't expect the magnitude to match real microburst at this point.
;***check to see that the dispersion is similar to real microburst 
tplot,['tmpcomb2','tmpcomb']



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
;options,'ub_spec_nonoise','title','Simulated Microburst flux (normalized)!Cno noise'

store_data,'tmpcomb3',data=['ub_spec_nonoise','fitlineHR']
options,'tmpcomb3','ytitle','Energy(keV)!Csimulated uB!CSpectral falloff!Capplied'

ylim,'tmpcomb3',0,1000
zlim,'ub_spec_nonoise',0.1,f0,1
loadct,39
tplot,['tmpcomb2','tmpcomb','tmpcomb3']

stop





;-------------------------------------------------------------------
;Add in contaminating microburst "noise" at sample rate cadence
;-------------------------------------------------------------------

;These are the sample rate bumps you see in real data that may correspond to 
;edge-detected microbursts?
;This "noise" is channel dependent. I'll set it at a fraction of the peak value 
;detected in each channel (noise_from_uB_fraction). NOTE: this value will need to be 
;abnormally high b/c the detectors integrate over many adjacent (1000x) energy bins and 
;the noise tends to disappear to unrealistically low values. For example, if you want the 
;noise in the realistic uB (say energy channel 500) to be 10% of the max value in that channel, 
;you would set noise_from_uB_fraction >> 0.1. 
;Just adjust so that things look real. 


noisetimesteps = tend_sec/(cadence_newdetector/1000.)
noisetimes = tend_sec*indgen(noisetimesteps)/(noisetimesteps-1)
noiseN = fltarr(nenergies,noisetimesteps)

for ee=0.,nenergies-1 do begin
  etmp = ee
  noiseamp = noise_from_uB_fraction * max(flux[ee,*],/nan)
  tmp = (maxx - minn) * randomu(etmp,noisetimesteps) + minn
  noise = noiseamp*tmp
  noiseN[ee,*] = noise
endfor
;Increase the cadence of noise array to the full time cadence. 
noiseF = fltarr(nenergies,ntsteps)
for i=0,nenergies-1 do noiseF[i,*] = interpol(reform(noiseN[i,*]),noisetimes,times_sec)
fluxN = flux + noiseF

;These variables are used in microburst_fit_slope.pro to determine noise level at which to ignore datapoints.
noise_max = max(noiseF)
noise_med = median(noiseF)
noise_mean = mean(noiseF)




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



;Plot a microburst spectrogram (Energy vs time with color indicating flux)
;****NOTE: NEED TO COMPARE THIS TO ACTUAL FIREBIRD MICROBURST FLUX.



store_data,'ub_spec_wnoise',data={x:tms,y:transpose(fluxN),v:energies_uB}
options,'ub_spec_wnoise','spec',1
options,'ub_spec_wnoise','ytitle','Energy (keV)'
options,'ub_spec_wnoise','xtitle','time (unitless)'
;options,'ub_spec_wnoise','title','Simulated Microburst flux (normalized)'

store_data,'tmpcomb4',data=['ub_spec_wnoise','fitlineHR']
options,'tmpcomb4','ytitle','Energy(keV)!Csimulated uB!Cnoise added'


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


;Flux array input is [nenergies, ntimes]  - idealized uB   e.g. [1000, 1000]
;fb_channel_response is [nenergies, nchannels]  e.g. [1000, 64] (step function array)
;Final output after detection by hypothetical detector is [nchannels, ntsteps]  ;e.g. [64, 1000]

for t=0,ntsteps-1 do begin
  fluxtmp = reform(fluxfin[*,t])  ;slice of energies
  for n=0,nchannels-1 do begin
    flux_after_fb_tmp = fluxtmp * reform(fb_channel_resp[*,n])
    ;Fch_integrated[n,t] = int_tabulated(energies_uB,flux_after_fb_tmp)/binwidth[n]
    Fch_integrated[n,t] = total(flux_after_fb_tmp)/binwidth[n]
  endfor
endfor


;loadct,39
;plot, Fch_integrated[*,0] * binwidth[0]


;;  ;plot output of FB detector at each time
;  !p.multi = [0,0,3]
;  plot,energies_uB,fluxtmp  ;total flux vs energy at current time
;  plot,energies_uB,fb_channel_resp[*,0],yrange=[0,2]  ;FB gain channels
;  for i=0,4 do oplot,energies_uB,fb_channel_resp[*,i],color=50*i
;  plot,energies_uB,flux_after_fb[*,0]
;  for i=0,4 do oplot,energies_uB,flux_after_fb[*,i],color=50*i
;stop


store_data,'ub_after_detection',tms,transpose(Fch_integrated)
store_data,'ub_spec_after_detection',data={x:tms,y:transpose(Fch_integrated),v:ecenter} ;spectral version
options,'ub_spec_after_detection','spec',1


;******************************
;Downsample the ideal microburst to the cadence of the new detector. 

t0 = tms[0]
t1 = tms[n_elements(tms)-1]
nelem = (t1 - t0)/(cadence_newdetector/1000.) + 2
newtimes = dindgen(nelem)*(cadence_newdetector/1000.) + t0
;print,time_string(newtimes,prec=3)

tinterpol_mxn,'ub_after_detection',newtimes,/overwrite
tinterpol_mxn,'ub_spec_after_detection',newtimes,/overwrite


;---------------------------------------------------
;Finally, add noise to the final values 
noise_from_uB_fraction = 2.

get_data,'ub_spec_after_detection',data=d,dlim=dlim,lim=lim


;;Scale the instrument noise values (nominally from 0-1) to this range

;Harlan has indicated that the noise level in "flux units" is 0.05. 
;I think he may mean 0.05/cm2-sr. So, if we divide this by dE then we get 1/cm2-sr-keV (time independent)
;---NO IDEA IF THIS IS RIGHT
noise_flux = 0.05/deltae  ;noise level in flux units
;noise_flux = 0.



noiseN = fltarr(n_elements(d.v),n_elements(d.x))

maxflux= max(d.y[*,*],/nan)

for ee=0.,n_elements(d.v)-1 do begin
  etmp = ee
  noise = randomu(etmp,n_elements(d.x))
  noise = 2.* (noise - median(noise))  ;+/- values
  noiseS = noise * noise_flux
  noiseN[ee,*] = noiseS
endfor

store_data,'ub_spec_after_detection',d.x,d.y + transpose(noiseN),d.v,dlim=dlim,lim=lim
store_data,'ub_after_detection',d.x,d.y + transpose(noiseN)

;store_data,'ub_spec_after_detectionN',d.x,d.y + transpose(noiseN),d.v,dlim=dlim,lim=lim
;tplot,['ub_spec_after_detectionN','ub_spec_after_detection']
;get_data,'ub_spec_after_detection',data=d1
;get_data,'ub_spec_after_detectionN',data=d2
;loadct,39
;plot,d1.y[*,40]
;oplot,d2.y[*,40],color=250















;options,'ub_spec_after_detection_interp','spec',1
;zlim,['ub_spec_after_detection','ub_spec_after_detection_interp'],0,20,0
;
;loadct,39
;tplot,['ub_spec_after_detection','ub_spec_after_detection_interp']



;*******FIX THIS???
;tshift = newtimes[0] - time_double('2014-01-01')





;Plot resultant uB after FIREBIRD detect

split_vec,'ub_after_detection',suffix='_'+chnameslow
chnames = 'ub_after_detection_'+chnameslow



;stop


;;ylim,'ub_after_detection_'+['220','283','384','520','721'],0,max(Fch_integrated)
;ylim,'ub_after_detection_220',0,200
;ylim,'ub_after_detection_283',0,150
;ylim,'ub_after_detection_384',0,80
;ylim,'ub_after_detection_520',0,25
;ylim,'ub_after_detection_721',0,4

;options,'ub_after_detection_220','colors',0
;options,'ub_after_detection_283','colors',50
;options,'ub_after_detection_384','colors',100
;options,'ub_after_detection_520','colors',150
;options,'ub_after_detection_721','colors',200

;
;RBSP_EFW> print,fblow
;      220.000      283.000      384.000      520.000      721.000
;RBSP_EFW> print,ecenter
;      251.500      333.500      452.000      620.500      853.000
;RBSP_EFW> print,fbhig
;      283.000      384.000      520.000      721.000      985.000

;tplot,'ub_spec_after_detection'


store_data,'tmpcomb5',data=['ub_spec_after_detection','fitlineHR']
options,'tmpcomb5','ytitle','Energy(keV)!Csimulated uB!Cafter detection!Cwith hypothetical!Cinstrument'


zlim,['ub_spec_nonoise','ub_spec_wnoise','ub_spec_after_detection'],0.001,150,1
ylim,'tmpcomb5',0,1000,0
zlim,'tmpcomb5',0.01,20,1
;ylim,['ub_spec_nonoise','ub_spec_wnoise','ub_spec_after_detection'],0,1000,0
;tplot,['ub_spec_nonoise','ub_spec_wnoise','uB_after_detection_721','uB_after_detection_520','uB_after_detection_384','uB_after_detection_283','uB_after_detection_220']
loadct,39
tplot,['tmpcomb2','tmpcomb4','tmpcomb5']
stop


;******TESTING*********
get_data,'ub_spec_after_detection',data=dd
loadct,39
plot,dd.y[*,4]


;150 counts for 97 keV channel and 64 bins

;93 counts for 141 keV channel and 8 bins
;80 counts for 141 keV channel and 64 bins

;*********************************************************
;When running the simulated microburst through the actual FIREBIRD channels:
;For each energy bin, compare the line plot from simulated to real microbursts. 
; *********************************************************
if fb_comparison eq 1 then begin 

  options,'ub_after_detection_???','color',250
  store_data,'tmp220',data=['ub_spec_after_detection_realuB_line_0','ub_after_detection_220']
  store_data,'tmp283',data=['ub_spec_after_detection_realuB_line_1','ub_after_detection_283']
  store_data,'tmp384',data=['ub_spec_after_detection_realuB_line_2','ub_after_detection_384']
  store_data,'tmp520',data=['ub_spec_after_detection_realuB_line_3','ub_after_detection_520']
  store_data,'tmp721',data=['ub_spec_after_detection_realuB_line_4','ub_after_detection_721']
  
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



;stop


;-------------------------------------------------------
;Save the data so I can load it up with microburst_fit_slope.pro 
tplot_save,'*',filename='~/Desktop/fb_ub_'+strtrim(nchannels,2)+'channel_cadence='+strtrim(floor(cadence_newdetector),2)+'msec_Emin='+strtrim(floor(ecenter[0]),2)+'_Emax='+strtrim(floor(max(ecenter)),2)+'_recreation_of_'+datetime_real


savevars = [fblow,fbhig,nchannels,noise_max,noise_med,noise_mean] 
fnn = '~/Desktop/fb_ub_'+strtrim(nchannels,2)+'channel_cadence='+strtrim(floor(cadence_newdetector),2)+'msec_Emin='+strtrim(floor(ecenter[0]),2)+'_Emax='+strtrim(floor(max(ecenter)),2)+'_recreation_of_'+datetime_real + '.sav'
save,savevars,filename=fnn,/variables


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
