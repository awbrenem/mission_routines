;Read in Endurance waveform data and run my_deg_pol.pro 
;

.compile /Users/abrenema/Desktop/code/Aaron/github/spectral-coherence-phaselag-analysis/aaron_chaston_polarization.pro
.compile /Users/abrenema/Desktop/code/Aaron/github/signal_analysis/plot_wavestuff.pro
rbsp_efw_init



;---------------------------------------------------------------
;----------------------------------
;Load uncalibrated B-field data 
;----------------------------------

path = '/Users/abrenema/Desktop/Research/Rocket_missions/BeamPIE/data/bfield_AC/'
fn = '52009_TM2_LFDSP_S1-DeltaBxyz_counts.sav'
restore,path+fn

bx = DMAGXCOUNTS   ;on 3-4 axis near V3
by = DMAGYCOUNTS   ;on 1-2 axis near V1
bz = DMAGZCOUNTS   ;on spin axis near V10 (aft)
btimes = tmag 

;Subtract baseline offset
bx = bx - 16384.
by = by - 16384.
bz = bz - 16384.

wf12 = by 
wf34 = bx 
wfz = bz
times = btimes

;Orient these so that the antenna coord system points along Bo.
;See:
;--boom_diagram.jpg
;This is a RH coord. system with Bo along the z-direction.
wf12 *= -1 ;x-hat
wfz *=  1  ;z-hat
wf34 *= 1  ;y-hat
;---------------------------------------------------------------
;---------------------------------------------------------------




;---------------------------------------------------------------
;-----------------------------------
;Load uncalibrated E-field data from Steve
;-----------------------------------

;path = '/Users/abrenema/Desktop/Research/Rocket_missions/BeamPIE/data/efield_AC/'
;fn = '52009_TM2_LFDSP_S1-VLF12_VLF12D_mvm.sav'
;restore,path+fn
;wf12 = DVLF_MVM
;times = tvlf

;fn = '52009_TM2_LFDSP_S1-VLF34_VLF34D_mvm.sav'
;restore,path+fn
;wf34 = DVLF_MVM

;fn = '52009_TM2_LFDSP_S1-VLF910_VLF910D_mvm.sav'
;restore,path+fn
;wfz = DVLF_MVM   ;V9-10


;;only have times starting t=0
;goo = where(times gt 0)
;wf12 = wf12[goo]
;wf34 = wf34[goo]
;wfz = wfz[goo]
;times = times[goo]

;Orient these so that the antenna coord system points along Bo. 
;See: 
;--BeamPIE-payload-sketch-Long_with_Magorient.png
;--BeamPIE-Cals-Spring-2023_v2.pdf

;wf12 *= -1 ;x-hat
;wfz *= -1  ;z-hat 
;wf34 *= 1  ;y-hat 
;---------------------------------------------------------------
;---------------------------------------------------------------


;Good for both Bw and Ew
sampfreq = 64000.

;By 200 sec the booms are deployed and system is aligned along Bo. 

;---------------------------------------------------------------------
;Analyze swooshes from 280-317 sec (5-10 kHz)

ssec = 260.
esec = 270.
;esec = 330.
too = where((times ge ssec) and (times le esec))


store_data,'wf12',data={x:time_double('1970-01-01/00:00:00') + times[too],y:wf12[too]}
store_data,'wf34',data={x:time_double('1970-01-01/00:00:00') + times[too],y:wf34[too]}
store_data,'wfz',data={x:time_double('1970-01-01/00:00:00') + times[too],y:wfz[too]}
store_data,'wf',data={x:time_double('1970-01-01/00:00:00') + times[too],y:[[wf12[too]],[wf34[too]],[wfz[too]]]}

plot_wavestuff,'wf',/spec,npts=256,/nodelete


ylim,['x_SPEC','y_SPEC','z_SPEC'],3000.,10000.
zlim,['x_SPEC','y_SPEC','z_SPEC'],1d-3,1d-2
tplot,['x_SPEC','y_SPEC','z_SPEC']



;************

totpoints = n_elements(where((times ge ssec) and (times le esec)))

start = time_string(time_double('1970-01-01/00:00:00') + ssec)
totpoints = round((esec - ssec) * sampfreq)
;totpointsDC = round((esec - ssec) * sampfreqDC)
rotatefield = 0
pol_lmt = 0.0
pow_lmt = 1d-12
pow_lmt_typ = 1
corse = 1
outps = 0


npts = 128.
aaron_chaston_polarization,'wf12','wf34','wfz',start,totpoints,rotatefield,pol_lmt,pow_lmt,pow_lmt_typ,corse,sampfreq,outps,npts=npts

ylim,['Power','Ellipticity','Helicity','Wavenormal$Angle','Degree$of$Polarization'],3000,10000
tplot,['Power','Ellipticity','Helicity','Wavenormal$Angle','Degree$of$Polarization']


;t0z = 260.
;t1z = 270.
;t0z2 = time_double('1970-01-01/00:00:00')+t0z
;t1z2 = time_double('1970-01-01/00:00:00')+t1z

;zlim,['x_SPEC','y_SPEC','z_SPEC'],1d-12,1d-7
;tplot,['x_SPEC','y_SPEC','z_SPEC']
;tlimit,t0z2,t1z2

;ylim,['Power','Wavenormal$Angle','Degree$of$Polarization','Ellipticity','Helicity'],0,30000,0
;zlim,'Ellipticity',-1,1
;tplot,['Power','Degree$of$Polarization','Ellipticity','Helicity','Wavenormal$Angle']
;tlimit,t0z2,t1z2

get_data,'Power',data=pow,dlim=dlimpow,lim=limppow
get_data,'Ellipticity',data=elip,dlim=dlime,lim=lime
get_data,'Helicity',data=hlip,dlim=dlimh,lim=limh
get_data,'Degree$of$Polarization',data=pol,dlim=dlimp,lim=limp
get_data,'Wavenormal$Angle',data=wna,dlim=dlimw,lim=limw

;Version of plot with low deg of pol values removed
pol_lmt2 = 0.7
elip2 = elip
hlip2 = hlip
pol2 = pol
wna2 = wna
goo = where(pol2.y le pol_lmt2)
elip2.y[goo] = !values.f_nan
hlip2.y[goo] = !values.f_nan
pol2.y[goo] = !values.f_nan
wna2.y[goo] = !values.f_nan

pow_lmt2 = 1d-9
pow2 = pow
goo = where(pow.y le pow_lmt2)
elip2.y[goo] = !values.f_nan
hlip2.y[goo] = !values.f_nan
pol2.y[goo] = !values.f_nan
wna2.y[goo] = !values.f_nan
pow2.y[goo] = !values.f_nan

store_data,'Power2',data=pow2,dlim=dlimpow,lim=limppow
store_data,'Wavenormal$Angle2',data=wna2,dlim=dlimw,lim=limw
store_data,'Ellipticity2',data=elip2,dlim=dlime,lim=lime
store_data,'Helicity2',data=hlip2,dlim=dlimh,lim=limh
store_data,'Degree$of$Polarization2',data=pol2,dlim=dlimp,lim=limp
zlim,'Ellipticity2',-1,1
zlim,'Helicity2',0,1
zlim,['Power','Power2'],1d-7,1d-6
ylim,['x_SPEC','y_SPEC','z_SPEC','Power','Power2','Ellipticity2','Helicity2','Degree$of$Polarization2','Wavenormal$Angle2'],3000,12000,0
;tplot_options,'title','BeamPIE AC magnetic fields!C t=261-269 sec'
tplot,['Power','Degree$of$Polarization2','Ellipticity2','Helicity2','Wavenormal$Angle2']


stop

zlim,'Degree$of$Polarization',0.7,1
tplot,['Power','Ellipticity','Degree$of$Polarization']
timebar,tz






;------------------------------------------------------------
;Make line plots of relevant quantities at selected times. 
;------------------------------------------------------------

sr_spec = 1/(pow.x[1]-pow.x[0])

tz = 266.88
tplot,['Power','Ellipticity','Degree$of$Polarization']
timebar,tz

;numavg_sec = 0.2
;numavg = round(numavg_sec*sr_spec)
;print,numavg
numavg = 5

title = string(numavg) + ' bins averaged starting at ' + strtrim(floor(tz),2) + ' sec'

goo = where(pow.x ge tz)
freqs = pow.v


powz = pow.y[goo[0]:goo[0]+numavg,*]
elipz = elip.y[goo[0]:goo[0]+numavg,*]
polz = pol.y[goo[0]:goo[0]+numavg,*]
powfin = fltarr(n_elements(freqs))
elipfin = fltarr(n_elements(freqs))
polfin = fltarr(n_elements(freqs))
for i=0,n_elements(freqs)-1 do begin $
  boo1 = where(finite(powz[*,i]) ne 0,c1) & $
  powfin[i] = total(powz[*,i],/nan)/c1 & $
  boo2 = where(finite(elipz[*,i]) ne 0,c2) & $
  elipfin[i] = total(elipz[*,i],/nan)/c2 & $
  boo3 = where(finite(polz[*,i]) ne 0,c3) & $
  polfin[i] = total(polz[*,i],/nan)/c3
;endfor



!p.charsize = 3
!p.multi = [0,0,3]
xr = [3000,12000]
plot,freqs,powfin,ytitle='Power',xtitle='freq (Hz)',xrange=xr,ylog=1,title=title,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
plot,freqs,elipfin,ytitle='Ellipticity',xtitle='freq (Hz)',xrange=xr,yrange=[-1.1,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0,0],color=250
plot,freqs,polfin,ytitle='Deg of Polarization',xtitle='freq (Hz)',xrange=xr,yrange=[0,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0.7,0.7],color=250


stop


end
