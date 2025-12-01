;Read in KiNETX waveform data and run my_deg_pol.pro. 
;
; Interesting Bernstein waves (Barium releases) at 594 sec and 629 sec
;


;Note that there was an issue with the VLF34 digital data and only analog data was available. 

;Weidong took the VLF12D data digitized at 128kHz and created a new data file sampled at 40 kHz that was sampled
;“as close as possible” (nearest neighbor) to the VLF34A data sampled at 40 ks/sec.

;Weidong's files are
;1. Full-flight VLF34A IDL savsets (VLF34A_counts_TMwords13579_40kHz_Flight52007.sav). The data is down-sampled to 40Khz (using only TM-words 1,3,5,7,9).
;2. Full-flight VLF12D IDL savsets (VLF12D_mvm_onVLF34ATimeTags_40kHz_Flight52007.sav). This data is the VLF12D data that is closest to the data point of VLF34A. 
;For convenience, the corresponding VLF12D time-tags are also included in this file.

;Additional information: I have checked the time of nearest VLF12D data point “shifted” to that of VLF34A. 
;The maximum of shifted amount is ~5.3 micro-seconds for the 1st burst time range (590 sec to 600 sec) 
;and is ~ 4.5 micro-seconds for the 2nd burst time range (625 sec to 635 sec). 
;So, for features near 10kHz or less, these shifted time amounts are small.

;--------------------------------------------------------------------------------------------------------------------------


.compile /Users/abrenema/Desktop/code/Aaron/github/spectral-coherence-phaselag-analysis/aaron_chaston_polarization.pro
rbsp_efw_init




;-----------------------------------
;Henry's despun (zonal, meridional) solution (DC, 8 kHz)
;-----------------------------------

fn = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/52007_full_Esoln.sav'
restore,fn

times = etime
wf12 = emer
wf34 = ezon
wfz = wf12*0.01

;-----------------------------------
;Steve's DC (8 kHz) calibrated data in SC coordinates
;-----------------------------------

fn = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/52007_TM2_S1-DCEfield_V12D_V34D_V56D_V15D_V62D_Volts_mvm.sav'
restore,fn

;times = etime
wf12 = V12D_MVM
wf34 = V34D_MVM
wfz = V56D_MVM



;-----------------------------------
;Load Weidong's modified VLF data
;-----------------------------------


path = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/efield_VLF/'
fn = 'VLF12D_mvm_onVLF34ATimeTags_40kHz_Flight52007.sav'
restore,path+fn

;DATA_VLF12D     FLOAT     = Array[33953740]
;FN              STRING    = 'VLF12D_mvm_onVLF34ATimeTags_40kHz_Flight52007.sav'
;NOTE1           STRING    = 'This set of Data_VLF12D here is obtained from original VLF12D data point which is closest to the data point in VLF34A. Time_VLF12D is the corresponding original VLF12D Time. The sampling'...
;PATH            STRING    = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/efield_VLF/'
;SAMPLERATE      STRING    = '40 kHz'
;TIME_VLF12D     DOUBLE    = Array[33953740]
;TIME_VLF34A     DOUBLE    = Array[33953740]


      
fn = 'VLF34A_counts_mvm_TMwords13579_40kHz_Flight52007.sav'
;fn = 'VLF34A_counts_TMwords13579_40kHz_Flight52007.sav'
restore,path+fn

wf12 = DATA_VLF12D
wf34 = DATA_VLF34A
times12 = TIME_VLF12D
times34 = TIME_VLF34A


;only have times starting t=0
goo = where(times12 gt 0)
wf12 = wf12[goo]
times12 = times12[goo]

goo = where(times34 gt 0)
wf34 = wf34[goo]
times34 = times34[goo]


;Weidong data only: Subtract off baseline for wf34 data (note: detrending not working well)
baseline = 512.  ;counts
wf34 = wf34 - baseline


;both are the same size
;n_elements(wf12)
;n_elements(wf34)

wfz = wf12*0.01

;Since the time differences are very small, as indicated by Weidong, let's adopt times12 for all components
;tdiff = times12 - times34
;plot,times12,tdiff
times = times12
;--------------------------------------------------------




;*********************
;*********************
;*********************
;get rid of NaN values in wf that screw up polarization analysis. 
;boo = where(finite(wf12) eq 0)
;;wf12[boo] = wf12[boo-1]
;wf12[boo] = 0
;boo = where(finite(wf34) eq 0)
;;wf34[boo] = wf34[boo-1]
;wf34[boo] = 0
;*********************
;*********************
;*********************
;get rid of NaN values in wf that screw up polarization analysis.
boo = where(finite(wf12) ne 0)
wf12 = wf12[boo]
wf34 = wf34[boo]
wfz = wfz[boo]
times = times[boo]

boo = where(finite(wf34) ne 0)
wf34 = wf34[boo]
wf12 = wf12[boo]
wfz = wfz[boo]
times = times[boo]



;ssec = 620.
;esec = 640.
ssec = 580.
esec = 640.


goo = where((times ge ssec) and (times le esec))





;Keep polarizations consistent before/after rocket flip maneuver
store_data,'wf12',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf12[goo]}
store_data,'wf34',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf34[goo]}
store_data,'wfz',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wfz[goo]}





;************

;sampfreq = 1/(times[1]-times[0])
sampfreq = 40000.
start = time_string(time_double('1970-01-01/00:00:00') + ssec)
totpoints = round((esec - ssec) * sampfreq)
;totpointsDC = round((esec - ssec) * sampfreqDC)
rotatefield = 0
pol_lmt = 0.2
pow_lmt = 1d-12
pow_lmt_typ = 1
corse = 1
outps = 0


npts = 1024
aaron_chaston_polarization,'wf12','wf34','wfz',start,totpoints,rotatefield,pol_lmt,pow_lmt,pow_lmt_typ,corse,sampfreq,outps,npts=npts



;Extract slices for plot 
get_data,'Power',data=pow
get_data,'Ellipticity',data=elip2
get_data,'Helicity',data=hlip2
get_data,'Degree$of$Polarization',data=pol2


save,pow,filename='~/Desktop/Bernstein_pow.sav'
save,elip2,filename='~/Desktop/Bernstein_elip.sav'
save,hlip2,filename='~/Desktop/Bernstein_hlip.sav'
save,pol2,filename='~/Desktop/Bernstein_degpol.sav'


;pow_lmt2 = 1d-11
pow_lmt2 = 1d-13
;pow_lmt2 = 1d-9
pow2 = pow
goo = where(pow.y le pow_lmt2)
elip2.y[goo] = !values.f_nan
hlip2.y[goo] = !values.f_nan
pol2.y[goo] = !values.f_nan
pow2.y[goo] = !values.f_nan

;tplot_options,'title','Bernstein at 595 sec - my_deg_pol_kinetx.pro'
tplot_options,'title','Bernstein at 629 sec - my_deg_pol_kinetx.pro'
store_data,'Ellipticity2',data=elip2,dlim=dlime,lim=lime
store_data,'Helicity2',data=hlip2,dlim=dlimh,lim=limh
store_data,'Degree$of$Polarization2',data=pol2,dlim=dlimp,lim=limp
zlim,'Ellipticity2',-1,1
zlim,'Helicity2',0,0.5
;zlim,'Power',1d-12,1d-11,1
zlim,'Power',1d-12,1d-9,1
zlim,'Degree$of$Polarization2',0.9,1
ylim,['Power','Ellipticity2','Helicity2','Degree$of$Polarization2'],2000,8000,0
tplot,['wf12','Power','Degree$of$Polarization2','Ellipticity2','Helicity2']

tlimit,'1970-01-01/00:09:52','1970-01-01/00:10:00'
tlimit,'1970-01-01/00:10:28.00','1970-01-01/00:10:32.000'


stop




;10:29.600 ;RH pol
;10:30.100 ;LH pol
;------------------------------------------------------------
;Make line plots of relevant quantities at selected times. 
;------------------------------------------------------------

sr_spec = 1/(pow.x[1]-pow.x[0])

;tz = 629.6  ;RH polarization
tz = 629.5  ;RH polarization
;tz = 629.65  ;LH polarization
tplot,['Power','Ellipticity','Degree$of$Polarization']
;timebar,tz

;numavg_sec = 0.2
;numavg = round(numavg_sec*sr_spec)
;print,numavg
numavg = 2

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
xr = [2000,8000]
plot,freqs,powfin,ytitle='Power',xtitle='freq (Hz)',xrange=xr,ylog=1,yrange=[1d-15,1d-9],title=title,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
plot,freqs,elipfin,ytitle='Ellipticity',xtitle='freq (Hz)',xrange=xr,yrange=[-1.1,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0,0],color=250
plot,freqs,polfin,ytitle='Deg of Polarization',xtitle='freq (Hz)',xrange=xr,yrange=[0,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0.7,0.7],color=250


stop


end
