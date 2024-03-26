;Read in Endurance waveform data and run my_deg_pol.pro 
;



rbsp_efw_init

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'
fn = '47001_TM1_LFDSP_S5_VLF_mvm.sav'
restore,path+fn

.compile /Users/abrenema/Desktop/code/Aaron/github/spectral-coherence-phaselag-analysis/aaron_chaston_polarization.pro

wf12 = DVLF12_MVM
wf34 = DVLF34_MVM 
wfz = wf34 
wfz[*] = 0.0001

times = TVLF


t0 = time_double('1970-01-01/00:00:00')


;goo = where((times ge 150) and (times le 160))
ssec = 840
esec = 900
goo = where((times ge ssec) and (times le esec))


;Keep polarizations consistent before/after rocket flip maneuver
fliptime = 600
if ssec gt fliptime then wf34 = -1*wf34
store_data,'wf12',data={x:t0 + times[goo],y:wf12[goo]}
store_data,'wf34',data={x:t0 + times[goo],y:wf34[goo]}
store_data,'wfz',data={x:t0 + times[goo],y:wfz[goo]}



;************


sampfreq = 1/(times[1]-times[0])

start = time_string(time_double(t0) + ssec)
totpoints = round((esec - ssec) * sampfreq)
rotatefield = 0
;pol_lmt = 0.5
pol_lmt = 0.0
;pow_lmt = 1d-7
pow_lmt = 1d-12
pow_lmt_typ = 1
corse = 1
outps = 0


npts = 1024
;npts = 512
aaron_chaston_polarization,'wf12','wf34','wfz',start,totpoints,rotatefield,pol_lmt,pow_lmt,pow_lmt_typ,corse,sampfreq,outps,npts=npts





ylim,'*',100,8000,0
zlim,'Ellipticity',-1,1
tplot,['Power','Degree$of$Polarization','Ellipticity','Helicity']

tlimit,'1970-01-01/00:14:30','1970-01-01/00:15:00'



;Extract slices for plot 
get_data,'Power',data=pow
get_data,'Ellipticity',data=elip,dlim=dlime,lim=lime
get_data,'Degree$of$Polarization',data=pol,dlim=dlimp,lim=limp






sr_spec = 1/(pow.x[1]-pow.x[0])

tz = 790. 
tplot,['Power','Ellipticity','Degree$of$Polarization']
timebar,tz

;numavg_sec = 0.2
;numavg = round(numavg_sec*sr_spec)
;print,numavg
numavg = 10

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
xr = [4000,8000]
plot,freqs,powfin,ytitle='Power',xtitle='freq (Hz)',xrange=xr,ylog=1,yrange=[1d-15,1d-11],title=title,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
plot,freqs,elipfin,ytitle='Ellipticity',xtitle='freq (Hz)',xrange=xr,yrange=[-1.1,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0,0],color=250
plot,freqs,polfin,ytitle='Deg of Polarization',xtitle='freq (Hz)',xrange=xr,yrange=[0,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0.7,0.7],color=250



;Version of plot with low deg of pol values removed

pol_lmt2 = 0.7

elip2 = elip
pol2 = pol
goo = where(pol2.y le pol_lmt2)
elip2.y[goo] = !values.f_nan
pol2.y[goo] = !values.f_nan

store_data,'Ellipticity2',data=elip2,dlim=dlime,lim=lime
store_data,'Degree$of$Polarization2',data=pol2,dlim=dlimp,lim=limp
tplot,['Power','Ellipticity2','Degree$of$Polarization2']
timebar,tz


zlim,'Degree$of$Polarization',0.7,1
tplot,['Power','Ellipticity','Degree$of$Polarization']
timebar,tz






;my_deg_pol,times[goo],wf12[goo],wf34[goo],wfz[goo]


store_data,'wfx',times[goo],wf12[goo]

ylim,'*',3000,9000,0

tplot,[36,37,39]








;openr,lun,'~/Desktop/file.txt',/get_lun
;vals = 
;while not eof(lun) do begin
;  readf,lun,jnk 
;endwhile


close,lun 


