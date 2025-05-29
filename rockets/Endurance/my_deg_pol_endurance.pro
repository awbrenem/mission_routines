;Read in Endurance waveform data and run my_deg_pol.pro 


.compile /Users/abrenema/Desktop/code/Aaron/github/spectral-coherence-phaselag-analysis/aaron_chaston_polarization.pro
rbsp_efw_init



;path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_DC/'
;fn = '47001_TM1_LFDSP_S5DCE_DCES5_calibrated.sav'
;restore,path+fn

;w12 = DV12_MVM
;w34 = DV34_MVM
;timesDC = times
;store_data,'w12',data={x:time_double('1970-01-01/00:00:00') + timesDC,y:w12}
;store_data,'w34',data={x:time_double('1970-01-01/00:00:00') + timesDC,y:w34}
;rbsp_detrend,'w12',0.001
;rbsp_detrend,'w34',0.001
;rbsp_detrend,'w12_smoothed',3.
;rbsp_detrend,'w34_smoothed',3.


;-----------------------------------
;Load calibrated data from Python (from save_data_as_netcdf.py)
;-----------------------------------


;NetCDF 
path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'
fn = 'efield_vlf_calibrated.netcdf'

ncid = NCDF_OPEN(path+fn)
E12id = NCDF_VARID(ncid, 'E12')
E34id = NCDF_VARID(ncid, 'E34')
timeid = NCDF_VARID(ncid, 'time')

;get a variable out of the data
NCDF_VARGET, ncid, E12id, E12
NCDF_VARGET, ncid, E34id, E34
NCDF_VARGET, ncid, timeid, time

;Close the netcdf file
NCDF_CLOSE, ncid

wf12 = E12
wf34 = E34
wfz = E34
wfz[*] = 0.0001
times = time

;First keep polarizations consistent before/after rocket flip maneuver
fliptime = 600
;if ssec gt fliptime then wf34 = -1*wf34
goo = where(times ge fliptime)
if goo[0] ne -1 then wf34[goo] = -1*wf34[goo] 

;Now keep only desired data (helps program to run much faster)
ssec = 120
esec = 170
;ssec = 550
;esec = 650
;ssec = 100
;esec = 900

goo = where((times ge ssec) and (times le esec))

;Keep polarizations consistent before/after rocket flip maneuver
store_data,'wf12',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf12[goo]}
store_data,'wf34',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf34[goo]}
store_data,'wfz',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wfz[goo]}



;;-----------------------------------
;;Load uncalibrated data from Steve
;;-----------------------------------
;
;path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'
;fn = '47001_TM1_LFDSP_S5_VLF_mvm.sav'
;restore,path+fn
;
;wf12 = DVLF12_MVM
;wf34 = DVLF34_MVM
;wfz = wf34 
;wfz[*] = 0.0001
;times = TVLF
;
;;**************
;;artificially offset the VLF12 data one point to the right to test for bad instrument timing
;;wf12 = shift(wf12,-1)
;;*************
;
;
;;ssec = 100
;;esec = 240
;ssec = 550
;esec = 650
;goo = where((times ge ssec) and (times le esec))
;
;
;;Keep polarizations consistent before/after rocket flip maneuver
;fliptime = 600
;if ssec gt fliptime then wf34 = -1*wf34
;store_data,'wf12',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf12[goo]}
;store_data,'wf34',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wf34[goo]}
;store_data,'wfz',data={x:time_double('1970-01-01/00:00:00') + times[goo],y:wfz[goo]}
;
;----------------------------------------------------------------------------------


sampfreq = 1/(times[1]-times[0])
start = time_string(time_double('1970-01-01/00:00:00') + ssec)
totpoints = round((esec - ssec) * sampfreq)
;totpointsDC = round((esec - ssec) * sampfreqDC)
rotatefield = 0
pol_lmt = 0.0
pow_lmt = 1d-12
pow_lmt_typ = 1
corse = 1
outps = 0


npts = 4096
aaron_chaston_polarization,'wf12','wf34','wfz',start,totpoints,rotatefield,pol_lmt,pow_lmt,pow_lmt_typ,corse,sampfreq,outps,npts=npts

stop



ylim,['Power','Degree$of$Polarization','Ellipticity','Helicity'],1000,9000,0
zlim,'Ellipticity',-1,1
zlim,'Power',1d-15,1d-11
zlim,'Degree$of$Polarization',0,1
tplot,['Power','Degree$of$Polarization','Ellipticity']


;Extract slices for plot 
get_data,'Power',data=pow
get_data,'Ellipticity',data=elip,dlim=dlime,lim=lime
get_data,'Degree$of$Polarization',data=pol,dlim=dlimp,lim=limp
get_data,'Helicity',data=hel

;Save data to IDL save file for loading in Python
times_pow = pow.x 
freqs_pow = pow.v
data_pow = pow.y

times_elip = elip.x
freqs_elip = elip.v
data_elip = elip.y

times_pol = pol.x
freqs_pol = pol.v
data_pol = pol.y

times_hel = hel.x
freqs_hel = hel.v
data_hel = hel.y

store_data,tnames(),/del

save,/variables,filename='~/Desktop/tst.sav'


stop



;Version of plot with low deg of pol values removed

pol_lmt2 = 0.5
elip2 = elip
pol2 = pol
goo = where(pol2.y le pol_lmt2)
elip2.y[goo] = !values.f_nan
pol2.y[goo] = !values.f_nan



store_data,'Ellipticity2',data=elip2,dlim=dlime,lim=lime
store_data,'Degree$of$Polarization2',data=pol2,dlim=dlimp,lim=limp
zlim,'Ellipticity2',-0.3,-0.3
ylim,['Power','Ellipticity2','Degree$of$Polarization2'],100,8000,0
;ylim,'Power',6000,15000,0
tplot,['Power','w12','w34','Ellipticity2','Degree$of$Polarization2']
;timebar,tz


stop

;zlim,'Degree$of$Polarization',0.7,1
;tplot,['Power','Ellipticity','Degree$of$Polarization']
;timebar,tz






;------------------------------------------------------------
;Make line plots of relevant quantities at selected times. 
;------------------------------------------------------------

sr_spec = 1/(pow.x[1]-pow.x[0])

tz = 805.
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
plot,freqs,powfin,ytitle='Power',xtitle='freq (Hz)',xrange=xr,ylog=1,yrange=[1d-15,1d-9],title=title,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
plot,freqs,elipfin,ytitle='Ellipticity',xtitle='freq (Hz)',xrange=xr,yrange=[-1.1,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0,0],color=250
plot,freqs,polfin,ytitle='Deg of Polarization',xtitle='freq (Hz)',xrange=xr,yrange=[0,1.1],ystyle=1,xticklen=1,xgridstyle=1,yticklen=1,ygridstyle=1
oplot,[0,12000],[0.7,0.7],color=250


stop


end
