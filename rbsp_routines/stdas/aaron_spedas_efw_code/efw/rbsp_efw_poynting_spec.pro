;Use this to produce Poynting flux spectrograms after running
;rbsp_efw_poynting_flux_crib


pro rbsp_efw_poynting_spec,Ename,Bname

rbsp_efw_init

rbx = 'rbspb_'

split_vec,Ename,suffix='_'+['p1','p2','p3']
split_vec,Bname,suffix='_'+['p1','p2','p3']

get_data,Ename,tt,dd
t0 = min(tt,/nan)
t1 = max(tt,/nan)

!p.charsize = 1.5
!x.margin = [5,5]


;minval = 1d-15
minvalperp = 0.00003
minvalpara = 0.00003

step = 0.1
step_overlap = 0.5


muo = 4d0*!DPI*1d-7     ; -Permeability of free space (N/A^2 or H/m)


;Choose timerange
;date = '2012-10-13'
;t0 = time_double(date + '/07:00')
;t1 = time_double(date + '/08:30')
;date = '2013-02-26'
;t0 = time_double(date + '/09:00')
;t1 = time_double(date + '/10:20')


ew1 = tsample(Ename+'_p1',[t0,t1],times=tx1)
ew2 = tsample(Ename+'_p2',[t0,t1],times=tx2)
ew3 = tsample(Ename+'_p3',[t0,t1],times=tx3)
efield = [[ew1],[ew2],[ew3]]
store_data,Ename+'_r',data={x:tx1,y:efield}
time = tx1

bw1 = tsample(Bname+'_p1',[t0,t1],times=tx1)
bw2 = tsample(Bname+'_p2',[t0,t1],times=tx2)
bw3 = tsample(Bname+'_p3',[t0,t1],times=tx3)
bfield = [[bw1],[bw2],[bw3]]
store_data,Bname+'_r',data={x:tx1,y:efield}


;pro rbsp_spec, tplot_var, $
;	tplot_var_spec=tplot_var_spec, $
;	npts=npts, n_ave=n_ave, $
;	tspec=tspec, spec=spec, freq=freq, df=df, $
;	nan_fill_gaps=nan_fill_gaps, $
;	median_subtract=median_subtract, median_width=median_width, median_val=median_val, $
;	mingap=mingap, $
;	verbose=verbose

stop
rbsp_spec,Ename+'_r',npts=32,freq=ff
rbsp_spec,Bname+'_r',npts=32,freq=ff
tplot,[Ename+'_r_SPEC',Bname+'_r_SPEC']



;pperp = tsample(rbx+'pflux_nospinaxis_perp',[t0,t1],times=tp1)
;ppara = tsample(rbx+'pflux_nospinaxis_para',[t0,t1],times=tp2)
;
;store_data,'pperp',data={x:tp1,y:pperp}
;store_data,'ppara',data={x:tp2,y:ppara}


;rbsp_spec,'pperp',npts=64,freq=ff
;rbsp_spec,'ppara',npts=64,freq=ff
;ylim,['pperp_SPEC','ppara_SPEC'],0.001,1,1
;tplot,['pperp','pperp_SPEC','ppara','ppara_SPEC']


rbsp_spec,rbx+'pflux_nospinaxis_perp',npts=64,freq=ff
rbsp_spec,rbx+'pflux_nospinaxis_para',npts=64,freq=ff
ylim,[rbx+'pflux_nospinaxis_p???_SPEC'],1,100,1
tplot,[rbx+'pflux_nospinaxis_perp',rbx+'pflux_nospinaxis_perp_SPEC',$
			 rbx+'pflux_nospinaxis_para',rbx+'pflux_nospinaxis_para_SPEC']


get_data,'ppara_SPEC',data=tmp
pparaw = tmp.y
get_data,'pperp_SPEC',data=tmp
pperpw = tmp.y



get_data,rbx+'pflux_nospinaxis_perp',data=d
pperp = d.y
get_data,rbx+'pflux_nospinaxis_para',data=d
ppara = d.y



;Separate positive and negative values
goop = where(ppara ge 0.)
goon = where(ppara lt 0.)
ppara_p = ppara
ppara_p[goon] = 0.
ppara_n = ppara
ppara_n[goop] = 0.

store_data,'ppara_p',d.x,ppara_p
store_data,'ppara_n',d.x,ppara_n


rbsp_spec,'ppara_p',npts=32
rbsp_spec,'ppara_n',npts=32

get_data,'ppara_p_SPEC',data=tmp
ppara_pw = tmp.y
get_data,'ppara_n_SPEC',data=tmp
ppara_nw = tmp.y


times = tmp.x
freq_bins = tmp.v
;ppara_pw = double(slide_spec(time,ppara_p,step,step_overlap))
;ppara_nw = double(slide_spec(time,ppara_n,step,step_overlap))



;----------------------------------
;set color scales for general plots
;----------------------------------


;cspperp = bytscl(pperpw,min=min(pperpw,/nan),max=max(pperpw,/nan))
;csppara = bytscl(pparaw,min=min(pparaw,/nan),max=max(pparaw,/nan))
;
;cspparap = bytscl(ppara_pw,min=min(ppara_pw,/nan),max=max(ppara_pw,/nan))
;cspparan = bytscl(ppara_nw,min=min(ppara_nw,/nan),max=max(ppara_nw,/nan))


;!p.multi = [0,0,3]
;stop
;plot_spec,Exw,times,freq_bins,/no_interp,/zlog;,/nocb
;plot_spec,Eyw,times,freq_bins,/no_interp,/zlog
;plot_spec,Ezw,times,freq_bins,/no_interp,/zlog
;
;plot_spec,Bxw,times,freq_bins,/no_interp,/zlog
;plot_spec,Byw,times,freq_bins,/no_interp,/zlog
;plot_spec,Bzw,times,freq_bins,/no_interp,/zlog


;ylim,'pperp_SPEC',0,35,0
;get_data,'pperp_SPEC',data=tmp
;times = tmp.x
;
;!p.multi = [0,0,2]
;plot_spec,pperpw,times,ff,/no_interp,/zlog,shades=cspperpw
;plot_spec,pparaw,times,ff,/no_interp,/zlog,shades=cspparaw






;--------------------------------------------
;create red, blue and puke-yellow colors for plots
;--------------------------------------------

loadct,39
tvlct,r,g,b,/get
r[3] = 218
g[3] = 165
b[3] = 32
r[1] = 178
g[1] = 34
b[1] = 34
r[2] = 72
g[2] = 61
b[2] = 139

modifyct,20,'mycolors',r,g,b
loadct,20

;-------------------
;remove small values
;-------------------

;Positive parallel values
pparabp = ppara_pw
csparap = fltarr(n_elements(times),n_elements(freq_bins))
tmpy = where(pparabp lt minvalpara^2)
if tmpy[0] ne -1 then pparabp[tmpy] = 0
if tmpy[0] ne -1 then csparap[tmpy] = 3


;Negative parallel values
pparabn = ppara_nw
csparan = fltarr(n_elements(times),n_elements(freq_bins))
tmpy = where(pparabn lt minvalpara^2)
if tmpy[0] ne -1 then pparabn[tmpy] = 0
if tmpy[0] ne -1 then csparan[tmpy] = 3


;Perpendicular values
pperpb = pperpw
csperp = fltarr(n_elements(times),n_elements(freq_bins))
tmpx = where(pperpb lt minvalperp^2)
if tmpx[0] ne -1 then pperpb[tmpx] = 0
if tmpx[0] ne -1 then csperp[tmpx] = 3




;---------------------------
;red - upwards Poynting flux
;---------------------------

tmpy = where(pparabp gt minvalpara^2)
if tmpy[0] ne -1 then pparabp[tmpy] = 1
if tmpy[0] ne -1 then csparap[tmpy] = 1

;------------------------------
;blue - downwards Poynting flux
;------------------------------

tmpy = where(pparabn gt minvalpara^2)
if tmpy[0] ne -1 then pparabn[tmpy] = 2
if tmpy[0] ne -1 then csparan[tmpy] = 2



;Combine the upwards and downwards parallel spectra
pparab = pparabp < pparabn
csparab = csparap < csparan

get_data,Ename+'_r_SPEC',data=ddd
store_data,'pparab',ddd.x,csparab,ddd.v

window,2,xsize=400,ysize=700
wset,2
!p.multi = [0,0,1]
plot_spec,pparab,t2,freq_bins,/no_interp,/nocb,shades=csparab,xtitle='time (msec)',ytitle='freq (kHz)',title='y-hat Poynting flux: blue=+y,red=-y'



end
