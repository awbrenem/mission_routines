;Use READ_MYCDF to easily read multiple CDF files and specific or all variables
;within them.
;See https://spdf.gsfc.nasa.gov/CDAWlib.html for documentation of the CDAWlib
;package

;**Crib sheet designed to be run by copy/paste


rbsp_efw_init

;Must include this to use the CDFWlib data package
@compile_cdaweb


;Load L3 SWEAP files (spp_swp_spc_l3i_20181106_v08.cdf)
cdfnames = dialog_pickfile(/multiple_files,path='~/Desktop/Research/other/Stuff_for_other_people/Wygant_John/SWEAP_SPC_L3_CDF/')


;**Read in all variables...can take a while
r = read_mycdf('',/all,cdfnames)
;**Read in only select variables
;r = read_mycdf(['np_fit','vp_fit_RTN','sc_pos_HCI','sc_vel_HCI'],cdfnames)


;See contents of structure you've just read in
help,r,/st


;Get the names of all the variables within the structure
varnames = tag_names(r)
varnames = varnames[where(varnames ne 'EPOCH')]


;get the times
t = r.epoch.dat


;**************************************************************
;;Use this if  CDF epoch times are milliseconds since 1-Jan-0000
;;(like the EFW CDFs)
;t = real_part(t)
;cdf_epoch,1000d*t,yr, mo, dy, hr, mn, sc, milli, /BREAK
;tunix = strarr(n_elements(yr))
;for i=0L,n_elements(tunix)-1 do tunix[i] = strtrim(yr[i],2)+'-'+strtrim(mo[i],2)+'-'+$
;strtrim(dy[i],2)+'/'+strtrim(hr[i],2)+':'+strtrim(mn[i],2)+':'+$
;strtrim(sc[i],2)+'.'+strtrim(milli[i],2)
;tunix = time_double(tunix)

;**************************************************************
;Use this if CDF times are in TT2000 times
;(like Justin Kasper's PSP CDFs)
t = long64(t)
CDF_TT2000, t, yr, mo, dy, hr, mn, sc, milli, /BREAK
tunix = strarr(n_elements(yr))
yr = strtrim(floor(yr),2) & mo = strtrim(floor(mo),2) & dy = strtrim(floor(dy),2) & hr = strtrim(floor(hr),2) & mn = strtrim(floor(mn),2) & sc = strtrim(floor(sc),2) & milli = strtrim(floor(milli),2)

;Pad with zeros
goo = where(mo lt 10) & if goo[0] ne -1 then mo[goo] = '0'+mo[goo]
goo = where(dy lt 10) & if goo[0] ne -1 then dy[goo] = '0'+dy[goo]
goo = where(hr lt 10) & if goo[0] ne -1 then hr[goo] = '0'+hr[goo]
goo = where(mn lt 10) & if goo[0] ne -1 then mn[goo] = '0'+mn[goo]
goo = where(sc lt 10) & if goo[0] ne -1 then sc[goo] = '0'+sc[goo]
goo = where(milli lt 10) & if goo[0] ne -1 then milli[goo] = '00'+milli[goo]
goo = where((milli ge 10) and (milli lt 100)) & if goo[0] ne -1 then milli[goo] = '0'+milli[goo]

tunix = yr+'-'+mo+'-'+$
dy+'/'+hr+':'+mn+':'+$
sc+'.'+milli
tunix = time_double(tunix)
;***************************************************


;now grab each quantity and store as tplot variable.
;This loop checks to see if the "data" is actually timeseries data
;by comparing the size of the data array to the time array
for j=0,n_elements(varnames)-1 do begin $
  strtmp = 'dat = r.'+varnames[j] + '.dat' & $
  void = execute(strtmp) & $
  sizetmp = size(dat) & $
  if sizetmp[0] eq 1 then sz = sizetmp[1] else sz = sizetmp[2] & $
  sizetst = n_elements(tunix) eq sz & $
  if sizetst then store_data,varnames[j],tunix,reform(transpose(dat))
endfor




;The SWEAP variables from above that are needed are:
;ylim,'VP_MOMENT_SC',-500,500
get_data,'VP_MOMENT_SC',data=d,dlim=dlim,lim=lim
goo = where(d.y eq -1.00000e+31)
if goo[0] ne -1 then d.y[goo] = !values.f_nan
store_data,'VP_MOMENT_SC',data=d,dlim=dlim,lim=lim

get_data,'VP_MOMENT_RTN',data=d,dlim=dlim,lim=lim
goo = where(d.y eq -1.00000e+31)
if goo[0] ne -1 then d.y[goo] = !values.f_nan
store_data,'VP_MOMENT_RTN',data=d,dlim=dlim,lim=lim

get_data,'WP_MOMENT',data=d,dlim=dlim,lim=lim
goo = where(d.y eq -1.00000e+31)
if goo[0] ne -1 then d.y[goo] = !values.f_nan
store_data,'WP_MOMENT',data=d,dlim=dlim,lim=lim

get_data,'NP_MOMENT',data=d,dlim=dlim,lim=lim
goo = where(d.y eq -1.00000e+31)
if goo[0] ne -1 then d.y[goo] = !values.f_nan
store_data,'NP_MOMENT',data=d,dlim=dlim,lim=lim


vars = ['VP_MOMENT_SC','VP_MOMENT_RTN',$
'CARR_LATITUDE','CARR_LONGITUDE',$
'WP_MOMENT','NP_MOMENT',$
'SC_POS_HCI','SC_VEL_HCI']

tplot,vars

tplot_save,vars,filename='~/Desktop/psp_sweap_vars_nov18'


;store_data,'density',tunix,r.np_fit.dat
;ylim,'density',1,10000,1

;vp_rtn = r.vp_fit_RTN.dat
;store_data,'vp_rtn',tunix,transpose(vp_rtn)
;store_data,'vmag',tunix,reform(sqrt(vp_rtn[0,*]^2 + vp_rtn[1,*]^2 + vp_rtn[2,*]^2))
;ylim,'vp_rtn',-1000,1000
;ylim,'vmag',0,10000

;pos = r.sc_pos_HCI.dat
;vel = r.sc_vel_HCI.dat
;store_data,'pos_hci',tunix,transpose(pos)
;store_data,'vel_hci',tunix,transpose(vel)
;store_data,'vsc_mag',tunix,reform(sqrt(vel[0,*]^2 + vel[1,*]^2 + vel[2,*]^2))


;tplot,['density','vp_rtn','vmag','pos_hci','vel_hci','vsc_mag']



end
