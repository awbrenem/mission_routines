;Cribsheet for loading various PSP L2 SWEAP data quantities.

;*****************************************************************
;********HOW TO FIND SWEAP FILES**********************************
;*****************************************************************

;1) New L2 and L3 files seem to be at

;LOW RESOLUTION MOMENTS (~15 sec cadence)
;-SPI L3:     http://sweap.cfa.harvard.edu/data/sci/sweap/spi/L3/spi_sf00/
;               e.g. psp_swp_spi_sf00_L3_mom_INST_20200128_v01.cdf
;             http://sweap.cfa.harvard.edu/data/sci/sweap/spi/L3/spi_sf01/
;               e.g. psp_swp_spi_sf01_L3_mom_INST_20200128_v01.cdf
;             http://sweap.cfa.harvard.edu/data/sci/sweap/spi/L3/spi_sf0a/
;               e.g. psp_swp_spi_sf0a_L3_mom_INST_20200128_v01.cdf
;All of the above 3 files seem basically the same and contain the moments
;               DENS
;               VEL
;               T_TENSOR
;               TEMP
;               MAGF_SC
;               MAGF_INST
;               EFLUX_VS_ENERGY
;               EFLUX_VS_THETA
;               EFLUX_VS_PHI
;
;********USE THESE!!!!!!!!!!!*********************
;********FOR CONSISTENCY WITH WHAT I'VE PREVIOUSLY SENT TO WYGANT
;*************************************************
;HIGH RESOLUTION CADENCE (~0.5 sec cadence)
;-SPC L3:     http://sweap.cfa.harvard.edu/data/sci/sweap/spc/L3/2020/01/
;               ;e.g. spp_swp_spc_l3i_20200128_v01.cdf
;               wp_fit_uncertainty
;               vp_fit_SC
;               vp_fit_SC_uncertainty
;               vp_fit_RTN
;               vp_fit_RTN_uncertainty
;               np1_fit
;               np1_fit_uncertainty
;               wp1_fit
;               wp1_fit_uncertainty
;               vp1_fit_SC
;               vp1_fit_SC_uncertainty
;               vp1_fit_RTN
;               vp1_fit_RTN_uncertainty
;               np_moment
;               np_moment_deltahigh
;               np_moment_deltalow
;               wp_moment
;               wp_moment_deltahigh
;               wp_moment_deltalow
;               vp_moment_SC
;               vp_moment_SC_deltahigh
;               vp_moment_SC_deltalow
;               vp_moment_RTN
;               vp_moment_RTN_deltahigh
;               vp_moment_RTN_deltalow
;               na_fit
;               na_fit_uncertainty
;               wa_fit
;               wa_fit_uncertainty
;               va_fit_SC
;               va_fit_SC_uncertainty
;               va_fit_RTN
;               va_fit_RTN_uncertainty
;               n3_fit
;               n3_fit_uncertainty
;               w3_fit
;               w3_fit_uncertainty
;               v3_fit_SC
;               v3_fit_SC_uncertainty
;               v3_fit_RTN
;               v3_fit_RTN_uncertainty
;               sc_pos_HCI
;               sc_vel_HCI
;               carr_latitude
;               carr_longitude
;
;
;-SPC L2:     http://sweap.cfa.harvard.edu/data/sci/sweap/spc/L2/2020/01/
;               e.g. spp_swp_spc_l2i_20200129_v02.cdf
;               a_current
;               b_current
;               c_current
;               d_current
;               flow_angle
;               diff_charge_flux_density



;**********OLD FILES -- DON'T SEEM TO BE UPDATED ANYMOREA
;2) L2 and L3 files at http://sweap.cfa.harvard.edu/Data.html
;***NOTE: The SWEAP team hasn't updated these files in a while.

;3) L2 Energy spectra files at http://sweap.cfa.harvard.edu/data/sci/
;jwygant
;wyg@nt4data!
;SPA:
  ;Low res (22 MB) Energy flux spectra: http://sweap.cfa.harvard.edu/data/sci/sweap/spa/L2/spa_sf1/
  ;High res (500 MB) Energy flux spectra: http://sweap.cfa.harvard.edu/data/sci/sweap/spa/L2/spa_sf0/
;SPB:
  ;Low res (22 MB) Energy flux spectra: http://sweap.cfa.harvard.edu/data/sci/sweap/spb/L2/spb_sf1/
  ;High res (500 MB) Energy flux spectra: http://sweap.cfa.harvard.edu/data/sci/sweap/spb/L2/spb_sf0/




;Use READ_MYCDF to easily read multiple CDF files and specific or all variables
;within them.
;See https://spdf.gsfc.nasa.gov/CDAWlib.html for documentation of the CDAWlib
;package

;**Crib sheet designed to be run by copy/paste



rbsp_efw_init

;Must include this to use the CDFWlib data package
@compile_cdaweb


;Load L3 SWEAP files (spp_swp_spc_l3i_20181106_v08.cdf)
cdfnames = dialog_pickfile(/multiple_files,path='~/Desktop/Research/OTHER/Stuff_for_other_people/Wygant_John/SWEAP_SPC_L3_CDF/')


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
;for i=0L,n_elements(tunix)-1 do tunix[i] = time_string(strtrim(yr[i],2)+'-'+strtrim(mo[i],2)+'-'+$
;strtrim(dy[i],2)+'/'+strtrim(hr[i],2)+':'+strtrim(mn[i],2)+':'+$
;strtrim(sc[i],2)+'.'+strtrim(milli[i],2))
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

tplot_save,vars,filename='~/Desktop/psp_sweap_vars_jan2020'



end
