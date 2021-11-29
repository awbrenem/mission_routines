;+
; NAME: rbsp_efw_make_pflux
; SYNTAX:
; PURPOSE: Create the EFW Pflux CDF file
; INPUT:
; OUTPUT:
; KEYWORDS: type -> hidden - version for creating hidden file with EMFISIS data
;                -> survey - version for long-duration survey plots
;                -> if not set defaults to standard L3 version
;           script -> set if running from script. The date is read in
;           differently if so
;           version -> 1, 2, 3, etc...Defaults to 1
;
; HISTORY: Created by Aaron W Breneman, May 2014
; VERSION:
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2015-01-08 16:00:07 -0600 (Thu, 08 Jan 2015) $
;   $LastChangedRevision: 16607 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/l1_to_l2/rbsp_efw_make_l3.pro $
;-


pro rbsp_efw_make_pflux,sc,date,periodshort,periodlong,$
                        cadence_mag=cadence_mag,$
                        spinfit=spinfit,$
                        folder=folder,version=version,testing=testing,$
                        script=script,boom_pair=bp


  if ~KEYWORD_SET(bp) then bp = '12'


  title_add = strtrim(floor(float(periodshort)),2)+'-'+strtrim(floor(float(periodlong)),2)+'sec'

  probe = sc
  rbspx = 'rbsp' + sc

  ;KEEP!!!!!! Necessary when running scripts
  if keyword_set(script) then date = time_string(time_double(date),prec=-3)



  rbsp_efw_init
  skip_plot = 1                 ;set to skip restoration of cdf file and test plotting at end of program


  starttime=systime(1)
  dprint,'BEGIN TIME IS ',systime()

  if ~keyword_set(version) then version = 1
  vstr = string(version, format='(I02)')



;__________________________________________________
;Get skeleton file
;__________________________________________________

  vskeleton = '02'
  skeleton='rbsp'+sc+'_efw-pflux_00000000_v'+vskeleton+'.cdf'


  sc=strlowcase(sc)
  if sc ne 'a' and sc ne 'b' then begin
     dprint,'Invalid spacecraft: '+sc+', returning.'
     return
  endif



  if ~keyword_set(folder) then folder = '~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/'
  ;make sure we have the trailing slash on folder
  if strmid(folder,strlen(folder)-1,1) ne path_sep() then folder=folder+path_sep()
  file_mkdir,folder



  ;Use local skeleton
  if keyword_set(testing) then begin
     source_file='~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/' + skeleton
  endif else source_file='/Volumes/UserA/user_homes/rbsp_efw/Code/tdas_svn_daily/general/missions/rbsp/efw/l1_to_l2/'+skeleton



  ;make sure we have the skeleton CDF
  source_file=file_search(source_file,count=found) ; looking for single file, so count will return 0 or 1
  if ~found then begin
     dprint,'Could not find l3 v'+vskeleton+' skeleton CDF, returning.'
     return
  endif
                                ; fix single element source file array
  source_file=source_file[0]

;__________________________________________________

  store_data,tnames(),/delete

  timespan,date


  ;On RBSPa determine whether I need to use L2 spinfit files (starting 2016-06-22)
  sftype = 0
  if sc eq 'a' and time_double(date) ge time_double('2016-06-22') then sftype = 1



  rbsp_efw_poynting_flux_survey_crib,date,sc,periodshort,periodlong,$
                                    cadence_mag=cadence_mag,$
                                    spinfit=spinfit,$
                                    use_l2_spinfitdata=sftype


;stop

  rbsp_efw_poynting_flux_survey_crib,date,sc,periodshort,periodlong,$
                                    cadence_mag=cadence_mag,$
                                    spinfit=spinfit,$
                                    use_l2_spinfitdata=sftype,$
                                    /edb,/use_again ;runs the EdotB=0 version

;stop


;***NOTE: SHOULD BE ABLE TO GET THIS FROM THE L3 FILES!!!!!!!



;  rbsp_load_efw_waveform,probe=probe,datatype='vsvy',type='calibrated'
;  split_vec,rbspx+'_efw_vsvy',suffix='_'+['v1','v2','v3','v4','v5','v6']
;  get_data,rbspx+'_efw_vsvy_v1',data=v1
;  get_data,rbspx+'_efw_vsvy_v2',data=v2
;  store_data,'sum12',data={x:v1.x,y:(v1.y+v2.y)/2.}



  date = strmid(date,0,10)




;Get the official times to which all quantities are interpolated to

  get_data,'Efield_mgse',data=tmp

  ;; if keyword_set(spinfit) then get_data,'rbsp'+probe+'_efw_efield_inertial_frame_mgse',data=tmp
  ;; if ~keyword_set(spinfit) then get_data,'Esvy_mgse_vxb_removed',data=tmp
  times = tmp.x
  epoch = tplot_time_to_epoch(times,/epoch16)


;Get all the flag values

  flag_str = rbsp_efw_get_flag_values(sc,times,boom_pair=bp)


  flag_arr = flag_str.flag_arr
  bias_sweep_flag = flag_str.bias_sweep_flag
  ab_flag = flag_str.ab_flag
  charging_flag = flag_str.charging_flag


  if keyword_set(spinfit) then get_data,'rbsp'+sc+'_efw_density',data=dens
  if ~keyword_set(spinfit) then get_data,'rbsp'+sc+'_density'+bp,data=dens


;--------------------------------------------------
;Nan out various values when global flag is thrown
;--------------------------------------------------

  ;;density
  goo = where(flag_arr[*,0] eq 1)
  if goo[0] ne -1 then dens.y[goo] = -1.e31


;--------------------------------------------------
;Set a 3D flag variable for the survey plots
;--------------------------------------------------

  ;charging, autobias and eclipse flags all in one variable for convenience
  flags = [[flag_arr[*,15]],[flag_arr[*,14]],[flag_arr[*,1]]]

  ;; store_data,'*OMNI*',/delete
  store_data,'omni_imf',/delete
  store_data,'rbsp?_efw_',/delete
  store_data,'*hsk*',/delete
  store_data,'*both*',/delete
  store_data,rbspx+'_efw_vsvy*',/delete
  store_data,'B_mgse!Cmodel_subtracted',/delete

  get_data,'E-dot-B!Cflag',data=EdB_flag
  newflag = fltarr(n_elements(EdB_flag.x))
  goo = where(EdB_flag.y[*,1] eq 1)
  if goo[0] ne -1 then newflag[goo] = 1
  store_data,'E-dot-B!Cflag',data={x:EdB_flag.x,y:newflag}


  get_data,'Efield_mgse',data=dat
  times = dat.x

  tn = tnames()

  for i=0,n_elements(tn)-1 do begin
     print,'tn[i] = ',tn[i]
     tinterpol_mxn,tn[i],times,/spline
  endfor


  get_data,rbspx+'_efw_VcoroxB_mgse_interp',data=VcoroxB_mgse
  get_data,rbspx+'_efw_VscxB_mgse_interp',data=VscxB_mgse
  if keyword_set(spinfit) then get_data,rbspx+'_efw_spinaxis_gse_interp',data=spinaxis_gse
  if ~keyword_set(spinfit) then get_data,rbspx+'_spinaxis_direction_gse_interp',data=spinaxis_gse
  get_data,rbspx+'_efw_mlt_lshell_mlat_interp',data=mlt_lshell_mlat
  get_data,rbspx+'_efw_flags_charging_bias_eclipse_interp',data=flags_charging_bias_eclipse




  get_data,rbspx+'_state_pos_gse_interp',data=pos_gse
  get_data,rbspx+'_state_vel_gse_interp',data=vel_gse
  get_data,rbspx+'_state_pos_mgse_interp',data=pos_mgse
  get_data,rbspx+'_state_vel_mgse_interp',data=vel_mgse
  get_data,rbspx+'_state_radius_interp',data=radial_distance
  get_data,rbspx+'_state_mlt_interp',data=mlt
  get_data,rbspx+'_state_lshell_interp',data=lshell
  get_data,rbspx+'_state_mlat_interp',data=mlat


  get_data,rbspx+'_state_pos_gsm_interp',data=pos_gsm
  get_data,rbspx+'_mag_gsm_t01_omni_interp',data=bfield_model_gsm
  get_data,rbspx+'_mag_gse_t01_omni_interp',data=bfield_model_gse
  get_data,rbspx+'_out_iono_foot_interp',data=pos_foot_gse


;  get_data,'sum12_interp',data=sum12


  ;;This is the same as the inertial efield
  get_data,'Efield_mgse',data=efield_mgse


  get_data,rbspx+'_mag_mgse_interp',data=bfield_mgse
  get_data,rbspx+'_mag_mgse_t01_omni_interp',data=bfield_model_mgse
  get_data,rbspx+'_mag_mgse_t01_omni_dif_interp',data=bfield_minus_model_mgse
  get_data,'|B|!CnT_interp',data=bfield_magnitude

  ;; get_data,rbspx+'-Bx-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=bxmod
  ;; get_data,rbspx+'-By-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=bymod
  ;; get_data,rbspx+'-Bz-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=bzmod
  ;; bfield_minus_model_mgse = {x:bxmod.x,y:[[bxmod.y],[bymod.y],[bzmod.y]]}
  ;; get_data,'Mag_mgse_interp',data=bfield_mgse
  ;; get_data,'Mag_mgse_mod_interp',data=bfield_model_mgse
  ;; get_data,'B_mgse!Cmodel_subtracted_interp',data=bfield_minus_model_mgse



  get_data,'E-dot-B!Cflag_interp',data=EdB_flag


  get_data,'Back-ground!CB-field_interp',data=bfield_mgse_smoothed
  get_data,rbspx+'-S-para!CS!imgse!n-dot-B!imgse!n!C'+title_add+'!Cergs/cm!e2!ns_interp',data=S_para
  get_data,rbspx+'-B-para-FAC!C'+title_add+'!CnT_interp',data=dB_para


  get_data,rbspx+'-Ex-mgse!CE-dot-B!C'+title_add+'!CmV/m_interp',data=dEfield_mgsex
  get_data,rbspx+'-Ey-mgse!C'+title_add+'!CmV/m_interp',data=dEfield_mgsey
  get_data,rbspx+'-Ez-mgse!C'+title_add+'!CmV/m_interp',data=dEfield_mgsez
  dEfield_mgse = [[dEfield_mgsex.y],[dEfield_mgsey.y],[dEfield_mgsez.y]]

  get_data,rbspx+'-Bx-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=dBfield_mgsex
  get_data,rbspx+'-By-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=dBfield_mgsey
  get_data,rbspx+'-Bz-mgse!Cmodel-sub!C'+title_add+'!CnT_interp',data=dBfield_mgsez
  dBfield_mgse = [[dBfield_mgsex.y],[dBfield_mgsey.y],[dBfield_mgsez.y]]



  get_data,rbspx+'-Sx-mgse!C'+title_add+'!Cergs/cm!e2!ns_interp',data=Sx_mgse
  get_data,rbspx+'-Sy-mgse!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',data=Sy_mgse_Ex0
  get_data,rbspx+'-Sz-mgse!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',data=Sz_mgse_Ex0
  get_data,rbspx+'-Sy-mgse!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',data=Sy_mgse_EdB0
  get_data,rbspx+'-Sz-mgse!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',data=Sz_mgse_EdB0
  get_data,'Poynting-flux!C'+title_add+'!CEx-eq-0!Cfield-aligned!CEarthward=positive!Cmapped-100km!Cergs/cm!e2!ns_interp',data=S_mapped_Ex0
  get_data,'Poynting-flux!C'+title_add+'!CE-dot-B!Cfield-aligned!CEarthward=positive!Cmapped-100km!Cergs/cm!e2!ns_interp',data=S_mapped_EdB0

  get_data,'Spatial-Int.-Poynting-flux!C'+title_add+'!CEx-eq-0!Cfield-aligned!Cmapped-100km!Cergs/cm-s_interp',$
           data=S_spatial_int_mapped_Ex0
  get_data,'Time-Int.-Poynting-flux!C'+title_add+'!CEx-eq-0!Cfield-aligned!Cmapped-100km!Cergs/cm!e2!n_interp',$
           data=S_time_int_mapped_Ex0

  get_data,'Spatial-Int.-Poynting-flux!C'+title_add+'!CE-dot-B!Cfield-aligned!Cmapped-100km!Cergs/cm-s_interp',$
           data=S_spatial_int_mapped_EdB0
  get_data,'Time-Int.-Poynting-flux!C'+title_add+'!CE-dot-B!Cfield-aligned!Cmapped-100km!Cergs/cm!e2!n_interp',$
           data=S_time_int_mapped_EdB0



  get_data,'Earthward-Poynting-flux!Cfield-aligned!Cmapped=100-km!C'+title_add+'!CEx-eq-0!Cergs/cm!e2!ns_interp',data=S_earthward_logscale_FA
  get_data,'Upwards-Poynting-flux!Cfield-aligned!Cmapped=100-km!C'+title_add+'!CEx-eq-0!Cergs/cm!e2!ns_interp',data=S_upwards_logscale_FA


  get_data,'dE-azimuthal!C(eastward)!CEx-eq-0!C'+title_add+'!CmV/m_interp',data=dE_azimuthal_Ex0
  get_data,'dE-azimuthal!C(eastward)!CE-dot-B!C'+title_add+'!CmV/m_interp',data=dE_azimuthal_EdB0
  get_data,'dE-radial!C(outward)!CEx-eq-0!C'+title_add+'!CmV/m_interp',data=dE_radial_Ex0
  get_data,'dE-radial!C(outward)!CE-dot-B!C'+title_add+'!CmV/m_interp',data=dE_radial_EdB0
  get_data,'dB-azimuthal!C(eastward)!C'+title_add+'!CnT_interp',data=dB_azimuthal
  get_data,'dB-radial!C(outward)!C'+title_add+'!CnT_interp',data=dB_radial


  get_data,'S!i||!n(in-situ)!CdE-perp-2XdB-perp-1!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp2_x_dBperp1_Ex0
  get_data,'S!i||!n(in-situ)!CdE-perp-1XdB-perp-2!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp1_x_dBperp2_Ex0
  get_data,'S!i||!n(in-situ)!CdE-perp-2XdB-perp-1!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp2_x_dBperp1_EdB0
  get_data,'S!i||!n(in-situ)!CdE-perp-1XdB-perp-2!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp1_x_dBperp2_EdB0

  get_data,'S-azimuthal!C(eastward)!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',data=S_azimuthal_eastward_Ex0
  get_data,'S-azimuthal!C(eastward)!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',data=S_azimuthal_eastward_EdB0

  get_data,'S-radial!C(outward)!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',data=S_radial_outward_Ex0
  get_data,'S-radial!C(outward)!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',data=S_radial_outward_EdB0

  get_data,'S!i||!n(in-situ)!CdEperp2XdBperp1-dEperp1XdBperp2!CEx-eq-0!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0
  get_data,'S!i||!n(in-situ)!CdEperp2XdBperp1-dEperp1XdBperp2!CE-dot-B!C'+title_add+'!Cergs/cm!e2!ns_interp',$
           data=S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0

  get_data,'Vx!CExB-drift!C5min-ave!Ckm/s_interp',data=ExB_vel_xMGSE_5minavg
  get_data,'Vy!CExB-drift!C5min-ave!CEx-eq-0!Ckm/s_interp',data=ExB_vel_yMGSE_5minavg_Ex0
  get_data,'Vy!CExB-drift!C5min-ave!CE-dot-B!Ckm/s_interp',data=ExB_vel_yMGSE_5minavg_EdB0
  get_data,'Vz!CExB-drift!Csmin-ave!CEx-eq-0!Ckm/s_interp',data=ExB_vel_zMGSE_5minavg_Ex0
  get_data,'Vz!CExB-drift!Csmin-ave!CE-dot-B!Ckm/s_interp',data=ExB_vel_zMGSE_5minavg_EdB0


  get_data,rbspx+'_mag_mgse_t01_omni_dif_interp',data=tmp
  mag_diff_magnitude = sqrt(tmp.y[*,0]^2 + tmp.y[*,1]^2 + tmp.y[*,2]^2)





;--------------------------------------------------
;Remove values under 2RE for the "noperigee" fields
;--------------------------------------------------

  goore = where(radial_distance.y le 2.)


  S_para_noperigee = S_para
  dB_para_noperigee = dB_para

  dEfield_mgsex_noperigee = dEfield_mgsex
  dEfield_mgsey_noperigee = dEfield_mgsey
  dEfield_mgsez_noperigee = dEfield_mgsez
  dEfield_mgse_noperigee = dEfield_mgse

  dBfield_mgsex_noperigee = dBfield_mgsex
  dBfield_mgsey_noperigee = dBfield_mgsey
  dBfield_mgsez_noperigee = dBfield_mgsez
  dBfield_mgse_noperigee = dBfield_mgse

  Sx_mgse_noperigee = Sx_mgse
  Sy_mgse_Ex0_noperigee = Sy_mgse_Ex0
  Sz_mgse_Ex0_noperigee = Sz_mgse_Ex0
  Sy_mgse_EdB0_noperigee = Sy_mgse_EdB0
  Sz_mgse_EdB0_noperigee = Sz_mgse_EdB0
  S_mapped_Ex0_noperigee = S_mapped_Ex0
  S_mapped_EdB0_noperigee = S_mapped_EdB0

  S_earthward_logscale_FA_noperigee = S_earthward_logscale_FA
  S_upwards_logscale_FA_noperigee = S_upwards_logscale_FA
  dE_azimuthal_Ex0_noperigee = dE_azimuthal_Ex0
  dE_azimuthal_EdB0_noperigee = dE_azimuthal_EdB0
  dE_radial_Ex0_noperigee = dE_radial_Ex0
  dE_radial_EdB0_noperigee = dE_radial_EdB0
  dB_azimuthal_noperigee = dB_azimuthal
  dB_radial_noperigee = dB_radial


  S_dEperp2_x_dBperp1_Ex0_noperigee = S_dEperp2_x_dBperp1_Ex0
  S_dEperp1_x_dBperp2_Ex0_noperigee = S_dEperp1_x_dBperp2_Ex0
  S_dEperp2_x_dBperp1_EdB0_noperigee = S_dEperp2_x_dBperp1_EdB0
  S_dEperp1_x_dBperp2_EdB0_noperigee = S_dEperp1_x_dBperp2_EdB0

  S_azimuthal_eastward_Ex0_noperigee = S_azimuthal_eastward_Ex0
  S_azimuthal_eastward_EdB0_noperigee = S_azimuthal_eastward_EdB0

  S_radial_outward_Ex0_noperigee = S_radial_outward_Ex0
  S_radial_outward_EdB0_noperigee = S_radial_outward_EdB0

  S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0_noperigee = S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0
  S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0_noperigee = S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0




;---------------------

  S_para_noperigee.y[goore] = !values.f_nan
  dB_para_noperigee.y[goore] = !values.f_nan

  dEfield_mgsex_noperigee.y[goore] = !values.f_nan
  dEfield_mgsey_noperigee.y[goore] = !values.f_nan
  dEfield_mgsez_noperigee.y[goore] = !values.f_nan
  dEfield_mgse_noperigee[goore,*] = !values.f_nan

  dBfield_mgsex_noperigee.y[goore] = !values.f_nan
  dBfield_mgsey_noperigee.y[goore] = !values.f_nan
  dBfield_mgsez_noperigee.y[goore] = !values.f_nan
  dBfield_mgse_noperigee[goore,*] = !values.f_nan

  Sx_mgse_noperigee.y[goore] = !values.f_nan
  Sy_mgse_Ex0_noperigee.y[goore] = !values.f_nan
  Sz_mgse_Ex0_noperigee.y[goore] = !values.f_nan
  Sy_mgse_EdB0_noperigee.y[goore] = !values.f_nan
  Sz_mgse_EdB0_noperigee.y[goore] = !values.f_nan
  S_mapped_Ex0_noperigee.y[goore] = !values.f_nan
  S_mapped_EdB0_noperigee.y[goore] = !values.f_nan

  S_earthward_logscale_FA_noperigee.y[goore] = !values.f_nan
  S_upwards_logscale_FA_noperigee.y[goore] = !values.f_nan
  dE_azimuthal_Ex0_noperigee.y[goore] = !values.f_nan
  dE_azimuthal_EdB0_noperigee.y[goore] = !values.f_nan
  dE_radial_Ex0_noperigee.y[goore] = !values.f_nan
  dE_radial_EdB0_noperigee.y[goore] = !values.f_nan
  dB_azimuthal_noperigee.y[goore] = !values.f_nan
  dB_radial_noperigee.y[goore] = !values.f_nan

  S_dEperp2_x_dBperp1_Ex0_noperigee.y[goore] = !values.f_nan
  S_dEperp1_x_dBperp2_Ex0_noperigee.y[goore] = !values.f_nan
  S_dEperp2_x_dBperp1_EdB0_noperigee.y[goore] = !values.f_nan
  S_dEperp1_x_dBperp2_EdB0_noperigee.y[goore] = !values.f_nan

  S_azimuthal_eastward_Ex0_noperigee.y[goore] = !values.f_nan
  S_azimuthal_eastward_EdB0_noperigee.y[goore] = !values.f_nan

  S_radial_outward_Ex0_noperigee.y[goore] = !values.f_nan
  S_radial_outward_EdB0_noperigee.y[goore] = !values.f_nan

  S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0_noperigee.y[goore] = !values.f_nan
  S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0_noperigee.y[goore] = !values.f_nan






;--------------------------------------------------
                                ;These are variables for the L3 survey plots
  mlt_lshell_mlat = [[mlt.y],[lshell.y],[mlat.y]]
  ;; location = [[mlt.y],[lshell.y],[mlat.y],$
  ;;             [pos_gse.y[*,0]],[pos_gse.y[*,1]],[pos_gse.y[*,2]],$
  ;;             [vel_gse.y[*,0]],[vel_gse.y[*,1]],[vel_gse.y[*,2]],$
  ;;             [spinaxis_gse.y[*,0]],[spinaxis_gse.y[*,1]],[spinaxis_gse.y[*,2]]]
  location = [[mlt.y],[lshell.y],[mlat.y],$
              [pos_gse.y[*,0]],[pos_gse.y[*,1]],[pos_gse.y[*,2]],$
              [vel_gse.y[*,0]],[vel_gse.y[*,1]],[vel_gse.y[*,2]],$
              [spinaxis_gse.y[*,0]],[spinaxis_gse.y[*,1]],[spinaxis_gse.y[*,2]]]
  bfield_data = [[bfield_mgse.y[*,0]],[bfield_mgse.y[*,0]],[bfield_mgse.y[*,0]],$
                 [bfield_model_mgse.y[*,0]],[bfield_model_mgse.y[*,0]],[bfield_model_mgse.y[*,0]],$
                 [bfield_minus_model_mgse.y[*,0]],[bfield_minus_model_mgse.y[*,0]],[bfield_minus_model_mgse.y[*,0]],$
                 [bfield_magnitude.y],[mag_diff_magnitude]]
                                ; density_potential = [[dens.y],[sum12],[v1.y],[v2.y],[v3.y],[v4.y],[v5.y],[v6.y]]
;--------------------------------------------------

  filename = 'rbsp'+sc+'_efw-pflux_'+strjoin(strsplit(date,'-',/extract))+'_'+title_add+'_v'+vstr+'.cdf'



  file_copy,source_file,folder+filename,/overwrite

  cdfid = cdf_open(folder+filename)









  ;----------------------------------------------------------------

  ;Rename some of the variables.

  varsave1 = ['S_para',$
    'Sx_mgse',$
    'Sy_mgse_Ex0',$
    'Sz_mgse_Ex0',$
    'Sy_mgse_EdB0',$
    'Sz_mgse_EdB0',$
    'S_mapped_Ex0',$
    'S_mapped_EdB0',$
    'S_spatial_int_mapped_Ex0',$
    'S_time_int_mapped_Ex0',$
    'S_spatial_int_mapped_EdB0',$
    'S_time_int_mapped_EdB0',$
    'S_earthward_logscale_FA',$
    'S_upwards_logscale_FA',$
    'S_dEperp2_x_dBperp1_Ex0',$
    'S_dEperp1_x_dBperp2_Ex0',$
    'S_dEperp2_x_dBperp1_EdB0',$
    'S_dEperp1_x_dBperp2_EdB0',$
    'S_azimuthal_eastward_Ex0',$
    'S_azimuthal_eastward_EdB0',$
    'S_radial_outward_Ex0',$
    'S_radial_outward_EdB0',$
    'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0',$
    'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0',$
    'dB_para',$
    'dEfield_mgse',$
    'dBfield_mgse',$
    'dE_azimuthal_Ex0',$
    'dE_azimuthal_EdB0',$
    'dE_radial_Ex0',$
    'dE_radial_EdB0',$
    'dB_azimuthal',$
    'dB_radial',$
    'S_para_noperigee',$
    'dB_para_noperigee',$
    'Sx_mgse_noperigee',$
    'Sy_mgse_Ex0_noperigee',$
    'Sz_mgse_Ex0_noperigee',$
    'Sy_mgse_EdB0_noperigee',$
    'Sz_mgse_EdB0_noperigee',$
    'S_mapped_Ex0_noperigee',$
    'S_mapped_EdB0_noperigee',$
    'S_earthward_logscale_FA_noperigee',$
    'S_upwards_logscale_FA_noperigee',$
    'S_dEperp2_x_dBperp1_Ex0_noperigee',$
    'S_dEperp1_x_dBperp2_Ex0_noperigee',$
    'S_dEperp2_x_dBperp1_EdB0_noperigee',$
    'S_dEperp1_x_dBperp2_EdB0_noperigee',$
    'S_azimuthal_eastward_Ex0_noperigee',$
    'S_azimuthal_eastward_EdB0_noperigee',$
    'S_radial_outward_Ex0_noperigee',$
    'S_radial_outward_EdB0_noperigee',$
    'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0_noperigee',$
    'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0_noperigee',$
    'dEfield_mgse_noperigee',$
    'dBfield_mgse_noperigee',$
    'dE_azimuthal_Ex0_noperigee',$
    'dE_azimuthal_EdB0_noperigee',$
    'dE_radial_Ex0_noperigee',$
    'dE_radial_EdB0_noperigee',$
    'dB_azimuthal_noperigee',$
    'dB_radial_noperigee',$
    'ExB_vel_xMGSE_5minavg',$
    'ExB_vel_yMGSE_5minavg_Ex0',$
    'ExB_vel_yMGSE_5minavg_EdB0',$
    'ExB_vel_zMGSE_5minavg_Ex0',$
    'ExB_vel_zMGSE_5minavg_EdB0']

  ;Rename the CDF variables in CDF skeleton to get rid of trailing "1"
  for i=0,n_elements(varsave1)-1 do cdf_varrename,cdfid,varsave1[i]+'1',varsave1[i]




  ;Variables to save that don't have the trailing "1"

  varsave2 = ['epoch',$
    'flags_all',$
    'flags_charging_bias_eclipse',$
    'density',$
    'efield_inertial_frame_mgse',$
    'VcoroxB_mgse',$
    'VscxB_mgse',$
    'spinaxis_gse',$
    'mlt_lshell_mlat',$
    'flags_charging_bias_eclipse',$
    'pos_gsm',$
    'pos_gse',$
    'vel_gse',$
    'bfield_model_gsm',$
    'bfield_model_gse',$
    'pos_foot_gse',$
    'radial_distance',$
    'mlt_lshell_mlat',$
    'vel_mgse',$
    'pos_mgse',$
    'Bfield',$
    'EdB_flag',$
    'bfield_mgse_smoothed',$
    'bfield_mgse',$
    'bfield_model_mgse',$
    'bfield_minus_model_mgse',$
    'bfield_magnitude',$
    'bfield_magnitude_minus_modelmagnitude']



  varsave = [varsave1,varsave2]





  ;Now that we have renamed some of the variables to our liking,
  ;get list of all the variable names in the CDF file.
  inq = cdf_inquire(cdfid)
  CDFvarnames = ''
  for varNum = 0, inq.nzvars-1 do begin $
    stmp = cdf_varinq(cdfid,varnum,/zvariable) & $
    if stmp.recvar eq 'VARY' then CDFvarnames = [CDFvarnames,stmp.name]
  endfor
  CDFvarnames = CDFvarnames[1:n_elements(CDFvarnames)-1]



  ;Delete all variables we don't want to save.
  for qq=0,n_elements(CDFvarnames)-1 do begin $
    tstt = array_contains(varsave,CDFvarnames[qq]) & $
    if not tstt then print,'Deleting var:  ', CDFvarnames[qq]
    if not tstt then cdf_vardelete,cdfid,CDFvarnames[qq]
  endfor





;----------------------------------------------------------------





  cdf_varput,cdfid,'epoch',epoch
  cdf_varput,cdfid,'flags_all',transpose(flag_arr)
  cdf_varput,cdfid,'flags_charging_bias_eclipse',transpose(flags)


;--------------------------------------------------
;Populate CDF file for L3 version
;--------------------------------------------------

  cdf_varput,cdfid,'density',transpose(dens.y)
  cdf_varput,cdfid,'efield_inertial_frame_mgse',transpose(efield_mgse.y)

  if keyword_set(spinfit) then cdf_varput,cdfid,'VcoroxB_mgse',transpose(VcoroxB_mgse.y) else cdf_vardelete,cdfid,'VcoroxB_mgse'
  if keyword_set(spinfit) then cdf_varput,cdfid,'VscxB_mgse',transpose(VscxB_mgse.y) else cdf_vardelete,cdfid,'VscxB_mgse'
  cdf_varput,cdfid,'spinaxis_gse',transpose(spinaxis_gse.y)
  ;cdf_varput,cdfid,'Vavg',transpose(sum12.y)
  cdf_varput,cdfid,'mlt_lshell_mlat',transpose(mlt_lshell_mlat)
  if keyword_set(spinfit) then cdf_varput,cdfid,'flags_charging_bias_eclipse',transpose(flags_charging_bias_eclipse.y) else $
     cdf_vardelete,cdfid,'flags_charging_bias_eclipse'

  cdf_varput,cdfid,'pos_gsm',transpose(pos_gsm.y)
  cdf_varput,cdfid,'pos_gse',transpose(pos_gse.y)
  cdf_varput,cdfid,'vel_gse',transpose(vel_gse.y)
  cdf_varput,cdfid,'bfield_model_gsm',transpose(bfield_model_gsm.y)
  cdf_varput,cdfid,'bfield_model_gse',transpose(bfield_model_gse.y)
  cdf_varput,cdfid,'pos_foot_gse',transpose(pos_foot_gse.y)

  cdf_varput,cdfid,'radial_distance',transpose(radial_distance.y)
  cdf_varput,cdfid,'mlt_lshell_mlat',transpose(mlt_lshell_mlat)

  cdf_varput,cdfid,'vel_mgse',transpose(vel_mgse.y)
  cdf_varput,cdfid,'pos_mgse',transpose(pos_mgse.y)

  cdf_varput,cdfid,'Bfield',transpose(bfield_data)

  cdf_varput,cdfid,'EdB_flag',EdB_flag.y
  cdf_varput,cdfid,'bfield_mgse_smoothed',transpose(bfield_mgse_smoothed.y)



  cdf_varput,cdfid,'S_para',transpose(S_para.y)
  cdf_varput,cdfid,'Sx_mgse',transpose(Sx_mgse.y)
  cdf_varput,cdfid,'Sy_mgse_Ex0',transpose(Sy_mgse_Ex0.y)
  cdf_varput,cdfid,'Sz_mgse_Ex0',transpose(Sz_mgse_Ex0.y)
  cdf_varput,cdfid,'Sy_mgse_EdB0',transpose(Sy_mgse_EdB0.y)
  cdf_varput,cdfid,'Sz_mgse_EdB0',transpose(Sz_mgse_EdB0.y)
  cdf_varput,cdfid,'S_mapped_Ex0',transpose(S_mapped_Ex0.y)
  cdf_varput,cdfid,'S_mapped_EdB0',transpose(S_mapped_EdB0.y)
  cdf_varput,cdfid,'S_spatial_int_mapped_Ex0',transpose(S_spatial_int_mapped_Ex0.y)
  cdf_varput,cdfid,'S_time_int_mapped_Ex0',transpose(S_time_int_mapped_Ex0.y)
  cdf_varput,cdfid,'S_spatial_int_mapped_EdB0',transpose(S_spatial_int_mapped_EdB0.y)
  cdf_varput,cdfid,'S_time_int_mapped_EdB0',transpose(S_time_int_mapped_EdB0.y)
  cdf_varput,cdfid,'S_earthward_logscale_FA',transpose(S_earthward_logscale_FA.y)
  cdf_varput,cdfid,'S_upwards_logscale_FA',transpose(S_upwards_logscale_FA.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_Ex0',transpose(S_dEperp2_x_dBperp1_Ex0.y)
  cdf_varput,cdfid,'S_dEperp1_x_dBperp2_Ex0',transpose(S_dEperp1_x_dBperp2_Ex0.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_EdB0',transpose(S_dEperp2_x_dBperp1_EdB0.y)
  cdf_varput,cdfid,'S_dEperp1_x_dBperp2_EdB0',transpose(S_dEperp1_x_dBperp2_EdB0.y)
  cdf_varput,cdfid,'S_azimuthal_eastward_Ex0',transpose(S_azimuthal_eastward_Ex0.y)
  cdf_varput,cdfid,'S_azimuthal_eastward_EdB0',transpose(S_azimuthal_eastward_EdB0.y)
  cdf_varput,cdfid,'S_radial_outward_Ex0',transpose(S_radial_outward_Ex0.y)
  cdf_varput,cdfid,'S_radial_outward_EdB0',transpose(S_radial_outward_EdB0.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0',$
             transpose(S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0',$
             transpose(S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0.y)

  cdf_varput,cdfid,'dB_para',transpose(dB_para.y)
  cdf_varput,cdfid,'dEfield_mgse',transpose(dEfield_mgse)
  cdf_varput,cdfid,'dBfield_mgse',transpose(dBfield_mgse)
  cdf_varput,cdfid,'dE_azimuthal_Ex0',transpose(dE_azimuthal_Ex0.y)
  cdf_varput,cdfid,'dE_azimuthal_EdB0',transpose(dE_azimuthal_EdB0.y)
  cdf_varput,cdfid,'dE_radial_Ex0',transpose(dE_radial_Ex0.y)
  cdf_varput,cdfid,'dE_radial_EdB0',transpose(dE_radial_EdB0.y)
  cdf_varput,cdfid,'dB_azimuthal',transpose(dB_azimuthal.y)
  cdf_varput,cdfid,'dB_radial',transpose(dB_radial.y)

  cdf_varput,cdfid,'S_para_noperigee',transpose(S_para_noperigee.y)
  cdf_varput,cdfid,'dB_para_noperigee',transpose(dB_para_noperigee.y)
  cdf_varput,cdfid,'Sx_mgse_noperigee',transpose(Sx_mgse_noperigee.y)
  cdf_varput,cdfid,'Sy_mgse_Ex0_noperigee',transpose(Sy_mgse_Ex0_noperigee.y)
  cdf_varput,cdfid,'Sz_mgse_Ex0_noperigee',transpose(Sz_mgse_Ex0_noperigee.y)
  cdf_varput,cdfid,'Sy_mgse_EdB0_noperigee',transpose(Sy_mgse_EdB0_noperigee.y)
  cdf_varput,cdfid,'Sz_mgse_EdB0_noperigee',transpose(Sz_mgse_EdB0_noperigee.y)
  cdf_varput,cdfid,'S_mapped_Ex0_noperigee',transpose(S_mapped_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_mapped_EdB0_noperigee',transpose(S_mapped_EdB0_noperigee.y)

  cdf_varput,cdfid,'S_earthward_logscale_FA_noperigee',transpose(S_earthward_logscale_FA_noperigee.y)
  cdf_varput,cdfid,'S_upwards_logscale_FA_noperigee',transpose(S_upwards_logscale_FA_noperigee.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_Ex0_noperigee',transpose(S_dEperp2_x_dBperp1_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_dEperp1_x_dBperp2_Ex0_noperigee',transpose(S_dEperp1_x_dBperp2_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_EdB0_noperigee',transpose(S_dEperp2_x_dBperp1_EdB0_noperigee.y)
  cdf_varput,cdfid,'S_dEperp1_x_dBperp2_EdB0_noperigee',transpose(S_dEperp1_x_dBperp2_EdB0_noperigee.y)
  cdf_varput,cdfid,'S_azimuthal_eastward_Ex0_noperigee',transpose(S_azimuthal_eastward_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_azimuthal_eastward_EdB0_noperigee',transpose(S_azimuthal_eastward_EdB0_noperigee.y)
  cdf_varput,cdfid,'S_radial_outward_Ex0_noperigee',transpose(S_radial_outward_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_radial_outward_EdB0_noperigee',transpose(S_radial_outward_EdB0_noperigee.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0_noperigee',$
             transpose(S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_Ex0_noperigee.y)
  cdf_varput,cdfid,'S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0_noperigee',$
             transpose(S_dEperp2_x_dBperp1_minus_dEperp1_x_dBperp2_EdB0_noperigee.y)

  cdf_varput,cdfid,'dEfield_mgse_noperigee',transpose(dEfield_mgse_noperigee)
  cdf_varput,cdfid,'dBfield_mgse_noperigee',transpose(dBfield_mgse_noperigee)
  cdf_varput,cdfid,'dE_azimuthal_Ex0_noperigee',transpose(dE_azimuthal_Ex0_noperigee.y)
  cdf_varput,cdfid,'dE_azimuthal_EdB0_noperigee',transpose(dE_azimuthal_EdB0_noperigee.y)
  cdf_varput,cdfid,'dE_radial_Ex0_noperigee',transpose(dE_radial_Ex0_noperigee.y)
  cdf_varput,cdfid,'dE_radial_EdB0_noperigee',transpose(dE_radial_EdB0_noperigee.y)
  cdf_varput,cdfid,'dB_azimuthal_noperigee',transpose(dB_azimuthal_noperigee.y)
  cdf_varput,cdfid,'dB_radial_noperigee',transpose(dB_radial_noperigee.y)

  cdf_varput,cdfid,'ExB_vel_xMGSE_5minavg',transpose(ExB_vel_xMGSE_5minavg.y)
  cdf_varput,cdfid,'ExB_vel_yMGSE_5minavg_Ex0',transpose(ExB_vel_yMGSE_5minavg_Ex0.y)
  cdf_varput,cdfid,'ExB_vel_yMGSE_5minavg_EdB0',transpose(ExB_vel_yMGSE_5minavg_EdB0.y)
  cdf_varput,cdfid,'ExB_vel_zMGSE_5minavg_Ex0',transpose(ExB_vel_zMGSE_5minavg_Ex0.y)
  cdf_varput,cdfid,'ExB_vel_zMGSE_5minavg_EdB0',transpose(ExB_vel_zMGSE_5minavg_EdB0.y)

  cdf_varput,cdfid,'bfield_mgse',transpose(bfield_mgse.y)
  cdf_varput,cdfid,'bfield_model_mgse',transpose(bfield_model_mgse.y)
  cdf_varput,cdfid,'bfield_minus_model_mgse',transpose(bfield_minus_model_mgse.y)
  cdf_varput,cdfid,'bfield_magnitude',transpose(bfield_magnitude.y)
  cdf_varput,cdfid,'bfield_magnitude_minus_modelmagnitude',transpose(mag_diff_magnitude)





  cdf_close, cdfid
  store_data,tnames(),/delete


end
