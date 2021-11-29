;+
;NAME: rbsp_efw_vxb_subtract_crib.pro
;PURPOSE:
; Creates vxB subtracted tplot variables for the 32 S/s EFW efield data
; using all four spin plane probes 
; ****NOTE: IF YOU DON'T NEED THE FULL 32 S/s RESOLUTION THEN CONSIDER USING
; ****rbsp_efw_spinfit_vxb_subtract_crib.pro
; ****THE DATA WILL BE AT SPINPERIOD CADENCE BUT WILL BE CLEANER
;
;CALLING SEQUENCE:
; timespan,'2014-01-01'
; rbsp_efw_vxb_subtract_crib,'a'
;
;INPUT: 'a' or 'b' for probe
;KEYWORDS:
; ql -> use EMFISIS quicklook UVW data instead of 4sec GSE L3. This will be
;       despun and spinfit. This has the advantage of not having to wait for
;       the EMFISIS L3 data set to be produced.
;       Note: Defaults to 4-sec EMFISIS data in GSE. This is actually superior
;       to the hires GSE data b/c it doesn't have spurious data spikes.
;       It works about the same as JBT's method of despinning and then
;       spinfitting the EMFISIS quicklook UVW data.
;
; hires -> use EMFISIS hires L3 GSE data instead of 4sec GSE L3.
; qa -> load the QA testing L1 EFW data instead of the usual L1 data
; _extra --> useful keywords are
;     no_spice_load
;     no_rbsp_efw_init
;
;OUTPUT:
;HISTORY:
; Created by Aaron Breneman, UMN, Dec 2012
;	email: awbrenem@gmail.com
;REQUIRES: THEMIS TDAS software
;			 http://themis.ssl.berkeley.edu/software.shtml
;		   as well as SPICE software
;
;$LastChangedBy: nikos $
;$LastChangedDate: 2020-05-21 22:36:46 -0500 (Thu, 21 May 2020) $
;$LastChangedRevision: 28720 $
;$URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/examples/rbsp_efw_vxb_subtract_crib.pro $
;-


pro rbsp_efw_vxb_subtract_crib,probe,$
  noplot=noplot,$
  ql=ql,$
  l2=l2,$
  hires=hires,$
  qa=qa,$
  bad_probe=bad_probe,$
  _extra=extra


  if ~KEYWORD_SET(ql) and ~KEYWORD_SET(l2) then level = 'l3'
  if KEYWORD_SET(ql) then level = 'ql'
  if KEYWORD_SET(l2) then level = 'l2'

  ;Set timerange if it's not already set
  x = timerange()
  date = strmid(time_string(x[0]),0,10)


  ;initialize RBSP environment
  rbsp_efw_init,_extra=extra

  ;set desired probe
  rbspx = 'rbsp'+probe

  ;Set other quantities
  suffix = ''



  ;Load definitive sc positions and velocities   
  rbsp_load_spice_cdf_file,probe


 
  rbsp_load_efw_esvy_mgse,probe=probe,bad_probe=bad_probe,_extra=extra
 



  ;Load EMFISIS data
  if keyword_set(hires) and ~keyword_set(l2) then rbsp_load_emfisis,probe=probe,coord='gse',cadence='hires',level='l3'
  if ~keyword_set(hires) and ~keyword_set(ql) and ~keyword_set(l2) then rbsp_load_emfisis,probe=probe,coord='gse',cadence='1sec',level='l3'
 ; if keyword_set(l2) then   rbsp_load_emfisis,probe=probe,coord='uvw',level='l2'
;  if keyword_set(ql) then rbsp_load_emfisis,probe=probe,/quicklook


  ;Check for data existence
  if keyword_set(hires) and ~keyword_set(l2) then get_data,rbspx+'_emfisis_l3_hires_gse_Mag',data=dd2
  if ~keyword_set(hires) and ~keyword_set(ql) and ~keyword_set(l2) then get_data,rbspx+'_emfisis_l3_1sec_gse_Mag',data=dd2
  ;if keyword_set(ql) then get_data,rbspx+'_emfisis_quicklook_Mag',data=dd2
  ;if keyword_set(l2) then get_data,rbspx+'_emfisis_l2_uvw_Mag',data=dd2



  if ~is_struct(dd2) then begin
     print,'******NO MAG DATA TO LOAD.....rbsp_efw_DCfield_removal_crib.pro*******'
     return
  endif



  ;Transform the Mag data to MGSE coordinates
  if ~keyword_set(ql) and ~keyword_set(l2) then begin

     if keyword_set(hires) then $
      get_data,rbspx+'_emfisis_l3_hires_gse_Mag',data=tmpp else $
      get_data,rbspx+'_emfisis_l3_1sec_gse_Mag',data=tmpp

     get_data,rbspx+'_spinaxis_direction_gse',data=wsc_GSE

     wsc_GSE_tmp = [[interpol(wsc_GSE.y[*,0],wsc_GSE.x,tmpp.x)],$
                    [interpol(wsc_GSE.y[*,1],wsc_GSE.x,tmpp.x)],$
                    [interpol(wsc_GSE.y[*,2],wsc_GSE.x,tmpp.x)]]


     if keyword_set(hires) then $
      rbsp_gse2mgse,rbspx+'_emfisis_l3_hires_gse_Mag',reform(wsc_GSE_tmp),$
      newname=rbspx+'_emfisis_l3_hires_mgse_Mag' else $
      rbsp_gse2mgse,rbspx+'_emfisis_l3_1sec_gse_Mag',reform(wsc_GSE_tmp),$
      newname=rbspx+'_emfisis_l3_1sec_mgse_Mag'

     if keyword_set(hires) then $
      copy_data,rbspx+'_emfisis_l3_hires_mgse_Mag',rbspx+'_mag_mgse' else $
      copy_data,rbspx+'_emfisis_l3_1sec_mgse_Mag',rbspx+'_mag_mgse'

  endif


;  if keyword_set(ql) then begin
;
 ;   ;Create the dlimits structure for the EMFISIS quantity. The spinfit
;    ;program needs to see that the coords are 'uvw'
;    data_att = {coord_sys:'uvw'}
;    dlim = {data_att:data_att}
;    store_data,rbspx +'_emfisis_quicklook_Mag',data=dd2,dlimits=dlim
;
;    ;spinfit the mag data and transform to MGSE
;    rbsp_decimate,rbspx +'_emfisis_quicklook_Mag', upper = 2
;    rbsp_spinfit,rbspx +'_emfisis_quicklook_Mag', plane_dim = 0
;    rbsp_cotrans,rbspx +'_emfisis_quicklook_Mag_spinfit',rbspx+'_mag_mgse', /dsc2mgse
;  endif

;  if keyword_set(l2) then begin
;
 ;   ;Create the dlimits structure for the EMFISIS quantity. Jianbao's spinfit program needs
  ;  ;to see that the coords are 'uvw'
   ; data_att = {coord_sys:'uvw'}
  ;  dlim = {data_att:data_att}
  ;  store_data,rbspx +'_emfisis_l2_uvw_Mag',data=dd2,dlimits=dlim
;
;
;    ;spinfit the mag data and transform to MGSE
;    rbsp_decimate,rbspx +'_emfisis_l2_uvw_Mag', upper = 2
;    rbsp_spinfit,rbspx +'_emfisis_l2_uvw_Mag', plane_dim = 0
;    rbsp_cotrans,rbspx +'_emfisis_l2_uvw_Mag_spinfit',rbspx+'_mag_mgse', /dsc2mgse
;
;  endif





  ;Determine corotation Efield
  rbsp_corotation_efield,probe,date,level=level,_extra=extra
  ;Get Vsc x B (motional) and Vcoro x B (corotation) electric fields
  get_data,rbspx+'_efw_esvy_mgse',data=dtmp     ;....times to interpolate to 
  rbsp_efw_vxb_create,rbspx+'_state_vel_mgse',rbspx+'_mag_mgse',dtmp.x,title = '(Vsc x B)!CmV/m!CMGSE'
  copy_data,'vxb',rbspx+'_vscxb_mgse'
  copy_data,rbspx+'_E_coro_mgse',rbspx+'_vcoroxb_mgse'







  ;Find total residual Efield. We need to subtract this off so that we can apply the effective 
  ;antenna length to the Efield measured only by the probes, and not any motional or corotation 
  ;field. This field consists of the Vsc x B field and the 
  ;Vcoro x B field (NOTE: Vcoro is the minus Vcoro field) 
  add_data,rbspx+'_vscxb_mgse',rbspx+'_vcoroxb_mgse',newname='Einertial+coro_mgse'




  ;create the rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed' variable
  ;Subtract off Ecoro + Emotional
  dif_data,rbspx+'_efw_esvy_mgse','Einertial+coro_mgse',newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed'


 



  ;Apply crude antenna effective length correction to minimize residual field. 
  ;Do this to the 
  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed', data = d
  if is_struct(d) then begin
     d.y[*, 1] *= 0.947d        ;found by S. Thaller
     d.y[*, 2] *= 0.947d
     store_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed', data = d
  endif





  ;Now that effective antenna length has been applied, add back in corotation field 
  ;to get the inertial frame Efield. NOTE: Vcoro has a negative sign, so we'll have to subtract here.
  dif_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed',rbspx+'_vcoroxb_mgse',newname=rbspx+'_efw_esvy_mgse_vxb_removed'





  options,rbspx+'_efw_esvy_mgse_vxb_removed','colors',[2,4,6]
  options,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed','colors',[2,4,6]
  options,rbspx+'_efw_esvy_mgse_vxb_removed','ysubtitle',''
  options,rbspx+'_mag_mgse','ytitle','Bfield MGSE!C[nT]'
  options,rbspx+'_mag_mgse','ysubtitle',''
  options,rbspx+'_efw_esvy_mgse','ztitle','RBSP'+probe+'!CEFW'
  options,rbspx+'_efw_esvy_mgse','ysubtitle',''
  options,rbspx+'_efw_esvy_mgse','colors',[0,1,2]




  ;Delete unnecessary variables
  store_data,['bfield_data_gei',$
              rbspx+'_state_vel_coro_gei',rbspx+'_E_coro_mgse',rbspx+'_E_coro_gse',rbspx+'_E_coro_gei',$
              rbspx+'_emfisis_l3_hires_gse_delta',rbspx+'_emfisis_l3_hires_gse_lambda',$
              rbspx+'_emfisis_l3_4sec_gse_rms',$
              rbspx+'_emfisis_l3_1sec_gse_rms',$
              rbspx+'_state_pos_gei',rbspx+'_efw_esvy_ccsds_data_BEB_config',$
              rbspx+'_efw_esvy_ccsds_data_DFB_config',$
              'vxb',$
              rbspx+'_spinaxis_direction_gse_interp',$
              rbspx+'_ew',rbspx+'_eu_fixed',rbspx+'_ev_fixed'],/delete




end
