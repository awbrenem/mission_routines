;+
; NAME: rbsp_efw_make_l3
; SYNTAX:
; PURPOSE: Create the EFW L3 CDF file
; INPUT:
; OUTPUT:
; KEYWORDS: type -> hidden - version for creating hidden file with EMFISIS data
;                -> survey - version for long-duration survey plots
;                -> if not set defaults to standard L3 version
;           script -> set if running from script. The date is read in
;           differently if so
;           version -> 1, 2, 3, etc...Defaults to 1
;           boom_pair -> defaults to '12' for spinfit data. Can change
;           to '34'
;
; HISTORY: Created by Aaron W Breneman, May 2014
; VERSION:
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2016-09-02 08:56:49 -0500 (Fri, 02 Sep 2016) $
;   $LastChangedRevision: 21782 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/l1_to_l2/rbsp_efw_make_l3.pro $
;-



;*****************
;TODO LIST
;*****************
;NECESSARY TO LOAD HSK DATA???
;WHAT TO RENAME "DIAG" VARIABLES? (include spin axis component somehow?)
;Remove data during eclipses?

;*****************

pro rbsp_efw_make_l3_OBSOLETE,sc,date,$
  folder=folder,$
  version=version,$
  type=type,$
  testing=testing,$
  script=script,$
  boom_pair=bp,$
  density_min=density_min

  print,date
  timespan,date


  ;KEEP!!!!!! Necessary when running scripts
  if keyword_set(script) then date = time_string(double(date),prec=-3)
  if ~keyword_set(testing) then begin
     openw,lun,'output.txt',/get_lun
     printf,lun,'date = ',date
     printf,lun,'date type: ',typename(date)
     printf,lun,'probe = ',sc
     printf,lun,'probe type: ',typename(sc)
     printf,lun,'bp = ',bp
     printf,lun,'bp type: ',typename(bp)
     close,lun
     free_lun,lun
  endif


  rbx = 'rbsp'+sc


  ;First (and only) load of these
;  rbsp_load_spice_kernels
  rbsp_efw_init


  ;Define keyword inheritance to pass to subroutines. This will ensure that
  ;subsequent routines don't reload spice kernels or rbsp_efw_init
  extra = create_struct('no_spice_load',1,$
                        'no_rbsp_efw_init',1)




  if ~keyword_set(type) then type = 'L3'
  if ~keyword_set(bp) then bp = '12'


  skip_plot = 1 ;set to skip restoration of cdf file and test plotting at end of program


  starttime=systime(1)
  dprint,'BEGIN TIME IS ',systime()

  if ~keyword_set(version) then version = 3
  vstr = string(version, format='(I02)')


;----------------------------------------------------------------
;Get skeleton file
;----------------------------------------------------------------

;Skeleton file
  vskeleton = 'XX'
  skeleton=rbx+'_efw-lX_00000000_v'+vskeleton+'.cdf'


  sc=strlowcase(sc)
  if sc ne 'a' and sc ne 'b' then begin
     dprint,'Invalid spacecraft: '+sc+', returning.'
     return
  endif






  ; Use local skeleton
  if keyword_set(testing) then begin
     folder = '~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/'
     source_file=folder + skeleton
;     if ~keyword_set(folder) then folder = '~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/'
     ; make sure we have the trailing slash on folder
     if strmid(folder,strlen(folder)-1,1) ne path_sep() then folder=folder+path_sep()
     file_mkdir,folder
  endif else source_file='/Volumes/UserA/user_homes/rbsp_efw/Code/tdas_svn_daily/general/missions/rbsp/efw/l1_to_l2/'+skeleton



  ; make sure we have the skeleton CDF
  source_file=file_search(source_file,count=found) ; looking for single file, so count will return 0 or 1
  if ~found then begin
     dprint,'Could not find l3 v'+vskeleton+' skeleton CDF, returning.'
     return
  endif
  ; fix single element source file array
  source_file=source_file[0]

;------------------------------------------------

  ;clean slate
  store_data,tnames(),/delete




  ;Load both the spinfit data and also the E*B=0 version
  rbsp_efw_edotb_to_zero_crib,date,sc,$
    /noplot,$
    suffix='edotb',$
    boom_pair=bp,$
    _extra=extra;,/noremove



  ;Get By/Bx and Bz/Bx from E*B=0 calculation
  get_data,'B2Bx_ratio',data=b2bx_ratio
  badyx = where(b2bx_ratio.y[*,0] gt 3.732)
  badzx = where(b2bx_ratio.y[*,1] gt 3.732)


  ;Get spinaxis component
  get_data,rbx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',data=diagEx
  diagEx = diagEx.y[*,0]


  ;Have two versions. First has all E*B=0 data, second has E*B=0 bad data removed
  diagEx1 = diagEx
  diagEx2 = diagEx
  if badyx[0] ne -1 then diagEx2[badyx,0] = !values.f_nan
  if badzx[0] ne -1 then diagEx2[badzx,0] = !values.f_nan


  ;Get the official times to which all quantities are interpolated to
  get_data,rbx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',data=tmp
  times = tmp.x
  epoch = tplot_time_to_epoch(times,/epoch16)


;Get all the flag values
  flag_str = rbsp_efw_get_flag_values(sc,times,boom_pair=bp,density_min=density_min)

  flag_arr = flag_str.flag_arr
  bias_sweep_flag = flag_str.bias_sweep_flag
  ab_flag = flag_str.ab_flag
  charging_flag = flag_str.charging_flag
  charging_flag_extreme = flag_str.charging_flag_extreme


  copy_data,rbx+'_density'+bp,rbx+'_density'
  get_data,rbx+'_density',data=dens





;--------------------------------------------------
;save all spinfit resolution Efield quantities
;--------------------------------------------------

  get_data,rbx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',data=tmp
  spinfit_vxb_edotb = tmp.y
  spinfit_vxb = tmp.y & spinfit_vxb[*,0] = -1.0E31


  get_data,rbx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',data=tmp
  spinfit_vxb_coro_edotb = tmp.y
  spinfit_vxb_coro = spinfit_vxb_coro_edotb & spinfit_vxb_coro[*,0] = -1.0E31






;--------------------------------------
;SUBTRACT OFF MODEL FIELD
;--------------------------------------


  if type eq 'hidden' then begin

    model = 't89'
    rbsp_efw_dcfield_removal_crib,sc,$
      /noplot,$
      model=model,$
      _extra=extra

    tinterpol_mxn,rbx+'_mag_mgse',times,/overwrite,/spline
    tinterpol_mxn,rbx+'_mag_mgse_model',times,/overwrite,/spline
    tinterpol_mxn,rbx+'_mag_mgse_t89_dif',times,/overwrite,/spline

    get_data,rbx+'_mag_mgse',data=mag_mgse
    get_data,rbx+'_mag_mgse_'+model+'',data=mag_model
    get_data,rbx+'_mag_mgse_t89_dif',data=mag_diff
    mag_model_magnitude = sqrt(mag_model.y[*,0]^2 + mag_model.y[*,1]^2 + mag_model.y[*,2]^2)
    mag_data_magnitude = sqrt(mag_mgse.y[*,0]^2 + mag_mgse.y[*,1]^2 + mag_mgse.y[*,2]^2)
    mag_diff_magnitude = mag_data_magnitude - mag_model_magnitude

    bfield_data = [[mag_mgse.y[*,0]],[mag_mgse.y[*,0]],[mag_mgse.y[*,0]],$
                   [mag_model.y[*,0]],[mag_model.y[*,0]],[mag_model.y[*,0]],$
                   [mag_diff.y[*,0]],[mag_diff.y[*,0]],[mag_diff.y[*,0]],$
                   [mag_data_magnitude],[mag_diff_magnitude]]

  endif


;--------------------------------------------------
;Nan out various values when global flag is thrown
;--------------------------------------------------

  ;density
  goo = where(flag_arr[*,0] eq 1)
  if goo[0] ne -1 then dens.y[goo] = -1.e31


;--------------------------------------------------
;Set a 3D flag variable for the survey plots
;--------------------------------------------------

  ;charging, extreme charging, autobias and eclipse flags all in one variable for convenience
  flags = [[flag_arr[*,15]],[flag_arr[*,16]],[flag_arr[*,14]],[flag_arr[*,1]]]




  ;Downsample variables to cadence of spinfit data
  vartmp = rbx+'_'+['E_coro_mgse','vscxb','state_vel_coro_mgse','state_pos_gse',$
  'state_vel_gse','state_mlt','state_mlat',$
  'state_lshell','ME_orbitnumber','spinaxis_direction_gse']
  for i=0,n_elements(vartmp)-1 do tinterpol_mxn,vartmp[i],times,/overwrite,/spline
  for i=1,6 do tinterpol_mxn,rbx +'_efw_vsvy_'+strtrim(i,2),times,/overwrite,/spline








  get_data,rbx+'_vxb',data=vxb
  get_data,rbx+'_state_pos_gse',data=pos_gse
  get_data,rbx+'_state_vel_gse',data=vel_gse
  get_data,rbx+'_E_coro_mgse',data=ecoro_mgse
  get_data,rbx+'_state_vel_coro_mgse',data=vcoro_mgse
  get_data,rbx+'_state_mlt',data=mlt
  get_data,rbx+'_state_mlat',data=mlat
  get_data,rbx+'_state_lshell',data=lshell
  get_data,rbx+'_spinaxis_direction_gse',data=sa
  get_data,rbx+'_angles',data=angles

  get_data,rbx +'_efw_vsvy_1',data=v1
  get_data,rbx +'_efw_vsvy_2',data=v2
  get_data,rbx +'_efw_vsvy_3',data=v3
  get_data,rbx +'_efw_vsvy_4',data=v4
  get_data,rbx +'_efw_vsvy_5',data=v5
  get_data,rbx +'_efw_vsvy_6',data=v6

  if bp eq '12' then vavg = (v1.y + v2.y)/2.
  if bp eq '34' then vavg = (v3.y + v4.y)/2.
  if bp eq '56' then vavg = (v5.y + v6.y)/2.
  if bp eq '24' then vavg = (v2.y + v4.y)/2.
  if bp eq '14' then vavg = (v1.y + v4.y)/2.
  if bp eq '13' then vavg = (v1.y + v3.y)/2.
  if bp eq '23' then vavg = (v2.y + v3.y)/2.



;--------------------------------------------------
  ;These are variables for the L3 survey plots
  mlt_lshell_mlat = [[mlt.y],[lshell.y],[mlat.y]]
  location = [[mlt.y],[lshell.y],[mlat.y],$
              [pos_gse.y[*,0]],[pos_gse.y[*,1]],[pos_gse.y[*,2]],$
              [vel_gse.y[*,0]],[vel_gse.y[*,1]],[vel_gse.y[*,2]],$
              [sa.y[*,0]],[sa.y[*,1]],[sa.y[*,2]]]
  density_potential = [[dens.y],[vavg],[v1.y],[v2.y],[v3.y],[v4.y],[v5.y],[v6.y]]

;--------------------------------------------------

  if type eq 'L3'     then filename = rbx+'_efw-l3_'+strjoin(strsplit(date,'-',/extract))+'_v'+vstr+'.cdf'
  if type eq 'hidden' then filename = rbx+'_efw-l3_'+strjoin(strsplit(date,'-',/extract))+'_v'+vstr+'_hidden.cdf'
  if type eq 'survey' then filename = rbx+'_efw-l3_'+strjoin(strsplit(date,'-',/extract))+'_v'+vstr+'_survey.cdf'



  ;Rename the CDF file with the date
  file_copy,source_file,folder+filename,/overwrite

  cdfid = cdf_open(folder+filename)




  ;Get list of all the variable names in the CDF file.
  inq = cdf_inquire(cdfid)
  CDFvarnames = ''
  for varNum = 0, inq.nzvars-1 do begin $
    stmp = CDF_VARINQ(cdfid,varnum,/ZVARIABLE ) & $
    if stmp.recvar eq 'VARY' then CDFvarnames = [CDFvarnames,stmp.name]
  endfor
  CDFvarnames = CDFvarnames[1:n_elements(CDFvarnames)-1]






  cdf_varput,cdfid,'epoch',epoch
  cdf_varput,cdfid,'flags_all',transpose(flag_arr)
  cdf_varput,cdfid,'flags_charging_bias_eclipse',transpose(flags)



;;--------------------------------------------------
;;Remove values during eclipse times
;;--------------------------------------------------

  goo = where(flags[*,2] eq 1)

  if goo[0] ne -1 then begin
     spinfit_vxb[goo,*] = !values.f_nan
     spinfit_vxb_coro[goo,*] = !values.f_nan
     dens.y[goo] = !values.f_nan
     vavg[goo] = !values.f_nan
  endif


;--------------------------------------------------
;Populate CDF file for L3 version
;--------------------------------------------------



  ;Fix labeling so that V? is changed to the proper antenna
  b1 = strmid(bp,0,1) & b2 = strmid(bp,1,1)


;  varnum = CDF_VARNUM(cdfid, 'efield_inertial_frame_mgse')
;  text = 'Inertial frame spinfit electric field in MGSE coord system. Derived from the V'+b1+' and V'+b2+' spinplane booms. The Vsc x B field is subtracted off, where Vsc is the spacecraft velocity and B is the measured ambient magnetic field'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','efield_inertial_frame_mgse', text
;
;  varnum = CDF_VARNUM(cdfid, 'efield_inertial_frame_mgse_edotb_zero')
;  text = 'Spin-fit electric field calculation (from the V'+b1+' and V'+b2+' spinplane booms) in the MGSE coordinate system. The Vsc x B field is subtracted off, where Vsc is the spacecraft velocity and B is the measured ambient magnetic field. The spin axis Efield is calculated with the E*B=0 assumption'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','efield_inertial_frame_mgse_edotb_zero', text
;
;  varnum = CDF_VARNUM(cdfid, 'efield_corotation_frame_mgse')
;  text = 'Corotation frame spinfit electric field in MGSE coord system. Derived from the V'+b1+' and V'+b2+' spinplane booms. The corotation and Vsc x B electric fields are subtracted off, where Vsc is the spacecraft velocity and B is the measured ambient magnetic field'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','efield_corotation_frame_mgse', text
;
;  varnum = CDF_VARNUM(cdfid, 'efield_corotation_frame_mgse_edotb_zero')
;  text = 'Spin-fit electric field calculation (from the V'+b1+' and V'+b2+' spinplane booms) in the MGSE coordinate system. The corotation and Vsc x B electric fields are subtracted off, where Vsc is the spacecraft velocity and B is the measured ambient magnetic field. The spin axis Efield is calculated with the E*B=0 assumption'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','efield_corotation_frame_mgse_edotb_zero', text
;
;  varnum = CDF_VARNUM(cdfid, 'density')
;  text = 'Density estimation from sc potential, determined from the form density = A1*exp(b*x)+A2*exp(c*x), where x=(V'+b1+' + V'+b2+')/2, where the constants are found by comparison to upper hybrid line from EMFISIS. The fit is updated approximately monthly. For a more accurate density (or verification) for a specific time period contact the EFW team. Values below 10cm-3 and above 3000 cm-3 are untrustworthy and are removed'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','density', text
;
;  varnum = CDF_VARNUM(cdfid, 'Vavg')
;  text = 'Average of opposing boom potentials (V'+b1+'+V'+b2+')/2. This is a proxy for plasma density. Positive excursions can indicate negative spacecraft charging.'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','Vavg', text
;
;  varnum = CDF_VARNUM(cdfid, 'density_potential')
;  text = 'Density calculated from antenna potential (calibrated to upper hybrid line)!C(V'+b1+'+V'+b2+')/2 - proxy for density!CV1....V6 - antenna potentials'
;  CDF_ATTPUT, cdfid, 'VAR_NOTES','density_potential', text



;; Walk through all of the ZVariable attributes
;inq = CDF_INQUIRE(cdfid)
; FOR attrNum = 0, inq.natts-1 DO BEGIN $
;    CDF_ATTGET_ENTRY, cdfid, attrNum, varNum, attType, value, status, /ZVARIABLE, CDF_TYPE=cdfType, ATTRIBUTE_NAME=attName & $
;    PRINT, "attr_name = ", attName  , ", ", value



  if type eq 'L3' then begin

     cdf_varput,cdfid,'efield_inertial_frame_mgse',transpose(spinfit_vxb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse',transpose(spinfit_vxb_coro)
     cdf_varput,cdfid,'VcoroxB_mgse',transpose(ecoro_mgse.y)
     cdf_varput,cdfid,'VscxB_mgse',transpose(vxb.y)
     cdf_varput,cdfid,'density',dens.y
     cdf_varput,cdfid,'Vavg',vavg
     cdf_varput,cdfid,'mlt_lshell_mlat',transpose(mlt_lshell_mlat)
     cdf_varput,cdfid,'position_gse',transpose(pos_gse.y)
     cdf_varput,cdfid,'velocity_gse',transpose(vel_gse.y)
     cdf_varput,cdfid,'spinaxis_gse',transpose(sa.y)
     cdf_varput,cdfid,'diagEx1',diagEx1
     cdf_varput,cdfid,'diagEx2',diagEx2
     cdf_varput,cdfid,'diagBratio',transpose(b2bx_ratio.y)

     cdf_varput,cdfid,'efield_inertial_frame_mgse_edotb_zero',transpose(spinfit_vxb_edotb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse_edotb_zero',transpose(spinfit_vxb_coro_edotb)
;     cdf_vardelete,cdfid,'efield_inertial_frame_mgse_edotb_zero'
;     cdf_vardelete,cdfid,'efield_corotation_frame_mgse_edotb_zero'

     cdf_vardelete,cdfid,'bfield_mgse'
     cdf_vardelete,cdfid,'bfield_model_mgse'
     cdf_vardelete,cdfid,'bfield_minus_model_mgse'
     cdf_vardelete,cdfid,'bfield_magnitude'
     cdf_vardelete,cdfid,'bfield_magnitude_minus_modelmagnitude'
     cdf_vardelete,cdfid,'Bfield'
     cdf_vardelete,cdfid,'density_potential'
     cdf_vardelete,cdfid,'ephemeris'
     cdf_vardelete,cdfid,'angle_spinplane_Bo'
     cdf_vardelete,cdfid,'bias_current'

  endif



;--------------------------------------------------
;Populate CDF file for survey version
;--------------------------------------------------


  if type eq 'survey' then begin

     cdf_varput,cdfid,'efield_inertial_frame_mgse',transpose(spinfit_vxb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse',transpose(spinfit_vxb_coro)
     cdf_varput,cdfid,'efield_inertial_frame_mgse_edotb_zero',transpose(spinfit_vxb_edotb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse_edotb_zero',transpose(spinfit_vxb_coro_edotb)
     cdf_varput,cdfid,'VcoroxB_mgse',transpose(ecoro_mgse.y)
     cdf_varput,cdfid,'VscxB_mgse',transpose(vxb.y)
     cdf_varput,cdfid,'Bfield',transpose(bfield_data)
     cdf_varput,cdfid,'density_potential',transpose(density_potential)
     cdf_varput,cdfid,'ephemeris',transpose(location)
     cdf_varput,cdfid,'angle_spinplane_Bo',transpose(angles.y)
     cdf_varput,cdfid,'diagEx1',diagEx1
     cdf_varput,cdfid,'diagEx2',diagEx2
     cdf_varput,cdfid,'diagBratio',transpose(b2bx_ratio.y)

     cdf_vardelete,cdfid,'density'
     cdf_vardelete,cdfid,'Vavg'
     cdf_vardelete,cdfid,'pos_gse'
     cdf_vardelete,cdfid,'vel_gse'
     cdf_vardelete,cdfid,'spinaxis_gse'
     cdf_vardelete,cdfid,'mlt_lshell_mlat'
     cdf_vardelete,cdfid,'bfield_mgse'
     cdf_vardelete,cdfid,'bfield_model_mgse'
     cdf_vardelete,cdfid,'bfield_minus_model_mgse'
     cdf_vardelete,cdfid,'bfield_magnitude_minus_modelmagnitude'
     cdf_vardelete,cdfid,'bfield_magnitude'

  endif


;--------------------------------------------------
;Populate CDF file for hidden version
;--------------------------------------------------

  if type eq 'hidden' then begin

     cdf_varput,cdfid,'efield_inertial_frame_mgse',transpose(spinfit_vxb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse',transpose(spinfit_vxb_coro)
     cdf_varput,cdfid,'efield_inertial_frame_mgse_edotb_zero',transpose(spinfit_vxb_edotb)
     cdf_varput,cdfid,'efield_corotation_frame_mgse_edotb_zero',transpose(spinfit_vxb_coro_edotb)
     cdf_varput,cdfid,'VcoroxB_mgse',transpose(ecoro_mgse.y)
     cdf_varput,cdfid,'VscxB_mgse',transpose(vxb.y)
     cdf_varput,cdfid,'bfield_magnitude',mag_data_magnitude
     cdf_varput,cdfid,'bfield_mgse',transpose(mag_mgse.y)
     cdf_varput,cdfid,'bfield_model_mgse',transpose(mag_model.y)
     cdf_varput,cdfid,'bfield_minus_model_mgse',transpose(mag_diff.y)
     cdf_varput,cdfid,'bfield_magnitude_minus_modelmagnitude',mag_diff_magnitude
     cdf_varput,cdfid,'density',dens.y
     cdf_varput,cdfid,'Vavg',vavg
     cdf_varput,cdfid,'mlt_lshell_mlat',transpose(mlt_lshell_mlat)
     cdf_varput,cdfid,'pos_gse',transpose(pos_gse.y)
     cdf_varput,cdfid,'vel_gse',transpose(vel_gse.y)
     cdf_varput,cdfid,'spinaxis_gse',transpose(sa.y)
     cdf_varput,cdfid,'angle_spinplane_Bo',transpose(angles.y)
     cdf_varput,cdfid,'diagEx1',diagEx1
     cdf_varput,cdfid,'diagEx2',diagEx2
     cdf_varput,cdfid,'diagBratio',transpose(b2bx_ratio.y)


     cdf_vardelete,cdfid,'Bfield'
     cdf_vardelete,cdfid,'ephemeris'
     cdf_vardelete,cdfid,'density_potential'

  endif


  cdf_close, cdfid
  store_data,tnames(),/delete


end
