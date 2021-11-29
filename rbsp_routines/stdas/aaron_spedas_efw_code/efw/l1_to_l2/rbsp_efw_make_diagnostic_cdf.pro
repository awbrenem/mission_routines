;+
; NAME:
;   rbsp_efw_make_diagnostic_cdf
;
; PURPOSE:
;   Generate diagnostic EFW CDF files 
;
;
; Two methods are used to calculate E23 and E24. (See rbsp_efw_edotb_to_zero_crib_v2.pro)
;
; CALLING SEQUENCE:
;   rbsp_efw_make_diagnostic_cdf, sc, date
;
; ARGUMENTS:
;   sc: IN, REQUIRED
;         'a' or 'b'
;   date: IN, REQUIRED
;         A date string in format like '2013-02-13'
;
; KEYWORDS:
;   folder: IN, OPTIONAL
;         Default is something like
;           !rbsp_efw.local_data_dir/rbspa/l2/spinfit/2012/
;
;   boom_pair -> specify for the spinfit routine. Defaults to '12' but
;   can be set to '34'
;
; OLD CDF files that are now obsolete are:
;      rbspa_efw-l2_e-spinfit-mgse_20130103_v01.cdf
;      rbspa_efw-l2_esvy_despun_20130105_v01.cdf
;      rbspa_efw-l2_vsvy-hires_20130105_v01.cdf
;      rbspa_efw-l2_combo_20130101_v03_hr.cdf
;      rbspa_efw-l2_combo_20130101_v03.cdf
;      rbspa_efw-l2_combo_pfaff_00000000_v01.cdf
;      rbspa_efw-l2_combo_wygant_00000000_v01.cdf
;
;
; HISTORY:
;   2014-12-02: Created by Aaron W Breneman, U. Minnesota
;				
;
; VERSION:
; $LastChangedBy: aaronbreneman $
; $LastChangedDate: 2016-02-04 11:23:28 -0600 (Thu, 04 Feb 2016) $
; $LastChangedRevision: 19899 $
; $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/l1_to_l2/rbsp_efw_make_l2.pro $
;
;-

pro rbsp_efw_make_diagnostic_cdf,sc,date,$
                                 folder=folder,$
                                 magExtra = magExtra,$
                                 version = version,$
                                 save_flags = save_flags,$
                                 no_spice_load = no_spice_load,$
                                 no_cdf = no_cdf,$
                                 testing=testing,$
                                 hires=hires,$
                                 ql=ql,$
                                 density_min=dmin

  

  if ~keyword_set(dmin) then dmin = 10.
  if ~keyword_set(ql) then ql = 0
  if ~keyword_set(testing) then begin
     openw,lun,'output.txt',/get_lun
     printf,lun,'date = ',date
     printf,lun,'probe = ',sc
;  printf,lun,'hires = ',hires
     close,lun
     free_lun,lun
  endif

  compile_opt idl2
  rbsp_efw_init

  if n_elements(version) eq 0 then version = 2
  vstr = string(version, format='(I02)')


  rbspx='rbsp' + strlowcase(sc[0])
  rbx = rbspx + '_'
  probe = sc

;------------ Set up paths. BEGIN. ----------------------------

  if ~keyword_set(no_cdf) then begin

     year = strmid(date, 0, 4)

     if ~keyword_set(folder) then folder = !rbsp_efw.local_data_dir + $
                                           'rbsp' + strlowcase(sc[0]) + path_sep() + $
                                           'l2' + path_sep() + $
                                           'spinfit' + path_sep() + $
                                           year + path_sep()

                                ; make sure we have the trailing slash on folder
     if strmid(folder,strlen(folder)-1,1) ne path_sep() then folder=folder+path_sep()
     if ~keyword_set(no_cdf) then file_mkdir, folder


     ;; Grab the skeleton file.
     skeleton='/Volumes/UserA/user_homes/kersten/Code/tdas_svn_daily/general/missions/rbsp/efw/l1_to_l2/'+$
              rbspx+'_efw-diagnostic_00000000_v01.cdf'

     found = 1
                                ; make sure we have the skeleton CDF
     if ~keyword_set(testing) then skeletonFile=file_search(skeleton,count=found)
     if keyword_set(testing) then $
        skeletonfile = '~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/rbsp'+$
                       sc+'_efw-diagnostic_00000000_v01.cdf'


     if ~found then begin
        dprint,'Could not find skeleton CDF, returning.'
        return
     endif
                                ; fix single element source file array
     skeletonFile=skeletonFile[0]

  endif

  if keyword_set(testing) then folder = '~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/l1_to_l2/'



;------------ Set up paths. END. ----------------------------


  skip = 'no'

  if skip eq 'no' then begin

     store_data,tnames(),/delete
     timespan,date
     rbsp_load_spice_kernels

     rbspx = 'rbsp'+probe
     rbsp_load_efw_waveform, probe=sc, datatype =['vsvy','esvy'], coord = 'uvw',/noclean


     get_data,'rbsp'+sc+'_efw_vsvy',data=vsvy
     epoch_v = tplot_time_to_epoch(vsvy.x,/epoch16)
     times_v = vsvy.x

     ;; full resolution (V1+V2)/2
     vsvy_vavg = [[(vsvy.y[*,0] - vsvy.y[*,1])/2.],$
                  [(vsvy.y[*,2] - vsvy.y[*,3])/2.],$
                  [(vsvy.y[*,4] - vsvy.y[*,5])/2.]]
     
     split_vec, 'rbsp'+sc+'_efw_vsvy', suffix='_V'+['1','2','3','4','5','6']
     get_data,'rbsp'+sc+'_efw_vsvy',data=vsvy

     
     get_data,'rbsp'+sc+'_efw_esvy',data=esvy
     epoch_e = tplot_time_to_epoch(esvy.x,/epoch16)
     times_e = esvy.x

     tinterpol_mxn,'rbsp'+sc+'_efw_esvy','rbsp'+sc+'_efw_vsvy',newname='rbsp'+sc+'_efw_esvy',/spline
     get_data,'rbsp'+sc+'_efw_esvy',data=esvy_v





;;-------------------------------------------------------------
;;Call the edotb crib for various boom pairs, as well as method 1 and
;;2 for boom pairs that are 90 deg apart (e.g. V2-V3)
;;-------------------------------------------------------------


     bp = ['12','34','23','23','24','24']
     method = ['1','1','1','2','1','2']
     for u=0,n_elements(bp)-1 do begin
        if u eq 0 then rerun = 0 else rerun = 1
        rbsp_efw_edotb_to_zero_crib_v2,$
           date,sc,/no_spice_load,/noplot,suffix='edotb',$
           boom_pair=bp[u],ql=ql,datatype='vsvy',method=method[u],$
           rerun=rerun

        copy_data,rbspx+'_efw_esvy_spinfit',$
                  'tmp'+method[u]+'_sf_'+bp[u]
        copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
                  'tmp'+method[u]+'_sf_vxb_'+bp[u]
        copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
                  'tmp'+method[u]+'_sf_vxb_coro_'+bp[u]
        copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
                  'tmp'+method[u]+'_sf_vxb_edotb_'+bp[u]
        copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
                  'tmp'+method[u]+'_sf_vxb_coro_edotb_'+bp[u]

        store_data,[rbspx+'_efw_esvy_spinfit',$
                    rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
                    rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
                    rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
                    rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete

     endfor

     

     ;;Rename 'tmp' variables
     for u=0,n_elements(bp)-1 do begin 
   
        if method[u] eq 1 then begin
           copy_data,'tmp1_sf_'+bp[u],$
                     rbspx+'_efw_esvy_spinfit_'+bp[u]
           copy_data,'tmp1_sf_vxb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_'+bp[u]
           copy_data,'tmp1_sf_vxb_coro_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_'+bp[u]
           copy_data,'tmp1_sf_vxb_edotb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_'+bp[u]
           copy_data,'tmp1_sf_vxb_coro_edotb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_'+bp[u]

        endif
        
        if method[u] eq 2 then begin
           copy_data,'tmp2_sf_'+bp[u],$
                     rbspx+'_efw_esvy_spinfit_'+bp[u]+'_m2'
           copy_data,'tmp2_sf_vxb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_'+bp[u]+'_m2'
           copy_data,'tmp2_sf_vxb_coro_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_'+bp[u]+'_m2'
           copy_data,'tmp2_sf_vxb_edotb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_'+bp[u]+'_m2'
           copy_data,'tmp2_sf_vxb_coro_edotb_'+bp[u],$
                     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_'+bp[u]+'_m2'
        endif

     endfor


     store_data,['*tmp1_sf*','*tmp2_sf*'],/delete



;Get the official times to which all quantities are interpolated to
     get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12',data=tmp
     times = tmp.x
     epoch = tplot_time_to_epoch(times,/epoch16)



;--------------------------------------------------
;Get flag values (also gets density values from v12 and v34)
;--------------------------------------------------


     ;; ;;Get hires density values
     ;; if type eq 'combo' then begin
     goo_str = rbsp_efw_get_flag_values(sc,times_v,density_min=dmin)
     copy_data,rbspx+'_density12',rbspx+'_density12_hires'
     copy_data,rbspx+'_density34',rbspx+'_density34_hires'
     store_data,[rbspx+'_density12',rbspx+'_density34'],/delete
     ;; endif


     flag_str = rbsp_efw_get_flag_values(sc,times,density_min=dmin)

     flag_arr = flag_str.flag_arr
     bias_sweep_flag = flag_str.bias_sweep_flag
     ab_flag = flag_str.ab_flag
     charging_flag = flag_str.charging_flag


     

;--------------------------------------------------
;save all spinfit resolution Efield quantities
;--------------------------------------------------

     for u=0,n_elements(bp)-1 do begin 


        tmp = 0.
        tinterpol_mxn,rbspx+'_efw_esvy_spinfit_'+bp[u],times,$
                      newname=rbspx+'_efw_esvy_spinfit_'+bp[u],/spline
        get_data,     rbspx+'_efw_esvy_spinfit_'+bp[u],data=tmp
        if is_struct(tmp) then begin
           tmp.y[*,0] = -1.0E31
           spinfit_esvy_12 = tmp.y
           tmp = 0.
        endif
        
                                ;Spinfit with corotation field
        tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12',times,$
                      newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12',/spline
        get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12',data=tmp
        if is_struct(tmp) then begin
           tmp.y[*,0] = -1.0E31
           spinfit_vxb_12 = tmp.y
           tmp = 0.
        endif
                                ;Spinfit with corotation field and E*B=0
        tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_12',times,$
                      newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_12',/spline
        get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_12',data=tmp
        if is_struct(tmp) then begin
           spinfit_vxb_edotb_12 = tmp.y
        endif

                                ;Spinfit without corotation field
        tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_12',times,$
                      newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_12',/spline
        get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_12',data=tmp
        if is_struct(tmp) then begin
           tmp.y[*,0] = -1.0E31
           spinfit_vxb_coro_12 = tmp.y
           tmp = 0.
        endif
                                ;Spinfit without corotation field and E*B=0
        tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_12',times,$
                      newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_12',/spline
        get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_12',data=tmp
        if is_struct(tmp) then begin
           spinfit_vxb_coro_edotb_12 = tmp.y
        endif


     endfor


                                ;----

     tmp = 0.
     tinterpol_mxn,rbspx+'_efw_esvy_spinfit_34',times,$
                   newname=rbspx+'_efw_esvy_spinfit_34',/spline
     get_data,     rbspx+'_efw_esvy_spinfit_34',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_esvy_34 = tmp.y
        tmp = 0.
     endif
                                ;Spinfit with corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_34',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_34',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_34',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_34 = tmp.y
        tmp = 0.
     endif
                                ;Spinfit with corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_34',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_34',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_34',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_edotb_34 = tmp.y
     endif
                                ;Spinfit without corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_34',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_34',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_34',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_coro_34 = tmp.y
        tmp = 0.
     endif
                                ;Spinfit without corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_34',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_34',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_34',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_coro_edotb_34 = tmp.y
     endif
;------





     tmp = 0.
     tinterpol_mxn,rbspx+'_efw_esvy_spinfit_23',times,$
                   newname=rbspx+'_efw_esvy_spinfit_23',/spline
     get_data,     rbspx+'_efw_esvy_spinfit_23',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_esvy_23 = tmp.y
        tmp = 0.
     endif
     tmp = 0.
     tinterpol_mxn,rbspx+'_efw_esvy_spinfit_23_m2',times,$
                   newname=rbspx+'_efw_esvy_spinfit_23_m2',/spline
     get_data,     rbspx+'_efw_esvy_spinfit_23_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_esvy_23_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit with corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_23 = tmp.y
        tmp = 0.
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_23_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit with corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_edotb_23 = tmp.y
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23_m2',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_edotb_23_m2 = tmp.y
     endif

                                ;Spinfit without corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_coro_23 = tmp.y
        tmp = 0.
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_coro_23_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit without corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_coro_edotb_23 = tmp.y
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23_m2',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_coro_edotb_23_m2 = tmp.y
     endif







;--------



     tmp = 0.
     tinterpol_mxn,rbspx+'_efw_esvy_spinfit_24',times,$
                   newname=rbspx+'_efw_esvy_spinfit_24',/spline
     get_data,     rbspx+'_efw_esvy_spinfit_24',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_esvy_24 = tmp.y
        tmp = 0.
     endif
     tmp = 0.
     tinterpol_mxn,rbspx+'_efw_esvy_spinfit_24_m2',times,$
                   newname=rbspx+'_efw_esvy_spinfit_24_m2',/spline
     get_data,     rbspx+'_efw_esvy_spinfit_24_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_esvy_24_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit with corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_24 = tmp.y
        tmp = 0.
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_24_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit with corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_edotb_24 = tmp.y
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_edotb_24_m2 = tmp.y
     endif

                                ;Spinfit without corotation field
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_coro_24 = tmp.y
        tmp = 0.
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24_m2',data=tmp
     if is_struct(tmp) then begin
        tmp.y[*,0] = -1.0E31
        spinfit_vxb_coro_24_m2 = tmp.y
        tmp = 0.
     endif

                                ;Spinfit without corotation field and E*B=0
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_coro_edotb_24 = tmp.y
     endif
     tinterpol_mxn,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24_m2',times,$
                   newname=rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24_m2',/spline
     get_data,     rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24_m2',data=tmp
     if is_struct(tmp) then begin
        spinfit_vxb_coro_edotb_24_m2 = tmp.y
     endif



;--------------------------------------
;SUBTRACT OFF MODEL FIELD
;--------------------------------------

     model = 't89'
     rbsp_efw_DCfield_removal_crib,sc,/no_spice_load,/noplot,model=model
     

;--------------------------------------
; SC potentials (V1+V2)/2 and (V3+V4)/2
;--------------------------------------

     get_data,rbspx +'_efw_vsvy_V1',data=v1
     get_data,rbspx +'_efw_vsvy_V2',data=v2
     get_data,rbspx +'_efw_vsvy_V3',data=v3
     get_data,rbspx +'_efw_vsvy_V4',data=v4
     get_data,rbspx +'_efw_vsvy_V5',data=v5
     get_data,rbspx +'_efw_vsvy_V6',data=v6

     sum12 = (v1.y + v2.y)/2.	
     sum34 = (v3.y + v4.y)/2.	
     sum56 = (v5.y + v6.y)/2.	

     sum56[*] = -1.0E31
     
     store_data,'sum12',data={x:v1.x,y:sum12}
     tinterpol_mxn,'sum12',times,newname='sum12',/spline
     get_data,'sum12',data=sum12
     sum12=sum12.y

     store_data,'sum34',data={x:v3.x,y:sum34}
     tinterpol_mxn,'sum34',times,newname='sum34',/spline
     get_data,'sum34',data=sum34
     sum34=sum34.y



                                ;Interpolate single-ended measurements
                                ;to low cadence for combo file
     tinterpol_mxn,rbspx+'_efw_vsvy',times,newname=rbspx+'_efw_vsvy_combo',/spline
     get_data,rbspx+'_efw_vsvy_combo',data=vsvy_spinres
     

;--------------------------------------------------
;Nan out various values when global flag is thrown
;--------------------------------------------------

     ;;density
     tinterpol_mxn,rbspx+'_density12',times,newname=rbspx+'_density12',/spline
     get_data,rbspx+'_density12',data=dens12
     goo = where(flag_arr[*,0] eq 1)
     if goo[0] ne -1 and is_struct(dens12) then dens12.y[goo] = -1.e31

     tinterpol_mxn,rbspx+'_density34',times,newname=rbspx+'_density34',/spline
     get_data,rbspx+'_density34',data=dens34
     goo = where(flag_arr[*,0] eq 1)
     if goo[0] ne -1 and is_struct(dens34) then dens34.y[goo] = -1.e31


     tinterpol_mxn,rbspx+'_density12_hires',times_e,newname=rbspx+'_density12_hires',/spline
     get_data,rbspx+'_density12_hires',data=dens12_hires

     tinterpol_mxn,rbspx+'_density34_hires',times_e,newname=rbspx+'_density34_hires',/spline
     get_data,rbspx+'_density34_hires',data=dens34_hires


;--------------------------------------------------
;Set a 3D flag variable for the survey plots
;--------------------------------------------------

     ;;charging, autobias and eclipse flags all in one variable for convenience
     flags = [[flag_arr[*,15]],[flag_arr[*,14]],[flag_arr[*,1]]]

     ;; if type ne 'spinfit_both_boompairs' then begin
     ;;    ;;Set the density flag based on the antenna pair
     ;;    ;;used. We don't want to do this if type =
     ;;    ;;'spinfit_both_boompairs' because I include density values
     ;;    ;;obtained from both V12 and V34 in these CDF files
     ;;    flag_arr[*,16] = 0
     ;;    if bp eq '12' then begin
     ;;       goo = where(dens12.y eq -1.e31)
     ;;       if goo[0] ne -1 then flag_arr[goo,16] = 1
     ;;    endif else begin
     ;;       goo = where(dens34.y eq -1.e31)
     ;;       if goo[0] ne -1 then flag_arr[goo,16] = 1
     ;;    endelse
     ;; endif     


     
;the times for the mag spinfit can be slightly different than the times for the
;Esvy spinfit. 
     tinterpol_mxn,rbspx+'_mag_mgse',times,newname=rbspx+'_mag_mgse',/spline
     get_data,rbspx+'_mag_mgse',data=mag_mgse


;Downsample the GSE position and velocity variables to cadence of spinfit data
     tinterpol_mxn,rbspx+'_E_coro_mgse',times,newname=rbspx+'_E_coro_mgse',/spline
     tinterpol_mxn,rbspx+'_vscxb',times,newname='vxb',/spline
     tinterpol_mxn,rbspx+'_state_vel_coro_mgse',times,newname=rbspx+'_state_vel_coro_mgse',/spline
     tinterpol_mxn,rbspx+'_state_pos_gse',times,newname=rbspx+'_state_pos_gse',/spline
     tinterpol_mxn,rbspx+'_state_vel_gse',times,newname=rbspx+'_state_vel_gse',/spline
     get_data,'vxb',data=vxb
     get_data,rbspx+'_state_pos_gse',data=pos_gse
     get_data,rbspx+'_state_vel_gse',data=vel_gse
     get_data,rbspx+'_E_coro_mgse',data=ecoro_mgse
     get_data,rbspx+'_state_vel_coro_mgse',data=vcoro_mgse
     
     tinterpol_mxn,rbspx+'_mag_mgse_'+model,times,newname=rbspx+'_mag_mgse_'+model,/spline
     tinterpol_mxn,rbspx+'_mag_mgse_t89_dif',times,newname=rbspx+'_mag_mgse_t89_dif',/spline
     get_data,rbspx+'_mag_mgse_'+model,data=mag_model
     get_data,rbspx+'_mag_mgse_t89_dif',data=mag_diff

     mag_model_magnitude = sqrt(mag_model.y[*,0]^2 + mag_model.y[*,1]^2 + mag_model.y[*,2]^2)
     mag_data_magnitude = sqrt(mag_mgse.y[*,0]^2 + mag_mgse.y[*,1]^2 + mag_mgse.y[*,2]^2)
     mag_diff_magnitude = mag_data_magnitude - mag_model_magnitude

     tinterpol_mxn,rbspx+'_state_mlt',times,newname=rbspx+'_state_mlt',/spline
     tinterpol_mxn,rbspx+'_state_mlat',times,newname=rbspx+'_state_mlat',/spline
     tinterpol_mxn,rbspx+'_state_lshell',times,newname=rbspx+'_state_lshell',/spline

     get_data,rbspx+'_state_mlt',data=mlt
     get_data,rbspx+'_state_mlat',data=mlat
     get_data,rbspx+'_state_lshell',data=lshell

     tinterpol_mxn,rbspx+'_spinaxis_direction_gse',times,newname=rbspx+'_spinaxis_direction_gse',/spline
     get_data,rbspx+'_spinaxis_direction_gse',data=sa

     tinterpol_mxn,'angles',times,newname='angles',/spline
     get_data,'angles',data=angles


  endif                         ;for skipping processing






;; ;;------
;; ;;RESTORE CURRENT PROGRESS FOR TESTING
;; ;Save all current files in session to a single file called barrel.tplot 
;; fileroot = '~/Desktop/'
;; tplot_save,'*',filename='rbsp_efw_make_diagnostic_cdf'    ;don't add .tplot 
;; tplot_restore,filenames=fileroot+'rbsp_efw_make_diagnostic_cdf.tplot'

;; rbspx = 'rbspa'
;; sc ='a'
;; ;;-----

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12',etimes,emgse

  ;;Rotate to field-aligned coord
  tinterpol_mxn,rbspx+'_spinaxis_direction_gse',etimes,/spline
  get_data,rbspx+'_spinaxis_direction_gse_interp',data=wsc_gse

  rbsp_gse2mgse,rbspx+'_state_pos_gse',reform(wsc_gse.y[0,*]),$
                newname=rbspx+'_state_pos_mgse'

  get_data,rbspx+'_state_pos_mgse',data=mgse_pos
  radial_vec = [[mgse_pos.y[*,0]],[mgse_pos.y[*,1]],[mgse_pos.y[*,2]]]


  tinterpol_mxn,rbspx+'_mag_mgse_for_subtract',etimes,/spline
  
  rbsp_detrend,rbspx+'_mag_mgse_for_subtract_interp',60.*20.
  tplot,[rbspx+'_mag_mgse_for_subtract_interp',rbspx+'_mag_mgse_for_subtract_interp_smoothed']
  
  get_data,rbspx+'_mag_mgse_for_subtract_interp_smoothed',btimes,bmgse
  get_data,rbspx+'_mag_mgse_for_subtract_interp',btimes,bwmgse
  ntimes = n_elements(etimes)
  perp1 = 'azimuthal!C(eastward)'
  perp2 = 'radial!C(outward)'
  

;;-----------------------------------



;define orthogonal perpendicular unit vectors based on smoothed
;background field
;;Radial (outwards), azimuthal (eastwards), Bo coord system
  perp1_dir = fltarr(ntimes,3)
  for xx=0L,ntimes-1 do perp1_dir[xx,*] = crossp(bmgse[xx,*],radial_vec[xx,*])
  perp2_dir = fltarr(ntimes,3)
  for xx=0L,ntimes-1 do perp2_dir[xx,*] = crossp(perp1_dir[xx,*],bmgse[xx,*])

  perp1_mag = sqrt(perp1_dir[*,0]^2 + perp1_dir[*,1]^2 + perp1_dir[*,2]^2)
  perp2_mag = sqrt(perp2_dir[*,0]^2 + perp2_dir[*,1]^2 + perp2_dir[*,2]^2)
  bmag = sqrt(bmgse[*,0]^2 + bmgse[*,1]^2 + bmgse[*,2]^2)

  perp1_dir = [[perp1_dir[*,0]/perp1_mag],[perp1_dir[*,1]/perp1_mag],[perp1_dir[*,2]/perp1_mag]]
  perp2_dir = [[perp2_dir[*,0]/perp2_mag],[perp2_dir[*,1]/perp2_mag],[perp2_dir[*,2]/perp2_mag]]
  par_dir = [[bmgse[*,0]/bmag],[bmgse[*,1]/bmag],[bmgse[*,2]/bmag]]

  ;; ang_test = fltarr(ntimes)
  ;; for i=0L,ntimes-1 do ang_test[i] = acos(total(perp1_dir[i,*]*par_dir[i,*]))/!dtor




  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_12',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]  
  Erab12 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab12',data = {x:etimes,y:Erab12},dlim={constant:[0],colors:[0]}


  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_12',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]  
  Erab12c = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab12c',data = {x:etimes,y:Erab12c},dlim={constant:[0],colors:[0]}


  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab23 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab23',data = {x:etimes,y:Erab23},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23_m2',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab23_m2 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab23_m2',data = {x:etimes,y:Erab23_m2},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab23c = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab23c',data = {x:etimes,y:Erab23c},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23_m2',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab23c_m2 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab23c_m2',data = {x:etimes,y:Erab23c_m2},dlim={constant:[0],colors:[0]}


  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab24 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab24',data = {x:etimes,y:Erab24},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab24_m2 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab24_m2',data = {x:etimes,y:Erab24_m2},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab24_m2 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab24_m2',data = {x:etimes,y:Erab24_m2},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab24c = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab24c',data = {x:etimes,y:Erab24c},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24_m2',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab24c_m2 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab24c_m2',data = {x:etimes,y:Erab24c_m2},dlim={constant:[0],colors:[0]}


  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_34',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab34 = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab34',data = {x:etimes,y:Erab34},dlim={constant:[0],colors:[0]}

  get_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_34',etimes,emgse
  E_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp1[xx] = emgse[xx,0]*perp1_dir[xx,0] + emgse[xx,1]*perp1_dir[xx,1] + emgse[xx,2]*perp1_dir[xx,2]
  E_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_perp2[xx] = emgse[xx,0]*perp2_dir[xx,0] + emgse[xx,1]*perp2_dir[xx,1] + emgse[xx,2]*perp2_dir[xx,2]
  E_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do E_par[xx] = emgse[xx,0]*par_dir[xx,0] + emgse[xx,1]*par_dir[xx,1] + emgse[xx,2]*par_dir[xx,2]
  Erab34c = [[E_perp2],[E_perp1],[E_par]]
  store_data,'Erab34c',data = {x:etimes,y:Erab34c},dlim={constant:[0],colors:[0]}
  ;; ylim,['Erab12'],-50,50
  ;; options,'Erab12','colors',[0,50,250]
  ;; tplot,['Erab12']



  
  B_perp1  = fltarr(ntimes)
  for xx=0L,ntimes-1 do B_perp1[xx] = bwmgse[xx,0]*perp1_dir[xx,0] + bwmgse[xx,1]*perp1_dir[xx,1] + bwmgse[xx,2]*perp1_dir[xx,2]
  B_perp2  = fltarr(ntimes)
  for xx=0L,ntimes-1 do B_perp2[xx] = bwmgse[xx,0]*perp2_dir[xx,0] + bwmgse[xx,1]*perp2_dir[xx,1] + bwmgse[xx,2]*perp2_dir[xx,2]
  B_par  = fltarr(ntimes)
  for xx=0L,ntimes-1 do B_par[xx] = bwmgse[xx,0]*par_dir[xx,0] + bwmgse[xx,1]*par_dir[xx,1] + bwmgse[xx,2]*par_dir[xx,2]

  ;;Radial (outwards), azimuthal (eastwards), Bo coord system
  Brab = [[B_perp2],[B_perp1],[B_par]]


  store_data,'Brab',data = {x:etimes,y:Brab},$
             dlim={constant:[0],colors:[0]}

  ylim,['Brab'],-50,50
  options,'Brab','colors',[0,50,250]
  tplot,['Brab']


















  year = strmid(date, 0, 4)
  mm   = strmid(date, 5, 2)
  dd   = strmid(date, 8, 2)

  if keyword_set(hires) then $
     datafile = folder + rbx + 'efw-diagnostic_' + year + mm + dd + '_v' + vstr + '_hr.cdf' else $
        datafile = folder + rbx + 'efw-diagnostic_' + year + mm + dd+ '_v' + vstr + '.cdf'

  file_copy, skeletonFile, datafile, /overwrite ; Force to replace old file.
  cdfid = cdf_open(datafile)



;;--------------------------------------------------
;;spinfit_both_boompairs
;;--------------------------------------------------


;;...quantities to add to CDF file
;;Efield and Bfield in radial, azimuthal, FA coord. Include both the
;;nonsmoothed and smoothed (over 20 min) versions.



  stop


  cdf_varput,cdfid,'epoch',epoch

  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e12',transpose(spinfit_vxb_edotb_12)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e23',transpose(spinfit_vxb_edotb_23)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e23_m2',transpose(spinfit_vxb_edotb_23_m2)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e24',transpose(spinfit_vxb_edotb_24)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e24_m2',transpose(spinfit_vxb_edotb_24_m2)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_mgse_e34',transpose(spinfit_vxb_edotb_34)

  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e12',transpose(spinfit_vxb_coro_edotb_12)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e23',transpose(spinfit_vxb_coro_edotb_23)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e23_m2',transpose(spinfit_vxb_coro_edotb_23_m2)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e24',transpose(spinfit_vxb_coro_edotb_24)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e24_m2',transpose(spinfit_vxb_coro_edotb_24_m2)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_mgse_e34',transpose(spinfit_vxb_coro_edotb_34)

  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e12',transpose(Erab12)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e23',transpose(Erab23)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e23_m2',transpose(Erab23_m2)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e24',transpose(Erab24)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e24_m2',transpose(Erab24_m2)
  cdf_varput,cdfid,'Einertial_spinfit_edotb_raf_e34',transpose(Erab34)

  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e12',transpose(Erab12c)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e23',transpose(Erab23c)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e23_m2',transpose(Erab23c_m2)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e24',transpose(Erab24c)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e24_m2',transpose(Erab24c_m2)
  cdf_varput,cdfid,'Ecoro_spinfit_edotb_raf_e34',transpose(Erab34c)


  cdf_varput,cdfid,'bfield_raf',transpose(brab)
  cdf_varput,cdfid,'bfield_mgse',transpose(mag_mgse.y)
  cdf_varput,cdfid,'bfield_model_mgse',transpose(mag_model.y)
  cdf_varput,cdfid,'bfield_minus_model_mgse',transpose(mag_diff.y)
  cdf_varput,cdfid,'bfield_magnitude_minus_modelmagnitude',mag_diff_magnitude
  cdf_varput,cdfid,'bfield_magnitude',mag_data_magnitude



  
  ;;Remove the highcadence version and rename the lowcadence one to vsvy_vavg
  cdf_vardelete,cdfid,'vsvy_vavg'
  cdf_varrename,cdfid,'vsvy_vavg_lowcadence','vsvy_vavg'

  cdf_varput,cdfid,'flags_charging_bias_eclipse',transpose(flags)
  cdf_varput,cdfid,'flags_all',transpose(flag_arr)

  if is_struct(dens12) then cdf_varput,cdfid,'density_v12',dens12.y
  if is_struct(dens34) then cdf_varput,cdfid,'density_v34',dens34.y

  cdf_varput,cdfid,'vsvy_vavg',transpose([[sum12],[sum34],[sum56]])
  cdf_varput,cdfid,'VxB_mgse',transpose(vxb.y)

  cdf_varput,cdfid,'efield_coro_mgse',transpose(ecoro_mgse.y)
  cdf_varput,cdfid,'mlt',reform(mlt.y)
  cdf_varput,cdfid,'mlat',reform(mlat.y)
  cdf_varput,cdfid,'lshell',reform(lshell.y)
  cdf_varput,cdfid,'pos_gse',transpose(pos_gse.y)
  cdf_varput,cdfid,'vel_gse',transpose(vel_gse.y)
  cdf_varput,cdfid,'spinaxis_gse',transpose(sa.y)
  cdf_varput,cdfid,'angle_Ey_Ez_Bo',transpose(angles.y)



;variables to delete




  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e12'
  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e23'
  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e23_m2'
  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e24'
  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e24_m2'
  cdf_vardelete,cdfid,'Einertial_spinfit_mgse_e34'

  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e12'
  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e23'
  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e23_m2'
  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e24'
  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e24_m2'
  cdf_vardelete,cdfid,'Ecoro_spinfit_mgse_e34'



  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e12'
  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e23'
  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e23_m2'
  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e24'
  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e24_m2'
  cdf_vardelete,cdfid,'Esc_spinfit_mgse_e34'


  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e12'
  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e23'
  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e23_m2'
  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e24'
  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e24_m2'
  cdf_vardelete,cdfid,'Einertial_spinfit_raf_e34'

  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e12'
  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e23'
  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e23_m2'
  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e24'
  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e24_m2'
  cdf_vardelete,cdfid,'Ecoro_spinfit_raf_e34'






  cdf_vardelete,cdfid,'density_v12_hires'
  cdf_vardelete,cdfid,'density_v34_hires'
  cdf_vardelete,cdfid,'e_spinfit_mgse_efw_qual'
  ;; cdf_vardelete,cdfid,'efw_qual'
;     cdf_vardelete,cdfid,'efield_spinfit_mgse_e12'
;     cdf_vardelete,cdfid,'efield_spinfit_mgse_e34'
;     cdf_vardelete,cdfid,'efield_spinfit_mgse_e24'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_mgse_e12'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_mgse_e34'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_mgse_e24'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_coro_e12'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_coro_e34'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_coro_e24'

  cdf_vardelete,cdfid,'e_spinfit_mgse_BEB_config'
  cdf_vardelete,cdfid,'e_spinfit_mgse_DFB_config'
  cdf_vardelete,cdfid,'sigma12_spinfit_mgse'
  cdf_vardelete,cdfid,'sigma34_spinfit_mgse'
  cdf_vardelete,cdfid,'npoints12_spinfit_mgse'
  cdf_vardelete,cdfid,'npoints34_spinfit_mgse'
  cdf_vardelete,cdfid,'efield_uvw'
  cdf_vardelete,cdfid,'efield_raw_uvw'
;     cdf_vardelete,cdfid,'density'
  cdf_vardelete,cdfid,'vsvy'
  cdf_vardelete,cdfid,'esvy'

  cdf_vardelete,cdfid,'vsvy_vavg_combo'
  cdf_vardelete,cdfid,'mag_model_mgse'
  cdf_vardelete,cdfid,'mag_minus_model_mgse'
  cdf_vardelete,cdfid,'mag_spinfit_mgse'
;     cdf_vardelete,cdfid,'efield_spinfit_vxb_mgse'
  cdf_vardelete,cdfid,'vel_coro_mgse'
  cdf_vardelete,cdfid,'esvy_vxb_mgse'
  cdf_vardelete,cdfid,'efield_mgse'
  cdf_vardelete,cdfid,'vsvy_combo'

;     cdf_vardelete,cdfid,'e12_spinfit_mgse'
;     cdf_vardelete,cdfid,'e34_spinfit_mgse'
;     cdf_vardelete,cdfid,'e12_vxb_spinfit_mgse'
;     cdf_vardelete,cdfid,'e34_vxb_spinfit_mgse'
                                ;   cdf_vardelete,cdfid,'e12_vxb_coro_spinfit_mgse'
                                ;   cdf_vardelete,cdfid,'e34_vxb_coro_spinfit_mgse'

;  endif



  cdf_close, cdfid


end







     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='12',ql=ql,datatype='vsvy',method=1

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp_sf_12'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp_sf_vxb_12'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp_sf_vxb_coro_12'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_edotb_12'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_coro_edotb_12'

     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete


     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='34',ql=ql,/rerun,datatype='vsvy',method=1

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp_sf_34'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp_sf_vxb_34'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp_sf_vxb_coro_34'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_edotb_34'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_coro_edotb_34'

     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete


     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='23',ql=ql,/rerun,datatype='vsvy',method=1

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp_sf_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp_sf_vxb_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp_sf_vxb_coro_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_edotb_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_coro_edotb_23'

     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete



     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='23',ql=ql,/rerun,datatype='vsvy',method=2

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp2_sf_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp2_sf_vxb_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp2_sf_vxb_coro_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp2_sf_vxb_edotb_23'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp2_sf_vxb_coro_edotb_23'

     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete



     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='24',ql=ql,/rerun,datatype='vsvy',method=1

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp_sf_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp_sf_vxb_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp_sf_vxb_coro_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_edotb_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp_sf_vxb_coro_edotb_24'


     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete



     ;; rbsp_efw_edotb_to_zero_crib_v2,$
     ;;    date,sc,/no_spice_load,/noplot,suffix='edotb',boom_pair='24',ql=ql,/rerun,datatype='vsvy',method=2

     ;; copy_data,rbspx+'_efw_esvy_spinfit',$
     ;;           'tmp2_sf_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;           'tmp2_sf_vxb_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;           'tmp2_sf_vxb_coro_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;           'tmp2_sf_vxb_edotb_24'
     ;; copy_data,rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb',$
     ;;           'tmp2_sf_vxb_coro_edotb_24'


     ;; store_data,[rbspx+'_efw_esvy_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb',$
     ;;             rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb'],/delete









;; ;;--------------------------------------------------
;; ;;Now copy to final names
;;      copy_data,'tmp1_sf_12',$
;;                rbspx+'_efw_esvy_spinfit_12'
;;      copy_data,'tmp1_sf_34',$
;;                rbspx+'_efw_esvy_spinfit_34'
;;      copy_data,'tmp1_sf_23',$
;;                rbspx+'_efw_esvy_spinfit_23'
;;      copy_data,'tmp2_sf_23',$
;;                rbspx+'_efw_esvy_spinfit_23_m2'
;;      copy_data,'tmp1_sf_24',$
;;                rbspx+'_efw_esvy_spinfit_24'
;;      copy_data,'tmp2_sf_24',$
;;                rbspx+'_efw_esvy_spinfit_24_m2'


;;      copy_data,'tmp1_sf_vxb_12',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_12'
;;      copy_data,'tmp1_sf_vxb_34',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_34'
;;      copy_data,'tmp1_sf_vxb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23'
;;      copy_data,'tmp2_sf_vxb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_23_m2'
;;      copy_data,'tmp1_sf_vxb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24'
;;      copy_data,'tmp2_sf_vxb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_24_m2'


;;      copy_data,'tmp1_sf_vxb_coro_12',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_12'
;;      copy_data,'tmp1_sf_vxb_coro_34',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_34'
;;      copy_data,'tmp1_sf_vxb_coro_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23'
;;      copy_data,'tmp2_sf_vxb_coro_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_23_m2'
;;      copy_data,'tmp1_sf_vxb_coro_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24'
;;      copy_data,'tmp2_sf_vxb_coro_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_24_m2'


;;      copy_data,'tmp1_sf_vxb_edotb_12',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_12'
;;      copy_data,'tmp1_sf_vxb_edotb_34',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_34'
;;      copy_data,'tmp1_sf_vxb_edotb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23'
;;      copy_data,'tmp2_sf_vxb_edotb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_23_m2'
;;      copy_data,'tmp1_sf_vxb_edotb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24'
;;      copy_data,'tmp2_sf_vxb_edotb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_spinfit_edotb_24_m2'


;;      copy_data,'tmp1_sf_vxb_coro_edotb_12',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_12'
;;      copy_data,'tmp1_sf_vxb_coro_edotb_34',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_34'
;;      copy_data,'tmp1_sf_vxb_coro_edotb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23'
;;      copy_data,'tmp2_sf_vxb_coro_edotb_23',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_23_m2'
;;      copy_data,'tmp1_sf_vxb_coro_edotb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24'
;;      copy_data,'tmp2_sf_vxb_coro_edotb_24',$
;;                rbspx+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_edotb_24_m2'







;; if ~keyword_set(no_spice_load) then rbsp_load_spice_kernels
;; rbsp_load_state,probe=sc,/no_spice_load,datatype=['spinper','spinphase','mat_dsc','Lvec'] 
;; get_data,rbspx+'_spinper',data=sp
;; spinper = median(sp.y)


;; trange = timerange()
;; cp0 = rbsp_efw_get_cal_params(trange[0])
;; cp = cp0.a

;; boom_length = cp.boom_length
;; boom_shorting_factor = cp.boom_shorting_factor

;; split_vec,rbspx+'_efw_vsvy',suffix = '_'+['1','2','3','4','5','6']

;; ;;Determine which antenna pairs you'd like to use
;; tplot,rbspx+'_efw_vsvy_[1-4]'



;; ;;Decide which antenna pair you want to mimic
;; if pair eq '13' then pairm = '12'
;; if pair eq '14' then pairm = '12'
;; if pair eq '24' then pairm = '21' ;-->shift V4 by 1/4 spinperiod to mimic V1
;; if pair eq '23' then pairm = '21'

;; plane_dim = 0


;; pa0 = strmid(pair,0,1)
;; pb0 = strmid(pair,1,1)
;; pa1 = strmid(pairm,0,1)
;; pb1 = strmid(pairm,1,1)

;; get_data,rbspx+'_efw_vsvy_'+pb0,data=v

;; ;;adjust the timing on one of the booms so that we always effectively
;; ;;have a V12 or V21 boom

;; ;;To mimic boom 2, add 1/4 spinperiod to V
;; if pair eq '13' or pair eq '24' then v_adj = v.x - spinper/4.
;; if pair eq '14' or pair eq '23' then v_adj = v.x + spinper/4.


;; store_data,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,data={x:v_adj,y:v.y}
;; ;;here the adjusted V3 should be similar to V2

;; ylim,[rbspx+'_efw_vsvy_'+pb0,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,rbspx+'_efw_vsvy_'+pa0],-10,10
;; ;;The second and third plotted quantities should be 180 deg out of phase
;; tplot,[rbspx+'_efw_vsvy_'+pb0,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,rbspx+'_efw_vsvy_'+pa0]


;; store_data,'v'+pa0+'_v'+pb0+'_to_v'+pa0+'_v'+pb1+'_comb',$
;;            data=[rbspx+'_efw_vsvy_'+pa0,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1]
;; store_data,'v'+pa1+'_v'+pb1+'_comb',$
;;            data=[rbspx+'_efw_vsvy_'+pa1,rbspx+'_efw_vsvy_'+pb1]
;; options,'v'+pa0+'_v'+pb0+'_to_v'+pa0+'_v'+pb1+'_comb','colors',[0,250]
;; options,'v'+pa1+'_v'+pb1+'_comb','colors',[0,250]
;; ;;These should be similar
;; tplot,['v'+pa0+'_v'+pb0+'_to_v'+pa0+'_v'+pb1+'_comb',$
;;        'v'+pa1+'_v'+pb1+'_comb']


;; ;create electric field
;; dif_data,rbspx+'_efw_vsvy_'+pa0,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,$
;;          newname=rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pa0+pb1
;; dif_data,rbspx+'_efw_vsvy_'+pa1,rbspx+'_efw_vsvy_'+pb1,newname=rbspx+'_efw_esvy_'+pa1+pb1


;; ;;These two should be very similar
;; ylim,[rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pa0+pb1,rbspx+'_efw_esvy_'+pa1+pb1],0,0
;; tplot,[rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pa0+pb1,$
;;        rbspx+'_efw_esvy_'+pa1+pb1]

;; ;;Now give these quantities units of electric field (keep the same
;; ;;names). If the pair quantities are reversed, then multiply by -1
;; mult = 1.
;; if pairm eq '21' or pairm eq '31' or pairm eq '32' or pairm eq '41' or pairm eq '42' or pairm eq '43' then mult = -1.

;; get_data,rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pa0+pb1,data=d
;; data_att = {coord_sys:'uvw'}
;; dlim = {data_att:data_att}

;; if mult eq 1 then ename = rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pa0+pb1
;; if mult eq -1 then ename = rbspx+'_efw_esvy_'+pa0+pb0+'_to_'+pb1+pa0

;; store_data,ename,$
;;            data={x:d.x,y:mult*(1000./boom_length[0])*[[d.y],[d.y],[d.y]]},dlim=dlim

;; if mult eq 1 then enamem = rbspx+'_efw_esvy_'+pa1+pb1
;; if mult eq -1 then enamem = rbspx+'_efw_esvy_'+pb1+pa1

;; get_data,rbspx+'_efw_esvy_'+pa1+pb1,data=d
;; store_data,enamem,data={x:d.x,y:mult*(1000./boom_length[0])*[[d.y],[d.y],[d.y]]},dlim=dlim


;; ;;now "pair" is mimicking "pairm". Note that the above
;; ;;two don't have to be very similar since they're measuring two
;; ;;different full-cadence efields and have different offsets, etc.
;; ylim,[ename,enamem],0,0
;; tplot,[ename,enamem]

;; rbsp_spinfit,enamem, plane_dim=plane_dim
;; rbsp_spinfit,ename,plane_dim=plane_dim

;; ylim,[enamem+'_spinfit',ename+'_spinfit'],-50,50
;; tplot,[enamem+'_spinfit',ename+'_spinfit']

;; rbsp_cotrans,enamem+'_spinfit',enamem+'_spinfit_mgse',/dsc2mgse
;; rbsp_cotrans,ename+'_spinfit',ename+'_spinfit_mgse', /dsc2mgse


;; ylim,[enamem+'_spinfit_mgse',ename+'_spinfit_mgse'],-100,100
;; ylim,[rbspx+'_efw_vsvy_'+pa1,$
;;       rbspx+'_efw_vsvy_'+pb1,$
;;       rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1],-20,20


;; split_vec,enamem+'_spinfit_mgse'
;; split_vec,ename+'_spinfit_mgse'

;; ;; tplot,[ename,enamem,$
;; ;;        ename+'_spinfit_mgse_y',$
;; ;;        enamem+'_spinfit_mgse_y',$
;; ;;        rbspx+'_efw_vsvy_'+pa1,$
;; ;;        rbspx+'_efw_vsvy_'+pb1,$
;; ;;        rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1]
;; ;; tplot,[ename,enamem,$
;; ;;        ename+'_spinfit_mgse_z',$
;; ;;        enamem+'_spinfit_mgse_z',$
;; ;;        rbspx+'_efw_vsvy_'+pa1,$
;; ;;        rbspx+'_efw_vsvy_'+pb1,$
;; ;;        rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1]


;; ylim,ename+'_spinfit_mgse_y',-40,40
;; ylim,ename+'_spinfit_mgse_z',-40,40
;; ylim,rbspx+'_efw_vsvy_'+pa1,-10,0
;; ylim,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,-10,0

;; ;;Plot for Wygant
;; options,ename+'_spinfit_mgse_y','ytitle','RBSP-A EFW!Cspinfit MGSEy!C[mV/m]'
;; options,ename+'_spinfit_mgse_z','ytitle','RBSP-A EFW!Cspinfit MGSEz!C[mV/m]'
;; options,ename,'ytitle','RBSP-A EFW!CEsvy'+pa0+pb0+' UVW!C[mV/m]'
;; options,rbspx+'_efw_vsvy_'+pa1,'ytitle','RBSP-A EFW!CVsvy'+pa0+'!C[V]'
;; options,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,'ytitle','RBSP-A EFW!CVsvy'+pb0+'!C[V]'

;; options,ename+'_spinfit_mgse_y','labels',''
;; options,ename+'_spinfit_mgse_z','labels',''
;; options,rbspx+'_efw_vsvy_'+pa1,'labels',''
;; options,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,'labels',''

;; options,ename+'_spinfit_mgse_y','ysubtitle',''
;; options,ename+'_spinfit_mgse_z','ysubtitle',''
;; options,rbspx+'_efw_vsvy_'+pa1,'ysubtitle',''
;; options,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,'ysubtitle',''


;; ;;Plot for Wygant
;; tplot_options,'title',''
;; tplot,[ename+'_spinfit_mgse_y',$
;;        ename+'_spinfit_mgse_z',$
;;        ename,$
;;        rbspx+'_efw_vsvy_'+pa1,$
;;        rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1]

;; ylim,ename+'_spinfit_mgse_y',-10,20
;; ylim,ename+'_spinfit_mgse_z',-60,20
;; ylim,rbspx+'_efw_vsvy_'+pa1,-10,0
;; ylim,rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1,-10,0

;; tlimit,'2015-03-17/19:53','2015-03-17/20:07' 

;; tplot_options,'title',''
;; tplot,[ename+'_spinfit_mgse_y',$
;;        ename+'_spinfit_mgse_z',$
;;        ename,$
;;        rbspx+'_efw_vsvy_'+pa1,$
;;        rbspx+'_efw_vsvy_'+pb0+'_to_'+pb1]



;; ;;--------------------------------------------------
;; ;;--------------------------------------------------










