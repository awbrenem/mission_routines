;+
; NAME: rbsp_efw_burst_fa_rotate_crib
; SYNTAX: 
; PURPOSE: Rotate RBSP EFW burst data to field-aligned coordinates
; INPUT: 
; OUTPUT: 
; KEYWORDS: 
; HISTORY: Created by Aaron W Breneman, Univ. Minnesota  4/10/2014
; VERSION: 
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2014-04-16 13:53:43 -0500 (Wed, 16 Apr 2014) $
;   $LastChangedRevision: 14832 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/examples/rbsp_efw_burst_fa_rotate_crib.pro $
;-

	rbsp_efw_init

	date = '2013-01-17'   ;Wave at 4 Hz
	timespan,date
	
	t0 = date + '/03:15'
	t1 = date + '/04:13'

	probe='b'
	rbspx = 'rbsp'+probe


;-----------
	;isolate a selected burst 
	;(CURRENTLY NEED AT LEAST 1 SEC OF BURST FOR MGSE TRANSFORMATION TO WORK!!!)

	;times for an RBSP-A burst
;	t0z = time_double(date + '/03:22:37.000')
;	t1z = time_double(date + '/03:22:39.000')

	;times for an RBSP-B burst
	t0z = time_double(date + '/03:43:32.200')
	t1z = time_double(date + '/03:43:34.000')
;-----------





	dt = time_double(t1) - time_double(t0)



	;Make tplot plots looks pretty
	charsz_plot = 0.8  ;character size for plots
	charsz_win = 1.2  
	!p.charsize = charsz_win
	tplot_options,'xmargin',[20.,15.]
	tplot_options,'ymargin',[3,6]
	tplot_options,'xticklen',0.08
	tplot_options,'yticklen',0.02
	tplot_options,'xthick',2
	tplot_options,'ythick',2
	tplot_options,'labflag',-1	

 
;--------------------------------------------------------------------------------
;Find the GSE coordinates of the sc spin axis. This will be used to transform the 
;Mag data from GSE -> MGSE coordinates
;--------------------------------------------------------------------------------


	;Load spice kernels and get antenna pointing direction
	rbsp_load_spice_kernels

	rbsp_efw_position_velocity_crib,/no_spice_load,/noplot
	get_data,rbspx+'_spinaxis_direction_gse',data=wsc_GSE	



;------------------------------------------------------
;Get EMFISIS DC mag data in GSE
;------------------------------------------------------

	;Load EMFISIS data (defaults to 'hires', but can also choose '1sec' or '4sec')
	rbsp_load_emfisis,probe=probe,coord='gse',cadence='hires',level='l3'
	;rbsp_load_emfisis,probe=probe,coord='gse',cadence='4sec',level='l3'


date = '2014-03-26'
timespan,date
rbsp_load_emfisis_quicklook,'a',datatype='uvw'
;pro rbsp_load_emfisis_quicklook, probe = probe, datatype = datatype $
;  , _extra = _extra, use_local = use_local, downloadonly = downloadonly


	;Transform the Mag data to MGSE coordinates
	get_data,rbspx+'_emfisis_l3_hires_gse_Mag',data=tmpp
	;get_data,rbspx+'_emfisis_l3_4sec_gse_Mag',data=tmpp

	wsc_GSE_tmp = [[interpol(wsc_GSE.y[*,0],wsc_GSE.x,tmpp.x)],$
				   [interpol(wsc_GSE.y[*,1],wsc_GSE.x,tmpp.x)],$
				   [interpol(wsc_GSE.y[*,2],wsc_GSE.x,tmpp.x)]]

	rbsp_gse2mgse,rbspx+'_emfisis_l3_hires_gse_Mag',reform(wsc_GSE_tmp),newname=rbspx+'_Mag_mgse'
	;rbsp_gse2mgse,rbspx+'_emfisis_l3_4sec_gse_Mag',reform(wsc_GSE_tmp),newname=rbspx+'_Mag_mgse'




;----------------------------------------------------------
;Get Esvy data in MGSE 
;----------------------------------------------------------

	timespan,t0,dt,/seconds


	;Load Esvy data in MGSE 
	;rbsp_load_efw_esvy_mgse,probe=probe,/no_spice_load

	rbsp_load_efw_waveform_partial,probe=probe,type='calibrated',datatype=['mscb2']
	rbsp_load_efw_waveform_partial,probe=probe,type='calibrated',datatype=['vb2']

	;plot all burst data
	tplot,[rbspx+'_efw_mscb2',rbspx+'_efw_vb2']
stop


	;plot zoomed-in time
	tplot,[rbspx+'_efw_mscb2',rbspx+'_efw_vb2']
	tlimit,t0z,t1z

stop

	;Reduce data to selected times
	dat1 = tsample(rbspx+'_efw_mscb2',[t0z,t1z],times=t1)
	store_data,rbspx+'_efw_mscb2_uvw',data={x:t1,y:dat1}
	dat1 = tsample(rbspx+'_efw_vb2',[t0z,t1z],times=t1)
	store_data,rbspx+'_efw_vb2_uvw',data={x:t1,y:dat1}
	dat1 = tsample(rbspx+'_Mag_mgse',[t0z,t1z],times=t1)
	store_data,rbspx+'_Mag_mgse',data={x:t1,y:dat1}
	
	
	;Create E-field variables (mV/m)
	trange = timerange()
	cp0 = rbsp_efw_get_cal_params(trange[0])

	if probe eq 'a' then cp = cp0.a else cp = cp0.b


	boom_length = cp.boom_length
	boom_shorting_factor = cp.boom_shorting_factor

	get_data,rbspx+'_efw_vb2_uvw',data=dd
	e12 = 1000.*(dd.y[*,0]-dd.y[*,1])/boom_length[0]
	e34 = 1000.*(dd.y[*,2]-dd.y[*,3])/boom_length[1]
	e56 = 1000.*(dd.y[*,4]-dd.y[*,5])/boom_length[2]
	
	eb = [[e12],[e34],[e56]]
	store_data,rbspx+'_efw_eb2_uvw',data={x:dd.x,y:eb}
	
	tplot,[rbspx+'_efw_eb2_uvw',rbspx+'_efw_mscb2_uvw']
	
	;Convert from UVW (spinning sc) to MGSE coord
	rbsp_uvw_to_mgse,probe,rbspx+'_efw_mscb2_uvw',/no_spice_load,/nointerp,/no_offset	
	rbsp_uvw_to_mgse,probe,rbspx+'_efw_eb2_uvw',/no_spice_load,/nointerp,/no_offset	
	
	copy_data,rbspx+'_efw_eb2_uvw_mgse',rbspx+'_efw_eb2_mgse'
	copy_data,rbspx+'_efw_mscb2_uvw_mgse',rbspx+'_efw_mscb2_mgse'

	tplot,[rbspx+'_efw_eb2_mgse',rbspx+'_efw_mscb2_mgse']

	split_vec,rbspx+'_efw_eb2_mgse'
	split_vec,rbspx+'_efw_mscb2_mgse'

	;Check to see how things look (MGSEx is spin axis)
	tplot,[rbspx+'_efw_eb2_mgse_x',rbspx+'_efw_eb2_mgse_y',rbspx+'_efw_eb2_mgse_z']
stop
	tplot,[rbspx+'_efw_mscb2_mgse_x',rbspx+'_efw_mscb2_mgse_y',rbspx+'_efw_mscb2_mgse_z']
stop



;---------------------------------------------------------
;Align searchcoil Bw to FA coord
;---------------------------------------------------------


	tplot,[rbspx+'_efw_mscb2_mgse',rbspx+'_Mag_mgse']

	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_mscb2_mgse',rbspx+'_Mag_mgse')

	copy_data,'theta_kb',rbspx+'_theta_kb'
	copy_data,'dtheta_kb',rbspx+'_dtheta_kb'
	copy_data,'minvar_eigenvalues',rbspx+'_minvar_eigenvalues'
	copy_data,'emax2eint',rbspx+'_emax2eint'
	copy_data,'eint2emin',rbspx+'_eint2emin'
  	copy_data,'emax_vec_minvar',rbspx+'_emax_vec_minvar'
  	copy_data,'eint_vec_minvar',rbspx+'_eint_vec_minvar'
  	copy_data,'emin_vec_minvar',rbspx+'_emin_vec_minvar'

	store_data,['theta_kb','dtheta_kb','minvar_eigenvalues','emax2eint','eint2emin','emax_vec_minvar','eint_vec_minvar','emin_vec_minvar'],/delete

	split_vec,rbspx+'_efw_mscb2_mgse_FA_minvar'
	ylim,[rbspx+'_efw_mscb2_mgse_FA_minvar_x',rbspx+'_efw_mscb2_mgse_FA_minvar_y',rbspx+'_efw_mscb2_mgse_FA_minvar_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_mscb2_mgse_FA_minvar_x',rbspx+'_efw_mscb2_mgse_FA_minvar_y',rbspx+'_efw_mscb2_mgse_FA_minvar_z']

	;plot wave normal angle from min variance analysis
	options,rbspx+'_theta_kb','ytitle',rbspx+ '!CWave normal!Cangle'
	options,rbspx+'_dtheta_kb','ytitle',rbspx+'!CWave normal!Cangle!Cuncertainty'
	
	tplot,[rbspx+'_theta_kb',rbspx+'_dtheta_kb',rbspx+'_minvar_eigenvalues',rbspx+'_emax2eint',rbspx+'_eint2emin']
	tplot,[rbspx+'_theta_kb',rbspx+'_dtheta_kb',rbspx+'_efw_mscb2_mgse_FA_minvar_x',rbspx+'_efw_mscb2_mgse_FA_minvar_y',rbspx+'_efw_mscb2_mgse_FA_minvar_z']

stop

	;------------------
	;Now try EFA coord
;	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_mscb2_mgse',rbspx+'_Mag_mgse',/efa)
;
;	split_vec,rbspx+'_efw_mscb2_mgse_EFA_coord'
;	ylim,[rbspx+'_efw_mscb2_mgse_EFA_coord_x',rbspx+'_efw_mscb2_mgse_EFA_coord_y',rbspx+'_efw_mscb2_mgse_EFA_coord_z'],-0.15,0.15
;
;	;Plot FA coord. z-hat is field direction
;	tplot,[rbspx+'_efw_mscb2_mgse_EFA_coord_x',rbspx+'_efw_mscb2_mgse_EFA_coord_y',rbspx+'_efw_mscb2_mgse_EFA_coord_z']



;---------------------------------------------------------
;Align Ew to FA coord
;---------------------------------------------------------


	tplot,[rbspx+'_efw_eb2_mgse',rbspx+'_Mag_mgse']


	;First try minimum variance option
	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_eb2_mgse',rbspx+'_Mag_mgse')

	store_data,['theta_kb','dtheta_kb','minvar_eigenvalues','emax2eint','eint2emin','emax_vec_minvar','eint_vec_minvar','emin_vec_minvar'],/delete


	split_vec,rbspx+'_efw_eb2_mgse_FA_minvar'
	ylim,[rbspx+'_efw_eb2_mgse_FA_minvar_x',rbspx+'_efw_eb2_mgse_FA_minvar_y',rbspx+'_efw_eb2_mgse_FA_minvar_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_eb2_mgse_FA_minvar_x',rbspx+'_efw_eb2_mgse_FA_minvar_y',rbspx+'_efw_eb2_mgse_FA_minvar_z']

stop

;	;------------------
;	;Now try EFA coord
;	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_eb2_mgse',rbspx+'_Mag_mgse',/efa)

;	split_vec,rbspx+'_efw_eb2_mgse_EFA_coord'
;	ylim,[rbspx+'_efw_eb2_mgse_EFA_coord_x',rbspx+'_efw_eb2_mgse_EFA_coord_y',rbspx+'_efw_eb2_mgse_EFA_coord_z'],-0.15,0.15
;
;	;Plot FA coord. z-hat is field direction
;	tplot,[rbspx+'_efw_eb2_mgse_EFA_coord_x',rbspx+'_efw_eb2_mgse_EFA_coord_y',rbspx+'_efw_eb2_mgse_EFA_coord_z']




;------------------------------------------------------------------------------
;Use Chasten's program to calculate spec, wave normal angle, ellipticity, etc..
;------------------------------------------------------------------------------


	twavpol,rbspx+'_efw_mscb2_mgse_FA_minvar',prefix=rbspx+'_efw_mscb2_mgse_FA_minvar'

	;change wave normal angle to degrees
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',data=wn
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',data={x:wn.x,y:wn.y/!dtor,v:wn.v}
	options,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle','spec',1
	zlim,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',0,90,0

	ylim,[rbspx+'_efw_mscb2_mgse_FA_minvar_powspec',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_degpol',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_elliptict',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_helicit',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_pspec3'],100,8000,1

	zlim,rbspx+'_efw_mscb2_mgse_FA_minvar_powspec',1d-9,1d-4,1
	zlim,rbspx+'_efw_mscb2_mgse_FA_minvar_pspec3',1d-9,1d-4,1

	;eliminate data under a certain deg of polarization threshold
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_degpol',data=degp
	goo = where(degp.y le 0.5)
	if goo[0] ne -1 then degp.y[goo] = !values.f_nan
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_degpol',data=degp
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_powspec',data=tmp
	if goo[0] ne -1 then tmp.y[goo] = !values.f_nan
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_powspec',data=tmp
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',data=tmp
	if goo[0] ne -1 then tmp.y[goo] = !values.f_nan
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',data=tmp
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_elliptict',data=tmp
	if goo[0] ne -1 then tmp.y[goo] = !values.f_nan
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_elliptict',data=tmp
	get_data,rbspx+'_efw_mscb2_mgse_FA_minvar_helict',data=tmp
	if goo[0] ne -1 then tmp.y[goo] = !values.f_nan
	store_data,rbspx+'_efw_mscb2_mgse_FA_minvar_helict',data=tmp

	tplot,[rbspx+'_efw_mscb2_mgse_FA_minvar_powspec',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_degpol',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_waveangle',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_elliptict',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_helict',$
			rbspx+'_efw_mscb2_mgse_FA_minvar_pspec3']

;----------------------------------------------------
;Plot FFTs
;----------------------------------------------------

	;Call program that plots FFTs, hodograms
	
	fftlims =  {xlog:1,xrange:[10,4000]} ;limits for the FFT
	hodlims = {xrange:[-20,20],yrange:[-20,20]}
		
	plot_wavestuff,rbspx+'_efw_eb2_mgse_FA_minvar',/psd,/hod,$
		/nodelete,$    ;keep tplot variables around
		extra_psd=fftlims,$
		extra_hod=hodlims

stop

;-----------------------------------------------------
;PLOT CROSS-CORRELATIONS
;-----------------------------------------------------
	
var1 = 'rbspb_efw_mscb2_mgse_FA_minvar'
;var2 = 'rbspb_efw_mscb2_mgse_FA_minvar'
var2 = 'rbspb_efw_eb2_mgse_FA_minvar'

vardim1 = 1  ;use this y-dim for var1 (x=0, y=1, z=2)
vardim2 = 2	 ;use this y-dim for var2 (x=0, y=1, z=2)


;get both variables on the same time cadence
tinterpol_mxn,var2,var1,newname=var2 + '_for_correlation',/spline
copy_data,var1,var1 + '_for_correlation'

var1 = var1 + '_for_correlation'
var2 = var2 + '_for_correlation'


;The following calculates the cross-phase and coherence between Ey-mgse and Ez-mgse 
;The function cross_spec_tplot return a multiple dimention array with the following structure. [[E_FREQUENCE_coordinate],[PHASE_E1_E2],[COHERENCE_E1_E2],[Aver_Pow_E_1],[Aver_Pow_E_2]]
Results=cross_spec_tplot(var1,vardim1,var2,vardim2,t0z,t1z,sub_interval=3,overlap_index=4)
Time_rbsp=strmid(time_string(t0z),0,10)+'_'+strmid(time_string(t0z),11,2)+strmid(time_string(t0z),14,2)+'UT'+'_to_'+strmid(time_string(t1z),0,10)+'_'+strmid(time_string(t1z),11,2)+strmid(time_string(t1z),14,2)+'UT'

freqrange = [0,4000]   ;Hz

;Plot the cross-phase and coherence in a .ps file
;Popen,rbspx+'_'+timerbsp
!p.multi = [0,0,4]
!p.charsize = 2
Plot,Results[*,0],Results[*,2],xtitle='f, Hz', ytitle='Coherence_Ey_Ez',title=rbspx+time_rbsp,xrange=freqrange
Plot,Results[*,0],Results[*,1]*180./3.14,xtitle='f, Hz', ytitle='Phase_Ey_Ez',title=rbspx+time_rbsp,xrange=freqrange
Plot,Results[*,0],Results[*,3],xtitle='f, Hz', ytitle='Power_Ey,mV^2/Hz',title=rbspx+time_rbsp,xrange=freqrange
Plot,Results[*,0],Results[*,4],xtitle='f, Hz', ytitle='Phase_Ez,MV^2/Hz',title=rbspx+time_rbsp,xrange=freqrange
;Pclose


end
