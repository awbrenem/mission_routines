;+
; NAME: rbsp_efw_survey_fa_rotate_crib
; SYNTAX: 
; PURPOSE: Rotate RBSP EFW survey data to field-aligned coordinates
; INPUT: 
; OUTPUT: 
; KEYWORDS: 
; HISTORY: Created by Aaron W Breneman, Univ. Minnesota  4/10/2014
; VERSION: 
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2014-04-10 14:44:47 -0500 (Thu, 10 Apr 2014) $
;   $LastChangedRevision: 14802 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/examples/rbsp_efw_burst_fa_rotate_crib.pro $
;-

	rbsp_efw_init

	date = '2013-01-17'   ;Wave at 4 Hz
	t0 = date + '/00:00'
	t1 = date + '/24:00'

	dt = time_double(t1) - time_double(t0)
	timespan,t0,dt,/sec


	probe='a'
	rbspx = 'rbsp'+probe

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

	;Load EMFISIS '4sec' data
	rbsp_load_emfisis,probe=probe,coord='gse',cadence='4sec',level='l3'


	;Transform the Mag data to MGSE coordinates
	get_data,rbspx+'_emfisis_l3_4sec_gse_Mag',data=tmpp

	wsc_GSE_tmp = [[interpol(wsc_GSE.y[*,0],wsc_GSE.x,tmpp.x)],$
				   [interpol(wsc_GSE.y[*,1],wsc_GSE.x,tmpp.x)],$
				   [interpol(wsc_GSE.y[*,2],wsc_GSE.x,tmpp.x)]]

	rbsp_gse2mgse,rbspx+'_emfisis_l3_4sec_gse_Mag',reform(wsc_GSE_tmp),newname='Mag_mgse'


;----------------------------------------------------------
;Get Esvy data in MGSE 
;----------------------------------------------------------


	;Load Esvy data in MGSE 
	rbsp_load_efw_esvy_mgse,probe=probe,/no_spice_load

	;plot all burst data
	tplot,[rbspx+'_efw_esvy_mgse']
stop

	;isolate a selected burst (CURRENTLY NEED AT LEAST 1 SEC OF BURST FOR MGSE TRANSFORMATION TO WORK!!!)
;	t0z = time_double(date + '/08:25:30')
;	t1z = time_double(date + '/08:25:34')
;	t0z = time_double(date + '/03:55:31.900')
;	t1z = time_double(date + '/03:55:33')
	t0z = time_double(date + '/03:22:37.000')
	t1z = time_double(date + '/03:22:39.000')

	;plot zoomed-in time
	tplot,[rbspx+'_efw_mscb2',rbspx+'_efw_vb2']
	tlimit,t0z,t1z

stop

	;Reduce data to selected times
	dat1 = tsample(rbspx+'_efw_mscb2',[t0z,t1z],times=t1)
	store_data,rbspx+'_efw_mscb2_uvw',data={x:t1,y:dat1}
	dat1 = tsample(rbspx+'_efw_vb2',[t0z,t1z],times=t1)
	store_data,rbspx+'_efw_vb2_uvw',data={x:t1,y:dat1}
	dat1 = tsample('Mag_mgse',[t0z,t1z],times=t1)
	store_data,'Mag_mgse',data={x:t1,y:dat1}
	
	
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


	tplot,[rbspx+'_efw_mscb2_mgse','Mag_mgse']

	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_mscb2_mgse','Mag_mgse')

	split_vec,rbspx+'_efw_mscb2_mgse_FA_minvar'
	ylim,[rbspx+'_efw_mscb2_mgse_FA_minvar_x',rbspx+'_efw_mscb2_mgse_FA_minvar_y',rbspx+'_efw_mscb2_mgse_FA_minvar_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_mscb2_mgse_FA_minvar_x',rbspx+'_efw_mscb2_mgse_FA_minvar_y',rbspx+'_efw_mscb2_mgse_FA_minvar_z']

stop

	;------------------
	;Now try EFA coord
	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_mscb2_mgse','Mag_mgse',/efa)

	split_vec,rbspx+'_efw_mscb2_mgse_EFA_coord'
	ylim,[rbspx+'_efw_mscb2_mgse_EFA_coord_x',rbspx+'_efw_mscb2_mgse_EFA_coord_y',rbspx+'_efw_mscb2_mgse_EFA_coord_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_mscb2_mgse_EFA_coord_x',rbspx+'_efw_mscb2_mgse_EFA_coord_y',rbspx+'_efw_mscb2_mgse_EFA_coord_z']


stop

;---------------------------------------------------------
;Align Ew to FA coord
;---------------------------------------------------------


	tplot,[rbspx+'_efw_eb2_mgse','Mag_mgse']


	;First try minimum variance option
	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_eb2_mgse','Mag_mgse')

	split_vec,rbspx+'_efw_eb2_mgse_FA_minvar'
	ylim,[rbspx+'_efw_eb2_mgse_FA_minvar_x',rbspx+'_efw_eb2_mgse_FA_minvar_y',rbspx+'_efw_eb2_mgse_FA_minvar_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_eb2_mgse_FA_minvar_x',rbspx+'_efw_eb2_mgse_FA_minvar_y',rbspx+'_efw_eb2_mgse_FA_minvar_z']

stop

	;------------------
	;Now try EFA coord
	fa = rbsp_rotate_field_2_vec(rbspx+'_efw_eb2_mgse','Mag_mgse',/efa)

	split_vec,rbspx+'_efw_eb2_mgse_EFA_coord'
	ylim,[rbspx+'_efw_eb2_mgse_EFA_coord_x',rbspx+'_efw_eb2_mgse_EFA_coord_y',rbspx+'_efw_eb2_mgse_EFA_coord_z'],-0.15,0.15

	;Plot FA coord. z-hat is field direction
	tplot,[rbspx+'_efw_eb2_mgse_EFA_coord_x',rbspx+'_efw_eb2_mgse_EFA_coord_y',rbspx+'_efw_eb2_mgse_EFA_coord_z']



stop

end
