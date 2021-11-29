;+
; NAME:
;   rbsp_efw_make_l2_vsvy_hires
;
; PURPOSE: Create L2 (level-2) CDF files from the Van Allen Probe EFW
;			Vsvy data product.
;
;
; CALLING SEQUENCE:
;
; ARGUMENTS:	sc -> Which probe? 'a' or 'b'
;				date -> ex, '2013-10-13'
;
; KEYWORDS:
;
;
; EXAMPLES:	rbsp_efw_make_l2_vsvy_hires,'a','2012-10-13'
;
; NOTES: This program stuffs the following quantities into the skeleton CDF file.
;		V1-V6 -> single-ended potential quantities
;		(Vx+Vy)/2 -> opposing boom averages
;
;
; HISTORY:
;   March 2020: Created by Aaron Breneman, University of Minnesota
;
;
; VERSION:
;	$LastChangedBy: aaronbreneman $
;	$LastChangedDate: 2020-07-08 10:38:26 -0500 (Wed, 08 Jul 2020) $
;	$LastChangedRevision: 28864 $
; 	$URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/cdf_file_production/rbsp_efw_make_l2_vsvy_hires.pro $
;-




pro rbsp_efw_make_l2_vsvy_hires,sc,date,$
	folder=folder,$
	version=version,$
	testing=testing


	rbsp_efw_init
	compile_opt IDL2  ;make IDL more friendly


	timespan,date





	if n_elements(version) eq 0 then version = 1
	vstr = string(version, format='(I02)')
	version = 'v'+vstr


	sc=strlowcase(sc)
	if sc ne 'a' and sc ne 'b' then begin
		dprint,'Invalid spacecraft: '+sc+', returning.'
		return
	endif
	rbx = 'rbsp' + strlowcase(sc[0]) + '_'



;	;Folder for locally stored CDF files
;	if ~keyword_set(folder) then folder = !rbsp_efw.local_data_dir + $
;									'rbsp' + strlowcase(sc[0]) + path_sep() + $
;									'l2' + path_sep() + $
;									'spinfit' + path_sep() + $
;									strmid(date,0,4) + path_sep()
;	; make sure we have the trailing slash on folder
;	if strmid(folder,strlen(folder)-1,1) ne path_sep() then folder=folder+path_sep()
;	if ~keyword_set(no_cdf) then file_mkdir, folder


	homedir = (file_search('~',/expand_tilde))[0]+'/'



	;Grab skeleton file
	if ~keyword_set(testing) then begin
		folder = '/Volumes/UserA/user_homes/rbsp_efw/Code/tdas_svn_daily/general/missions/rbsp/efw/cdf_file_production/'
		skeletonfile=file_search(folder + rbx+'efw-lX_00000000_vXX.cdf',count=found)
		if ~found then begin
	     dprint,'Could not find skeleton CDF, returning.'
	     return
	  endif
	endif else begin		
		folder = homedir+'Desktop/code/Aaron/RBSP/TDAS_trunk_svn/idl/general/missions/rbsp/efw/cdf_file_production/'
		skeletonfile = folder +rbx+'efw-lX_00000000_vXX.cdf'
	endelse
	





	;Load the survey data
	rbsp_load_efw_waveform,probe=sc,type='calibrated',datatype='vsvy',/noclean
	get_data,rbx+'efw_vsvy',data=vsvy
	if ~is_struct(vsvy) then begin
  		dprint,rbx+'efw_vsvy unavailable, returning.'
  		return
	endif
	times_v = vsvy.x



;	;Not ready to release the V5, V6 values.
;	vsvy.y[*,4] = -1.0e31
;	vsvy.y[*,5] = -1.0e31


;Define keyword inheritance to pass to subroutines. This will ensure that
;subsequent routines don't reload spice kernels or rbsp_efw_init
;extra = create_struct('no_rbsp_efw_init',1)



;------------------------------------------------------------------
;Load ephemeris data and put on a once/min cadence
;------------------------------------------------------------------

	;Load SPICE data
	rbsp_load_spice_cdf_file,sc



	;Put ephemeris on a once/sec cadence 
	newt = dindgen(1440)*60. + time_double(date)




	;Grab official times for once/min values
	times_ephem = newt
	epoch_ephem = tplot_time_to_epoch(times_ephem,/epoch16)


	;interpolate ephemeris values to this once/min cadence
	tinterpol_mxn,rbx+'state_pos_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbx+'state_vel_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbx+'spinaxis_direction_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbx+'state_mlt',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbx+'state_mlat',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbx+'state_lshell',times_ephem,/overwrite,/spline

	get_data,rbx+'state_pos_gse',data=pos_gse
	get_data,rbx+'state_vel_gse',data=vel_gse
	get_data,rbx+'spinaxis_direction_gse',data=wgse
	get_data,rbx+'state_mlt',data=mlt
	get_data,rbx+'state_mlat',data=mlat
	get_data,rbx+'state_lshell',data=lshell






  ;--------------------------------------------------
  ;Get flag values (also gets density values from v12 and v34)
  ;--------------------------------------------------



   flag_str = rbsp_efw_get_flag_values(sc,times_ephem,density_min=dmin,boom_pair=bp,_extra=extra,flag_names=fn)
   flag_arr = flag_str.flag_arr
   bias_sweep_flag = flag_str.bias_sweep_flag
   ab_flag = flag_str.ab_flag
   charging_flag = flag_str.charging_flag




;	;not ready to release V5 or V6 values
;	v5v6[*] = -1.0e31


;;*****TEMPORARY********
;;Set bad vavg values to ISTP compliant value of -1.0E31
;goo = where((finite(vsvy.y[*,0]) eq 0) or (finite(vsvy.y[*,1]) eq 0))
;if goo[0] ne -1 then v1v2[goo] = -1.0e31
;goo = where((finite(vsvy.y[*,2]) eq 0) or (finite(vsvy.y[*,3]) eq 0))
;if goo[0] ne -1 then v3v4[goo] = -1.0e31
;;**********************





	;Get official 16 S/sec times
	epoch_v = tplot_time_to_epoch(times_v,/epoch16)



;;;charging, autobias, eclipse, and extreme charging flags all in one variable for convenience
;flags = [[flag_arr[*,15]],[flag_arr[*,14]],[flag_arr[*,1]],[flag_arr[*,16]]]


	;Rename generic skeleton file with date 
	datafile = folder+rbx+'efw-l2_vsvy-hires_'+strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)+'_v'+vstr+'.cdf'
	file_copy, skeletonFile, datafile, /overwrite ; Force to replace old file.
	cdfid = cdf_open(datafile)


stop



	;Final list of variables to NOT delete
	varsave_general = ['epoch','epoch_v','global_flag',$
					   'vsvy',$
					   'flags_all',$
					   'mlt','mlat','lshell',$
					   'position_gse','velocity_gse',$
					   'spinaxis_gse']
;	'flags_charging_bias_eclipse']



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
		tstt = array_contains(varsave_general,CDFvarnames[qq]) & $
		if not tstt then print,'Deleting var:  ', CDFvarnames[qq] & $
		if not tstt then cdf_vardelete,cdfid,CDFvarnames[qq]
	endfor




	;Stuff values into CDF file
	cdf_varput,cdfid,'epoch',epoch_ephem
	cdf_varput,cdfid,'epoch_v',epoch_v
	cdf_varput,cdfid,'vsvy',transpose(vsvy.y)



	
	cdf_varput,cdfid,'flags_all',transpose(flag_arr)
  	cdf_varput,cdfid,'global_flag',reform(flag_arr[*,0])
	cdf_varput,cdfid,'mlt',transpose(mlt.y)
	cdf_varput,cdfid,'mlat',transpose(mlat.y)
	cdf_varput,cdfid,'lshell',transpose(lshell.y)
	cdf_varput,cdfid,'position_gse',transpose(pos_gse.y)
	cdf_varput,cdfid,'velocity_gse',transpose(vel_gse.y)
	cdf_varput,cdfid,'spinaxis_gse',transpose(wgse.y)


stop

	;Close out...
	cdf_close, cdfid
	dprint,'END TIME IS: ',systime()
	store_data,tnames(),/delete

end
