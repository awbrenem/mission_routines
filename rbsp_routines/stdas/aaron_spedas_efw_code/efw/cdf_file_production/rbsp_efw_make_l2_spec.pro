;+
;rbsp_efw_spec_L2_crib
;
;Loads and plots RBSP (Van Allen probes) spectral data
;used for the L2 CDF files
;
; SPEC returns 7 channels, with nominal data selection:
;		SPEC0: E12AC
;		SPEC1: E56AC
;		SPEC2: SCMpar
;		SPEC3: SCMperp
;		SPEC4: SCMW
;		SPEC5: V1AC
;		SPEC6: V2AC
;
;			Select 7 of: E12dc,E34dc,E56dc
;						 E12ac,E34ac,E56ac
;						 Edcpar,Edcprp
;						 Eacpar,Eacprp
;						 V1ac,V2ac,V3ac,V4ac,V5ac,V6ac
;						 SCMU,SCMV,SCMW
;						 SCMpar,SCMprp,
;						 (V1dc+V2dc+V3dc+V4dc)/4,
;						 Edcprp2, Eacprp2, SCMprp2
;
;	notes:
;
;
;	Aaron Breneman, UMN, Feb 2013
;	email: awbrenem@gmail.com
;
; VERSION:
;	$LastChangedBy: aaronbreneman $
;	$LastChangedDate: 2020-07-08 10:38:26 -0500 (Wed, 08 Jul 2020) $
;	$LastChangedRevision: 28864 $
;	$URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/cdf_file_production/rbsp_efw_make_l2_spec.pro $
;
;-



pro rbsp_efw_make_l2_spec,sc,date,$
	folder=folder,$
	testing=testing,$
	boom_pair=bp,$
	version=version


	timespan, date


	if ~KEYWORD_SET(bp) then bp = '12'



	if n_elements(version) eq 0 then version = 3
	vstr = string(version, format='(I02)')



	sc=strlowcase(sc)
	if sc ne 'a' and sc ne 'b' then begin
		dprint,'Invalid spacecraft: '+sc+', returning.'
		return
	endif
	rbspx = 'rbsp'+sc




    homedir = (file_search('~',/expand_tilde))[0]+'/'


	;Grab skeleton file
	if ~keyword_set(testing) then begin
		folder = '/Volumes/UserA/user_homes/rbsp_efw/Code/tdas_svn_daily/general/missions/rbsp/efw/cdf_file_production/'
		skeletonfile=file_search(folder + rbspx+'efw-l2_spec_00000000_v'+vstr+'.cdf',count=found)
		if ~found then begin
	     dprint,'Could not find skeleton CDF, returning.'
	     return
	  endif
	endif else begin		
		folder = homedir+'Desktop/code/Aaron/RBSP/TDAS_trunk_svn/idl/general/missions/rbsp/efw/cdf_file_production/'
		skeletonfile = folder +rbspx+'_efw-l2_spec_00000000_v'+vstr+'.cdf'
	endelse


	file_mkdir,folder




;	; Grab the skeleton file.
;	skeleton=rbspx+'/l2/spec/0000/'+ $
;		rbspx+'_efw-l2_spec_00000000_v'+vstr+'.cdf'
;	if ~keyword_set(testing) then source_file=file_retrieve(skeleton,_extra=!rbsp_efw)
;
;
 ; if keyword_set(testing) then begin
  ;   skeleton = 'rbspa_efw-l2_spec_00000000_v'+vstr+'.cdf'
   ;  source_file='~/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/idl/general/missions/rbsp/efw/cdf_file_production/' + skeleton
  ;endif




;	source_file=file_retrieve(skeleton,_extra=!rbsp_efw)
	;source_file=file_retrieve(source_file,_extra=!rbsp_efw)

	; use skeleton from the staging dir until we go live in the main data tree
	;source_file='/Volumes/DataA/user_volumes/kersten/data/rbsp/'+skeleton


	; make sure we have the skeleton CDF
	source_file=file_search(skeletonfile,count=found) ; looking for single file, so count will return 0 or 1
	if ~found then begin
		dprint,'Could not find spec v'+vstr+' skeleton CDF, returning.'
		return
	endif
	source_file=source_file[0]



	;Load the spectrogram data
	rbsp_load_efw_spec,probe=sc,type='calibrated'




	get_data,rbspx+'_efw_64_spec0',data=spec0,dlimits=dlim0
	get_data,rbspx+'_efw_64_spec1',data=spec1,dlimits=dlim1
	get_data,rbspx+'_efw_64_spec2',data=spec2,dlimits=dlim2
	get_data,rbspx+'_efw_64_spec3',data=spec3,dlimits=dlim3
	get_data,rbspx+'_efw_64_spec4',data=spec4,dlimits=dlim4
	get_data,rbspx+'_efw_64_spec5',data=spec5,dlimits=dlim5
	get_data,rbspx+'_efw_64_spec6',data=spec6,dlimits=dlim6

	chn0 = strlowcase(dlim0.data_att.channel)
	chn1 = strlowcase(dlim1.data_att.channel)
	chn2 = strlowcase(dlim2.data_att.channel)
	chn3 = strlowcase(dlim3.data_att.channel)
	chn4 = strlowcase(dlim4.data_att.channel)
	chn5 = strlowcase(dlim5.data_att.channel)
	chn6 = strlowcase(dlim6.data_att.channel)

	;Get official times at spec cadence
	datatimes = spec0.x
	epoch_spec = tplot_time_to_epoch(datatimes,/epoch16)


	;Get the time structure for the flag values. These are not necessarily at the cadence
	;of physical data.
;	epoch_flag_times,date,5,epoch_flags,timevals
;	epoch_flag_times,date,5,epoch_flags,datatimes












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
	tinterpol_mxn,rbspx+'_state_pos_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbspx+'_state_vel_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbspx+'_spinaxis_direction_gse',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbspx+'_state_mlt',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbspx+'_state_mlat',times_ephem,/overwrite,/spline
	tinterpol_mxn,rbspx+'_state_lshell',times_ephem,/overwrite,/spline

	get_data,rbspx+'_state_pos_gse',data=pos_gse
	get_data,rbspx+'_state_vel_gse',data=vel_gse
	get_data,rbspx+'_spinaxis_direction_gse',data=wgse
	get_data,rbspx+'_state_mlt',data=mlt
	get_data,rbspx+'_state_mlat',data=mlat
	get_data,rbspx+'_state_lshell',data=lshell
















	;Get all the flag values
;	flag_str = rbsp_efw_get_flag_values(sc,timevals,density_min=dmin,boom_pair=bp,_extra=extra)
	flag_str = rbsp_efw_get_flag_values(sc,times_ephem,density_min=dmin,boom_pair=bp,_extra=extra)
	flag_arr = flag_str.flag_arr
	bias_sweep_flag = flag_str.bias_sweep_flag
	ab_flag = flag_str.ab_flag
	charging_flag = flag_str.charging_flag






	;--------------------------------------------------
  ;Get burst times
	;This is a bit complicated for spinperiod data b/c the short
	;B2 snippets can be less than the spinperiod.
	;So, I'm padding the B2 times by +/- a half spinperiod so that they don't
	;disappear upon interpolation to the spinperiod data.
  ;--------------------------------------------------

	b1_flag = intarr(n_elements(times_ephem))
	b2_flag = b1_flag


	;get burst 1 and burst 2 times
	b1t = rbsp_get_burst_times_rates_list(sc)
	b2t = rbsp_get_burst2_times_list(sc)



  ;Pad B2 by +/- half spinperiod.
	b2t.startb2 -= 6.
  	b2t.endb2   += 6.

  for q=0,n_elements(b1t.startb1)-1 do begin
    goodtimes = where((times_ephem ge b1t.startb1[q]) and (times_ephem le b1t.endb1[q]))
    if goodtimes[0] ne -1 then b1_flag[goodtimes] = b1t.samplerate[q]
  endfor
  for q=0,n_elements(b2t.startb2[*,0])-1 do begin $
    goodtimes = where((times_ephem ge b2t.startb2[q]) and (times_ephem le b2t.endb2[q])) & $
    if goodtimes[0] ne -1 then b2_flag[goodtimes] = 1
  endfor





	;Rename the skeleton file
	filename = 'rbsp'+sc+'_efw-l2_spec_'+strjoin(strsplit(date,'-',/extract))+'_v'+vstr+'.cdf'
	file_copy,source_file,folder+filename,/overwrite

	;Open the new skeleton file
	cdfid = cdf_open(folder+filename)
	cdf_control, cdfid, get_var_info=info, variable='epoch'





	;Final list of variables to NOT delete
	varsave_general = ['epoch','epoch_spec',$
					   'spec64_'+chn0,'spec64_'+chn1,'spec64_'+chn2,$
					   'spec64_'+chn3,'spec64_'+chn4,'spec64_'+chn5,$
					   'spec64_'+chn6,$
					   'burst1_avail','burst2_avail',$
					   'flags_all','global_flag',$
					   'mlt','mlat','lshell',$
					   'position_gse','velocity_gse',$
					   'spinaxis_gse']







	;Get list of all the variable names in the CDF file.
	inq = cdf_inquire(cdfid)
	CDFvarnames = ''
	for varNum = 0, inq.nzvars-1 do begin $
		stmp = cdf_varinq(cdfid,varnum,/zvariable) & $
		if stmp.recvar eq 'VARY' then CDFvarnames = [CDFvarnames,stmp.name]
	endfor
	CDFvarnames = CDFvarnames[1:n_elements(CDFvarnames)-1]

stop


	cdf_varput,cdfid,'epoch',epoch_ephem
	cdf_varput,cdfid,'epoch_spec',epoch_spec

;	cdf_varput,cdfid,'efw_flags_all',transpose(flag_arr)
	cdf_varput,cdfid,'burst1_avail',b1_flag
	cdf_varput,cdfid,'burst2_avail',b2_flag

	cdf_varput,cdfid,'flags_all',transpose(flag_arr)
  	cdf_varput,cdfid,'global_flag',reform(flag_arr[*,0])

	cdf_varput,cdfid,'mlt',transpose(mlt.y)
	cdf_varput,cdfid,'mlat',transpose(mlat.y)
	cdf_varput,cdfid,'lshell',transpose(lshell.y)
	cdf_varput,cdfid,'position_gse',transpose(pos_gse.y)
	cdf_varput,cdfid,'velocity_gse',transpose(vel_gse.y)
	cdf_varput,cdfid,'spinaxis_gse',transpose(wgse.y)






	if is_struct(spec0) then cdf_varput,cdfid,'spec64_'+chn0,transpose(spec0.y)
	if is_struct(spec1) then cdf_varput,cdfid,'spec64_'+chn1,transpose(spec1.y)
	if is_struct(spec2) then cdf_varput,cdfid,'spec64_'+chn2,transpose(spec2.y)
	if is_struct(spec3) then cdf_varput,cdfid,'spec64_'+chn3,transpose(spec3.y)
	if is_struct(spec4) then cdf_varput,cdfid,'spec64_'+chn4,transpose(spec4.y)
	if is_struct(spec5) then cdf_varput,cdfid,'spec64_'+chn5,transpose(spec5.y)
	if is_struct(spec6) then cdf_varput,cdfid,'spec64_'+chn6,transpose(spec6.y)



	;Now delete unused CDF variables
	for qq=0,n_elements(CDFvarnames)-1 do begin $
		tstt = array_contains(varsave_general,CDFvarnames[qq]) & $
		if not tstt then print,'Deleting var:  ', CDFvarnames[qq] & $
		if not tstt then cdf_vardelete,cdfid,CDFvarnames[qq]
	endfor



stop

	cdf_close, cdfid



	store_data,tnames(),/delete



end
