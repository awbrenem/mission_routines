;+
; NAME:	RBSP_LOAD_EFW_ESVY_MGSE
;
; SYNTAX:
;   rbsp_load_efw_esvy_mgse,probe='a'
;   rbsp_load_efw_esvy_mgse,probe='a',/no_spice_load
;
; PURPOSE:	Loads EFW ESVY data and despins using SPICE via rbsp_uvw_to_mgse.pro
;
;			The MGSE coordinate system is defined:
;				Y_MGSE=-W_SC(GSE) x Z_GSE
;				Z_MGSE=W_SC(GSE) x Y_MGSE
;				X_MGSE=Y_MGSE x Z_MGSE
;			where W_SC(GSE) is the spin axis direction in GSE.
;
;			This is equivalent to the GSE coordinate system if the spin axis
;			lies along the X_GSE direction.
;
; KEYWORDS:
;	probe = 'a' or 'b'  NOTE: single spacecraft only, does not accept ['a b']
;		NOTE: defaults to probe='a'
;	/no_spice_load - skip loading/unloading of SPICE kernels
;		NOTE: This assumes spice kernels have been manually loaded using:
;			rbsp_load_spice_predict ; (optional)
;			rbsp_load_spice_kernels ; (required)
;	/debug - prints debugging info
;	/qa - load the QA test file instead of standard L1 file
; bad_probe -> integer indicating a bad probe.
;
;
; NOTES:
;
; HISTORY:
;	1. Created Nov 2012 - Kris Kersten, kris.kersten@gmail.com
;
; VERSION:
;   $LastChangedBy: aaronbreneman $
;   $LastChangedDate: 2017-04-05 12:09:53 -0500 (Wed, 05 Apr 2017) $
;   $LastChangedRevision: 23110 $
;   $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/rbsp_load_efw_esvy_mgse.pro $
;
;-


pro rbsp_load_efw_esvy_mgse,probe=probe,$
	debug=debug,qa=qa,bad_probe=bad_probe,$
	load_from_wake_file=load_from_wake_file,$
	_extra=extra


	etype='esvy'

	date = timerange()
	date = strmid(time_string(date[0]),0,10)

	; check probe keyword
	if ~keyword_set(probe) then begin
		message,"Probe keyword not set. Using default probe='a'.",/continue
		probe='a'
	endif else begin
		probe=strlowcase(probe) ; this turns any data type to a string
		if probe ne 'a' and probe ne 'b' then begin
			message,"Invalid probe keyword. Using default probe='a'.",/continue
			probe='a'
		endif
	endelse


	;Two options for loading Esvy data
	;1) load the CDF files Sheng created that have the wake effects removed.
	;2) create the Esvy products using rbsp_load_efw_waveform, etc.

	;--Method 1 - load Sheng's CDF file (only good for boom pairs 12, 34)

	;***************
	;***************
	;***************
	;***************
	;THIS NEEDS TO BE FIXED

	wakecdf = 1
;	if tag_exist(extra,'sheng_cdf') then begin
;		if extra.sheng_cdf then wakecdf = 1 else wakecdf = 0
;	endif else wakecdf = 0
	;***************
	;***************
	;***************
	;***************


	if wakecdf then begin


		;Possibly load the wake flag
		tr = timerange()
		if ~tdexists('rbsp'+probe+'_eu_wake_flag',tr[0],tr[1]) then $
			rbsp_load_wake_effect_cdf_file,probe


		get_data,'rbsp'+probe+'_eu_fixed',data=eu
		get_data,'rbsp'+probe+'_ev_fixed',data=ev
		ew = eu.y
		ew[*] = 0.

		store_data,'rbsp'+probe+'_efw_'+etype,eu.x,[[eu.y],[ev.y],[ew]]


		rbsp_uvw_to_mgse,probe,'rbsp'+probe+'_efw_'+etype,debug=debug,_extra=extra


	;--Method 2
	endif else begin


		; force reload of esvy in uvw coordinates without cleaning
		if ~keyword_set(bad_probe) then begin
			rbsp_load_efw_waveform, probe=probe, datatype=etype, type='cal',coord='uvw',/noclean
		endif else rbsp_efw_create_esvy_uvw_from_vsvy,date,probe,bad_probe


		rbsp_load_spice_kernels,_extra=extra

		rbsp_uvw_to_mgse,probe,'rbsp'+probe+'_efw_'+etype,debug=debug,_extra=extra

		rbsp_load_spice_kernels,/unload,_extra=extra


	endelse




end
