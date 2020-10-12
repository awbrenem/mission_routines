;+
;PROGRAM:   wave_quality.pro
;  STEREO TDS or LRSBurst Low Frequency Wave Finder (uses waves in SW coord)
;  version 2.0.0, 20110809 (remember to change version ID in startup output)
;
;ARGUMENTS:
;       TAG -> STRING, prepended to output filename
;
;KEYWORDS:
;       SC        ->  SC='a' or SC='b'
;       HQ        ->  HQ='h' or HQ='q', to designate honesty or quality stream
;       STARTTIME ->  starting time (format: '2007-04-08 00:00:00.000')
;		ENDTIME   ->  ending time
;		nopatdsgain -> don't use the preamp gain boost for low freqs. When this is selected
;						the program also won't remove the three lowest freqs. WARNING: WITHOUT THIS
;						CORRECTION YOU CAN'T TRUST ANYTHING < 200 HZ.
;		which_calib = '?'      ;defaults to '1'
;			'1' = antenna half-lengths
;			'2' = Zaslavsky 2011 calibrations
;			'3' = Bale 2008 calibrations (with 1 for pseudo dipole)
;		ROTSC -> rotate to sc coordinate system before calculating amp, freq, quality, etc. Note
;				that the pseudo-dipole channel will remain in SWAVES coordinates.
;
;
;RETURNS:
;		WRITES OUTPUT FILES:
;			tag+'STEREO_TDS_STB_q_waveclass.txt' --> for TDS quality data on STB
;
;
;NOTES:
;       typing 'q' at any point during execution will end program and output results to save files
;
;		Because of the short timespan of the TDS captures I use wavelet transforms to find
;		the freq at peak power. FFTs just don't do too well.
;		There will be a min detectable freq based on the length of the TDS capture.
;		Need to take this into account when analyzing bursts with different sample lengths.
;		Ex. for the 0.131s captures min freq is ~32.71 Hz.
;
;		quality=0 indicated cal_state, thrust_state
;
;
;---------------------------------------------------------------------------------------------------
;		Notes on gain corrections:
;			Method1: using gain correction
;				-The timeseries should be     err=err+tm_get_item_r4(sid,"Time_Series_Raw",raw_temp,npoints,n_items)
;				-Apply the preamp gain correction (patdsgain.pro). This primarily affects freqs at f<200 Hz which get
;					dramatically boosted. The curve flattens out at f>200 Hz.
;				-Apply the A/D conversion factor to mV/m  (2.5/65.5360)
;				-Remove the lowest three freq bins which we can't trust. (~anything less than 20 Hz)
;
;
;			Method2: ignoring gain correction
;				-The timeseries should be     err=err+tm_get_item_r4(sid,"Time_Series_mV",raw_temp,npoints,n_items)
;				 	This timeseries includes the A/D conversion and also has a simplified preamp gain correction
;					which is equal to the actual PA gain correction at f>200 Hz where the curve is flat.
;				-The power at f<200 Hz should be less than actual value and therefore there is no need to remove
;				 the lowest 3 freq bins
;---------------------------------------------------------------------------------------------------
;
;
;
;INCLUDED MODULES:
;       swtds_lff
;
;LIBS USED:
;       (none)
;
;DEPENDENCIES
;       *REQUIRES TMlib SIDL <http://homepage.mac.com/swaves>
;       sw_errordump.pro (display error dumps)
;       make_scet.pro    (conversion from yyyymmdd,hhmmss to UR8)
;-
;
;
;-KEY for APM and Burst data:
;   nchannel=0 - Ex APM
;   nchannel=1 - Ey APM
;   nchannel=2 - Ez APM
;   nchannel=3 - Ex LRS
;   nchannel=4 - Ex LRS
;   nchannel=5 - Ex LRS
;   nchannel=6 - Ex-Ey LRS
;   nchannel=7 - Ex-Ez LRS
;
;-KEY for TDS data:
;   nchannel=0 - empty
;   nchannel=1 - Ex
;   nchannel=2 - Ey
;   nchannel=3 - Ez
;   nchannel=4 - Ex-Ey
;
;
;ANTENNA CALIBRATIONS (Zaslavsky, 2011 for all antennas except for Y-monopole)
;					  (Personal email from Zaslavsky for Y-monopole)
;
;
;Note that the y-monopole can be from ch1 or ch2, whereas x is always ch1 and
;z is always ch2
;
;-----------------------------------
;STA          Gamma*Leff (m)
;X         1.40 +- 0.079
;Y  Ch1    1.83 +- 0.09
;Y  CH2    1.91 +- 0.17
;Z  	   1.19 +- 0.15
;XY 	   2.02 +- 0.074
;YZ		   2.06 +- 0.12
;
;STB
;X		   1.51 +- 0.073
;Y  CH1    1.86 +- 0.04
;Y  CH2    1.73 +- 0.11
;Z		   1.08 +- 0.10
;XY		   2.04 +- 0.060
;YZ		   1.85 +- 0.15
;-----------------------------------
;
;
;
;
;CREATED BY:    Aaron Breneman - modified heavily from original program by Kris Kersten (thresh_lff.pro)
;
;HISTORY:
;		Created: 12/02/09 - AWB
;				 01/20/11 - AWB - many modifications
;				 04/04/11 - AWB - changed scet_ur8 call to scet_start_ur8. This is the actual start time
;									of the TDS capture, whereas scet_ur8 is when the burst capture is sent
;									to memory, or something.
;								  Gain correction removed. We decided that it is better to use the raw signals
;								  b/c nobody really understands the gain correction yet.
;				04/11/11 - AWB -  Updated effective antenna lengths
;				07/11/11 - AWB -  Changed format so that all data for current event is loaded, then things like
;									amp, freq, quality, etc are determined. This allows the rotation of waveforms
;									to other coord systems.
;								  Added "nopatdsgain" and "rotsc" keywords
;				08/09/11 - AWB -  Most of the calculation of various wave quantities is done outside
;								  of this program in get_wave_attributes.pro. Also, the wavelet spec
;								  is no longer used to calculate these quantities. Instead the FFT is used,
;								  speeding up the program
;
;#############
    ;APPLY CORRECT CALIBRATIONS TO THE RAW_TEMP FOR THE LRSBURST DATA.
;#############
;-


;wave_quality,sc='a',hq='q',starttime='2007-04-08 00:00:00.000',endtime='2007-04-08 05:00:00.000',which_calib='3'
;wave_quality,sc='a',hq='q',starttime='2011-06-01 00:00:00.000',endtime='2012-01-01 00:00:00.000',which_calib='2',/rotsc

;.compile /Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/stereo_wave_id/IDL> wave_quality,sc='a',hq='q',starttime='2007-04-08 00:00:00.000',endtime='2007-04-09 00:00:00.000',which_calib='2'


pro wave_quality,tag,sc=sc,hq=hq,starttime=starttime,endtime=endtime,$
	nopatdsgain=nopatdsgain,which_calib=which_calib,rotsc=rotsc

if not keyword_set(tag) then tag = '~/Desktop/'
if not keyword_set(which_calib) then which_calib = '1'

if which_calib eq '1' then begin
	;Half-length calibrations (maybe this is the best for the VLF range?)
	version='wave_quality.pro v2.0.0, 20110809 - antenna half-lengths used for calibrations'
	exel = 3.
;	eyel_ch1 = 3.
;	eyel_ch2 = 3.
	eyel = 3.
	ezel = 3.
	exeyel = 3.
	eyezel = 3.
endif
if which_calib eq '2' then begin
	version='wave_quality.pro v2.0.0, 20110809 - Zaslavsky 2011 calibrations used'
	if sc eq 'a' then begin
		exel = 1.40
;		eyel_ch1 = 1.83
;		eyel_ch2 = 1.91

		eyel = 1.91
		ezel = 1.19
		exeyel = 2.02
		eyezel = 2.06
	endif else begin
		exel = 1.51
;		eyel_ch1 = 1.86
;		eyel_ch2 = 1.73
		eyel = 1.73
		ezel = 1.08
		exeyel = 2.04
		eyezel = 1.85
	endelse
endif
if which_calib eq '3' then begin
	version='wave_quality.pro v2.0.0, 20110809 - Bale08 calibrations'
	exel = 1.17
;	eyel_ch1 = 1.44
;	eyel_ch2 = 1.44
	eyel = 1.44
	ezel = 0.97
	exeyel = 1.
	eyezel = 1.
endif



;----GLOBAL VARS----
err=0
sid=0 ; StreamID
mission="STEREO"
inst="SWAVES"
read_spcr=''
spcr=''
read_ev_type=''
ev_type=''
yyyymmdd=0D
hhmmss=0D
UR8start=0d
scetur8=0D
SCETstring=''
nchannel=0 ; what channel is TMlib returning?
read_quit=''
lowcount=0L
midcount=0L
highcount=0L
dustcount=0L
max_amplitude=0.
oldevid=0L
newevid=0L
count_a = 0L ;total number of regular antenna TDS captures
count_p = 0L ;  '      '    ' pseudo dipole    '      '
counter = 0.
qtmp = 0.

nullstr = '                       '

print,""
print,version
print,""


; check input keywords
if keyword_set(sc) then begin
    if strlowcase(sc) eq 'a' then read_spcr='a'
    if strlowcase(sc) eq 'b' then read_spcr='b'
endif
if keyword_set(hq) then begin
    if strlowcase(hq) eq 'h' then read_ev_type='h'
    if strlowcase(hq) eq 'q' then read_ev_type='q'
endif


;---------------------------------------------------
;Direct stream selection method - get the basic info
;---------------------------------------------------

while read_spcr ne 'a' and read_spcr ne 'b' do read,read_spcr,prompt="Enter: a for Ahead, b for Behind:  "
if read_spcr eq 'a' then spcr="STEREO_A" else spcr="STEREO_B"

while read_ev_type ne 'h' and read_ev_type ne 'q' and read_ev_type ne 'b' do read,read_ev_type,$
	prompt="Enter: h for Honesty, q for Quality, b for LRSBurst/APM:  "
if read_ev_type eq 'h' then ev_type="TDS_a"
if read_ev_type eq 'q' then ev_type="TDS_b"
if read_ev_type eq 'b' then ev_type="LRSBurst"

if ev_type eq 'TDS_a' or ev_type eq 'TDS_b' then ev_type2 = 'TDS'
if ev_type eq 'LRSBurst' then ev_type2 = 'LRSBurst'




;-----------------------
;-try to open the stream
;-----------------------

print,"Selecting Domain:  "+mission+"/"+spcr+"/"+inst+"/"+ev_type
stereo_open_stream,sid,mission,spcr,inst,ev_type,starttime,nchannel,endtime=endtime,$
                   evid,ev_err,UR8start,UR8end

print,evid
print,UR8start
print,UR8end







if not keyword_set(rotsc) then outfile1 = tag + mission + '_' + ev_type2 + '_' + 'ST' + strupcase(read_spcr) + '_' + read_ev_type + '_waveclass_SWCOORD.txt'
if keyword_set(rotsc) then outfile1 = tag + mission + '_' + ev_type2 + '_' + 'ST' + strupcase(read_spcr) + '_' + read_ev_type + '_waveclass_SCCOORD.txt'


;format for printed and regular (screen) output
format1 =    '(a24,5x,a3,2x,i5,2x,i5,2x,i3,2x,i5,2x,f7.2,2x,i3,2x,i5,2x,f4.2,2x,i3,2x,i9,2x,a3,2x,i5)'


openw,outlun1,outfile1,/get_lun

printf,outlun1,'Run started: ' + systime()
printf,outlun1,'TDS waveform capture data for: ' + ev_type2
printf,outlun1,'Spacecraft: ' + read_spcr
if read_ev_type eq 'h' or read_ev_type eq 'q' then printf,outlun1,'Quality or honesty: ' + read_ev_type else printf,outlun1,'Quality or honesty: --'
printf,outlun1,'Instrument: ' + inst
printf,outlun1,'Mission: ' + mission
printf,outlun1,'Version: ' + version
if ev_type2 eq 'TDS_STA' or ev_type2 eq 'TDS_STB' then printf,outlun1,'Channels = [Ex,Ey,Ez,Ex-Ey]'
if ev_type2 eq 'LRSBurst_APM' then printf,outlun1,'Channels = [Ex_APM,Ey_APM,Ez_APM,Ex_LRSBurst,Ey_LRSBurst,Ez_LRSBurst,Ex-Ey_LRSBurst,Ex-Ez_LRSBurst]'
printf,outlun1,'Freq = Frequency in Hz'
printf,outlun1,'Qual = arbitrary measure of quality'
printf,outlun1,'Amp = Wave amplitude in mV/m (zero-peak)'
printf,outlun1,'Bw = bandwidth of peak power in Hz'
printf,outlun1,'df/f = normalized bandwidth'
printf,outlun1,'drate = data rate'
printf,outlun1,'npts = number of samples in burst'
printf,outlun1,'blen = burst length in sec'
printf,outlun1,'maxf = data removed above this frequency '
printf,outlun1,'ID = official event ID number'
printf,outlun1,'ts = channel TDS triggered off of'
printf,outlun1,'TDSQ = quality of event as defined onboard the sc'
printf,outlun1,'format: ' + '(a24,5x,a3,2x,i5,2x,i5,2x,i3,2x,i5,2x,f7.2,2x,i3,2x,i5,2x,f4.2,2x,i3,2x,i9,2x,a3,2x,i5)'

printf,outlun1,''
print,'Opened output files: '+outfile1
print,''


printf,outlun1,'        Time                  Src   Freq  Qual  Amp    BW     df/f  drate npts   blen  maxf    ID       ts   TDSQ'

eventquality = replicate(!values.f_nan,5)
eventsource = ['--','--','--','--','--']
eventfreq = replicate(!values.f_nan,5)
event_max_amp = replicate(!values.f_nan,5)
eventbandwidth = replicate(!values.f_nan,5)
;eventbv = replicate(!values.f_nan,5)
event_deltaf_to_f = replicate(!values.f_nan,5)
event_drate = ['--','--','--','--','--'] ;data rate
event_npoints = replicate(!values.f_nan,5)    ;number of samples
event_stime = replicate(!values.f_nan,5)      ;burst sample time
event_filter = ['--','--','--','--','--']      ;Filter
event_id = replicate(!values.f_nan,5)         ;EventID
event_triggersource = ['--','--','--','--','--']
event_TDS_Q = replicate(!values.f_nan,5)         ;Event quality as defined onboard the sc
event_deltat = replicate(!values.f_nan,5)


scetur8 = UR8start


;----------------------------
;----MAIN PROGRAM LOOP----
;----------------------------

while ev_err eq 0 and read_quit ne 'q' and scetur8 le UR8end do begin

	currentevid = evid

	err=err+tm_get_item_i4(sid,"Nsamples",npoints,1l,n_items)
	timeseries_tmp = fltarr(npoints,5)


	;Get all the data for current event
	while (evid eq currentevid) and (ev_err eq 0) do begin

		thrust=0.
		err=0
		err=err+tm_get_item_char(sid,"Source",source,10,n_items)
		err=err+tm_get_item_i4(sid,"Nsamples",npoints,1l,n_items)
		err=err+tm_get_item_i4(sid,"Channel",nchannel,1l,n_items)
		if ~keyword_set(nopatdsgain) then $
    err=err+tm_get_item_r4(sid,"Time_Series_Raw",raw_temp,npoints,n_items)   ;use when applying patdsgain and removing lowest 3 freqs
	  if keyword_set(nopatdsgain) then $
    err=err+tm_get_item_r4(sid,"Time_Series_mV",raw_temp,npoints,n_items)   ;use when not applying patdsgain or removing 3 lowest freqs
		err=err+tm_get_item_r4(sid,"Filter",ifilt,1L,n_items) ;needed for patdsgain.pro
		err=err+tm_get_item_char(sid,"Filter",ifilt2,10,n_items)
		err=err+tm_get_item_r8(sid,"scet_start_ur8",scetur8,1l,n_items)
		err=err+tm_get_item_r4(sid,"Sample_Period_Seconds",deltat,1l,n_items)
		err=err+tm_get_item_r4(sid,"Sample_Speed_SPS",nsps,1l,n_items)
		err=err+tm_get_item_char(sid,"Sample_Speed",ssps,10,n_items)
		err=err+tm_get_item_char(sid,"Trigger_Source",ts,10,n_items)
		err=err+tm_get_item_i4(sid,"Thruster_State",thrust,1l,n_items)
		err=err+tm_get_item_i4(sid,"Event_Quality",tds_quality,1l,n_items)
		err=err+tm_get_item_i4(sid,"CAL_State",cal_state,1l,n_items)
		err=err+tm_ur8_to_string_absolute(scetur8,scetstr)

		evlength=deltat*npoints


		if source ne 'Ex-Ey' then source = strmid(source,1,1) else source = 'x-y'
		if ts ne 'Ex-Ey' then ts = strmid(ts,1,1) else ts = 'x-y'
		ssps = strmid(ssps,0,strpos(ssps,'k'))
		ifilt2 = strmid(ifilt2,0,strpos(ifilt2,'k'))


		timeseries_tmp[*,nchannel] = raw_temp
		eventsource[nchannel] = source
		event_drate[nchannel] = ssps
		event_npoints[nchannel] = npoints
		event_stime[nchannel] = evlength
		event_filter[nchannel] = ifilt2
		event_id[nchannel] = currentevid
		event_triggersource[nchannel] = ts
		event_TDS_Q[nchannel] = tds_quality
		event_deltat[nchannel] = deltat




		;turn month string into number
		monthtmp = strmid(scetstr,3,3)

    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    gooo = where(months eq monthtmp)
    monthtmp = gooo[0] + 1
    if monthtmp lt 10 then monthtmp = '0' + strtrim(monthtmp,2) else monthtmp = strtrim(monthtmp,2)



		scetstr = strmid(scetstr,7,4) + '-' + monthtmp + '-' + strmid(scetstr,0,2) + '/' + strmid(scetstr,12,12)

		if thrust eq 2 then begin
			thrust=0
			if spcr eq 'STEREO_A' and ((9112.521423612d le scetur8 and scetur8 le 9112.522754640d) $
				or (9148.806319450d le scetur8 and scetur8 le 9148.808761578d) $
				or (9190.813946766d le scetur8 and scetur8 le 9190.816342595d) $
				or (9232.793599547d le scetur8 and scetur8 le 9232.795601863d) $
				or (9273.606944454d le scetur8 and scetur8 le 9273.608750005d) $
				or (9316.545648153d le scetur8 and scetur8 le 9316.547581025d) $
				or (9358.948125007d le scetur8 and scetur8 le 9358.949803248d) $
				or (9398.541990741d le scetur8 and scetur8 le 9398.543344918d) $
				or (9442.729456019d le scetur8 and scetur8 le 9442.731192130d) $
				or (9484.583553246d le scetur8 and scetur8 le 9484.585254632d) $
				or (9532.573576394d le scetur8 and scetur8 le 9532.575057874d) $
				or (9580.687766208d le scetur8 and scetur8 le 9580.689537038d)) then thrust=1
			if spcr eq 'STEREO_B' and ((9183.584560188d le scetur8 and scetur8 le 9183.586203713d) $
				or (9225.584814825d le scetur8 and scetur8 le 9225.586909725d) $
				or (9275.836168985d le scetur8 and scetur8 le 9275.837997695d) $
				or (9336.791990752d le scetur8 and scetur8 le 9336.793530096d) $
				or (9398.750046305d le scetur8 and scetur8 le 9398.752268519d) $
				or (9462.667372688d le scetur8 and scetur8 le 9462.669039363d) $
				or (9518.646562505d le scetur8 and scetur8 le 9518.648009271d) $
				or (9573.893067141d le scetur8 and scetur8 le 9573.894583335d)) then thrust=1
		endif


		ev_err=tm_find_event(sid)
		err=err+tm_get_item_i4(sid,"EventID",evid,1l,n_items)
		undefine,raw_temp

	endwhile  ;for evid eq currentevid


	;Get rid of first channel, which is blank
	timeseries_tmp = timeseries_tmp[*,1:4]
	eventsource = eventsource[1:4]
	event_drate = event_drate[1:4]
	event_npoints = event_npoints[1:4]
	event_stime = event_stime[1:4]
	event_filter = event_filter[1:4]
	event_id = event_id[1:4]
	event_triggersource = event_triggersource[1:4]
	event_TDS_Q = event_TDS_Q[1:4]
	event_deltat = event_deltat[1:4]


	;-----------------------------------------------------------------------------------------
	;APPLY GAIN CORRECTION AND EFFECTIVE ANTENNA LENGTHS
	;-----------------------------------------------------------------------------------------


	for chn=0,3 do begin


		if eventsource[chn] ne '--' then begin

		if not keyword_set(nopatdsgain) then begin

				;-calculate gain corrected mV from raw TM
				;-set up gain correction array
			tempf=1.0/(event_deltat[chn]*event_npoints[chn])
			tds_corr=complexarr(event_npoints[chn])
			mvc_temp=fltarr(event_npoints[chn])
			tds_corr[0]=complex(1.0,0.0)
			tds_corr[event_npoints[chn]/2]=tds_corr[0]
			ipa=1
			for ccount=1,event_npoints[chn]/2-1 do begin
			  f=tempf*ccount
			  patdsgain,f,ipa,ifilt,cgaintds,gaintds,phasetds
			  tds_corr[ccount]=cgaintds
			  tds_corr[event_npoints[chn]-ccount]=conj(cgaintds)
			endfor


			;The 2.5/65.5360 is the A/D conversion factor and is only for the raw (not mV) data
			fft_temp=(2.5/65.5360)*fft(float(timeseries_tmp[0:event_npoints[chn]-1,chn]))/tds_corr[0:event_npoints[chn]-1]

			;-zero out three lowest frequency bands
			fft_temp[0:3]=0
			fft_temp[event_npoints[chn]-(3+1):event_npoints[chn]-1]=0
			mvc_temp[0:event_npoints[chn]-1]=float(fft(fft_temp[0:event_npoints[chn]-1],1))

			timeseries_tmp[*,chn] = mvc_temp[0:event_npoints[chn]-1]     ;select if using the gain correction and patdsgain

		endif


		if eventsource[chn] eq 'Ex' then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/exel
	;	if eventsource[chn] eq 'Ey' and event_nchannel[chn] eq 1 then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/eyel_ch1
	;	if eventsource[chn] eq 'Ey' and event_nchannel[chn] eq 2 then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/eyel_ch2

		if eventsource[chn] eq 'Ey' then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/eyel
		if eventsource[chn] eq 'Ez' then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/ezel
		if eventsource[chn] eq 'Ex-Ey' then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/exeyel
		if eventsource[chn] eq 'Ey-Ez' then timeseries_tmp[*,chn] = timeseries_tmp[*,chn]/eyezel


		endif
	endfor


	;------------------------------------
	;ROTATE DATA TO SPACECRAFT COORD SYSTEM
	;------------------------------------

	if keyword_set(rotsc) then begin

		for rotcount=0,event_npoints[1]-1 do begin
;			tempvec = swrot_Vant2Esc([timeseries_tmp[rotcount,1],timeseries_tmp[rotcount,2],timeseries_tmp[rotcount,3]])

;Changed by Zac on 2017-06-27
			tempvec = swrot_Vant2Esc([timeseries_tmp[rotcount,0],timeseries_tmp[rotcount,1],timeseries_tmp[rotcount,2]])
			timeseries_tmp[rotcount,1] = tempvec[0]
			timeseries_tmp[rotcount,2] = tempvec[1]
			timeseries_tmp[rotcount,3] = tempvec[2]
		  endfor

	endif

;-----------------------------------------------------------------------------------
;Now that we have all the data for every channel read in let's do some calculations
;-----------------------------------------------------------------------------------

	for chn=0,3 do begin

		if eventsource[chn] ne '--' then begin

			;--------------------------------------
			;CALCULATE VARIOUS WAVE QUALITIES
			;--------------------------------------


			maxtime = event_deltat[chn]*event_npoints[chn]
			time = evlength*indgen(event_npoints[chn])/(event_npoints[chn]-1)

			wa = get_wave_attributes(time,timeseries_tmp[*,chn],wssn='neither')


			;-------------------------------------------
			;DETERMINE WAVE QUALITY FROM AUTOCORRELATION
			;-------------------------------------------

			if ev_type ne 'LRSBURST' THEN BEGIN

				;all possible lags (slow)
				lag2 = 2d*indgen(event_npoints[chn]) - (event_npoints[chn] - 1d)
				;subset of all possible lags
				lag2 = lag2[event_npoints[chn]/2.-(event_npoints[chn]/3.28):event_npoints[chn]/2.+(event_npoints[chn]/3.28)]
				ac = a_correlate(timeseries_tmp[*,chn]/max(timeseries_tmp[*,chn],/nan),lag2,/double)

				;find ratio of second largest peak to largest peak in autocorrelation.
				;This method is not susceptible to high frequency autocorrelation signals
				;which can create local small peaks
				half = n_elements(ac)/2.
				actmp = reverse(ac[0:half-1])
				jj = 0l
				maxim = 0l

				gt0 = where(actmp gt 0.)

				if gt0[0] ne -1 then begin
				;find first zero crossing of autocorrelation
				while (actmp[jj] ge 0.) and (jj lt n_elements(actmp)-2) do jj++
				;get back to positive autocorrelation values
				while (actmp[jj] lt maxim) and (jj lt n_elements(actmp)-2) do jj++
				;find max positive peak from rest of autocorrelation array
				endif

				maxim = max(actmp[jj:n_elements(actmp)-1])

				ratio_ac = abs(maxim/actmp[0])  ;ratio of max peak (value=1 for autocorrelation) to the second peak

				quality = 100*ratio_ac/wa.deltaf_to_f


			endif else quality = !values.f_nan


			;-----------------
			;TEST FOR DUST
			;-----------------
			;ratio_ac will be small
			;peak value on one side will be much larger than peak on other

			dust = 'no'

			goo1 = where(timeseries_tmp[*,chn] gt 0.)
			goo2 = where(timeseries_tmp[*,chn] lt 0.)

			if goo1[0] ne -1 then max1 = max(timeseries_tmp[goo1,chn]) else max1 = 0.
			if goo2[0] ne -1 then max2 = min(timeseries_tmp[goo2,chn]) else max2 = 0.

			if max1 ne 0. and max2 ne 0. then begin
				if abs(max1) ge abs(max2) then wr = abs(max1/max2) else wr = abs(max2/max1)
			endif else wr = 10. ;arbitrary number to classify as dust


			t1 = where(ac gt 0.)
			t2 = where(ac lt 0.)

			if t1[0] ne -1 and t2[0] ne -1 then begin
				max1 = max(ac(t1))
				max2 = min(ac(t2))
				if abs(max1) ge abs(max2) then ar = abs(max1/max2) else ar = abs(max2/max1)
			endif else begin
				ar = 1000.  ;arbitrary ratio
			endelse

			;four conditions that classify wave as dust
			if quality le 1. then dust = 'yes'
			if wr ge 4 and ar ge 2 then dust = 'yes'
			if ar ge 4 and wr ge 2 then dust = 'yes'
			if ar eq 1000 then dust = 'yes'


			;------------------------------------------------------------------------------
			;RESET WAVE QUALITY IF SIGNAL IS DUST OR FREQUENCY IS OUTSIDE OF ALLOWED RANGE
			;------------------------------------------------------------------------------

			if thrust eq 1 or cal_state eq 1 then quality=0.
			if keyword_set(qtemp) then if qtemp eq -1 then quality = -1
			if dust eq 'yes' then quality = -1


			;--------------------
			;SAVE VALUES TO ARRAY
			;--------------------


			eventquality[chn] = quality
			eventfreq[chn] = wa.freq_of_max_power
			event_max_amp[chn] = max(abs(timeseries_tmp[*,chn]))
			eventbandwidth[chn] = wa.bandwidth
			;eventbv[chn] = wa.bandwidth_stdev
			event_deltaf_to_f[chn] = wa.deltaf_to_f


		endif
	endfor


event_drate = floor(float(event_drate))
event_filter = floor(float(event_filter))


	;-------------
	;SAVE TO FILE
	;-------------

	if counter eq 10. then print,'PRESS "q" TO QUIT AND SAVE CURRENT RESULTS'


		printf,outlun1,format=format1,scetstr,eventsource[0],eventfreq[0],eventquality[0],event_max_amp[0],eventbandwidth[0],event_deltaf_to_f[0],event_drate[0],event_npoints[0],event_stime[0],event_filter[0],event_id[0],event_triggersource[0],event_TDS_Q[0]
		printf,outlun1,format=format1,nullstr,eventsource[1],eventfreq[1],eventquality[1],event_max_amp[1],eventbandwidth[1],event_deltaf_to_f[1],event_drate[1],event_npoints[1],event_stime[1],event_filter[1],event_id[1],event_triggersource[1],event_TDS_Q[1]
		printf,outlun1,format=format1,nullstr,eventsource[2],eventfreq[2],eventquality[2],event_max_amp[2],eventbandwidth[2],event_deltaf_to_f[2],event_drate[2],event_npoints[2],event_stime[2],event_filter[2],event_id[2],event_triggersource[2],event_TDS_Q[2]
		printf,outlun1,format=format1,nullstr,eventsource[3],eventfreq[3],eventquality[3],event_max_amp[3],eventbandwidth[3],event_deltaf_to_f[3],event_drate[3],event_npoints[3],event_stime[3],event_filter[3],event_id[3],event_triggersource[3],event_TDS_Q[3]

		print,'CAPTURE NUMBER ' + strtrim(count_a,2) + ' : ' + scetstr

		count_a+=1



		timeseries_tmp = fltarr(npoints,5)
		eventquality = replicate(!values.f_nan,5)
		eventsource = ['--','--','--','--','--']
		eventfreq = replicate(!values.f_nan,5)
		event_max_amp = replicate(!values.f_nan,5)
		eventbandwidth = replicate(!values.f_nan,5)
		event_deltaf_to_f = replicate(!values.f_nan,5)
		event_drate = ['--','--','--','--','--'] ;data rate
		event_npoints = replicate(!values.f_nan,5)    ;number of samples
		event_stime = replicate(!values.f_nan,5)      ;burst sample time
		event_filter = ['--','--','--','--','--']      ;Filter
		event_id = replicate(!values.f_nan,5)         ;EventID
		event_triggersource = ['--','--','--','--','--']
		event_TDS_Q = replicate(!values.f_nan,5)         ;Event quality as defined onboard the sc
		event_deltat = replicate(!values.f_nan,5)



		qtemp = 0.

		if counter eq 10. then counter = 0.
		counter = counter + 1

		read_quit=get_kbrd(0)

endwhile  ;no more TDS captures





;------------
;FINISH UP
;------------


printf,outlun1,'Run terminated: ' + systime()
print,""
close,outlun1,exit_status=fstatus1
print,'Wrote ' + strcompress(string(count_a),/remove_all) + ' events to file: ' + outfile1
print,'Output files closed with status: ',fstatus1 ;,fstatus2
print,'Final time saved: ' + scetstr
err = tm_close(StreamId)

sttm = strmid(starttime,0,4) + strmid(starttime,5,2) + strmid(starttime,8,2)
goo = TM_UR8_to_string_comparison(scetUR8,edtm)
edtm = strmid(edtm,0,4) + strmid(edtm,5,2) + strmid(edtm,8,2)

goo = strpos(outfile1,'.txt')
outfilet = strmid(outfile1,0,goo)
outfile2 = outfilet + '_' + sttm + '_to_' + edtm + '.txt'


file_move,outfile1,outfile2,/overwrite


return
end














;###
	;for dust impacts
;	if zerocrossings eq 4 then quality=-1
;	if laindex[0] ne -1 then quality=-1
;	if ratio ge 5 then quality = -1
;###


;plot the original and removed values
;	!p.multi=[0,0,2,0,0]
;	plot,abs(timeseries),xstyle=1,xrange=[0,2e4],title = 'channel ' + strtrim(nchannel,2) + ', time= ' + scetstr
;	if xright[0] ne -1 or xleft[0] ne -1 then begin
;		numb = n_elements(timeseries)-(2*n_elements(tmpright))
;		xxx = replicate(0,numb)
;		plot,[tmpleft,xxx,tmpright],xstyle=1,xrange=[0,2e4]
;	endif else plot,[0,0]
;	help,xright,xleft
;	if keyword_set(numb) then print,'numb = ',numb else print,'numb = none'

	;laindex = where(abs(timeseries) gt 100.)
	;tmpp = where(timeseries ge 0.)
	;maxpos = max(timeseries[tmpp])
	;tmpp = where(timeseries lt 0.)
	;maxneg = max(abs(timeseries[tmpp]))

	;ratio = maxneg/maxpos
;#######


;#################################################
;#################################################
;#################################################
;------------------------
;COUNT ZERO CROSSINGS
;------------------------

;do_zc = 'no'

;if do_zc eq 'yes' then begin
;
;stop
;	meanval=mean(timeseries)
;	peakval=max([abs(timeseries-meanval)])
;	if not keyword_set(thres) then threshold=0.25*peakval else threshold=thres*peakval
;
;	gtindex=where(timeseries gt meanval+threshold)
;	ltindex=where(timeseries lt meanval-threshold,complement=ltcomplement)
;	midindex=where(timeseries[ltcomplement] lt meanval+threshold)
;
;
;	tmpright = abs(timeseries[(npoints/2. + 0.05*n_elements(timeseries)):(n_elements(timeseries)-1)])
;	tmpleft =  abs(timeseries[0:(npoints/2. - 0.05*n_elements(timeseries))])
;
;	xright = where(tmpright ge threshold)
;	xleft = where(tmpleft ge threshold)
;
;	if xright[0] eq -1 and xleft[0] eq -1 then qtemp = -1
;
;
;
;	keyarray=intarr(size(timeseries,/n_elements))
;
;	if gtindex[0] ne -1 then keyarray[gtindex]=1
;	if ltcomplement[0] ne -1 and midindex[0] ne -1 then keyarray[ltcomplement[midindex]]=0
;	if ltindex[0] ne -1 then keyarray[ltindex]=-1
;
;	zerocrossings=0
;	zc_pairs=intarr(size(timeseries,/n_elements),2)
;	keyend=size(keyarray,/n_elements)-1
;	keycount=0
;
;	while keyarray[keycount] eq 0 do keycount+=1
;    lastgood=keyarray[keycount]
;	lastgoodindex=keycount
;	to = deltat*lastgoodindex
;    keycount+=1
;    while keycount le keyend do begin
;    	if lastgood*keyarray[keycount] eq -1 then begin
;    		zc_pairs[zerocrossings,*]=[lastgoodindex,keycount]
;    		zerocrossings+=1
;    		lastgood=keyarray[keycount]
;    		lastgoodindex=keycount
;    		tf = deltat*lastgoodindex
;    	endif else if keyarray[keycount] ne 0 then lastgoodindex=keycount
;    	keycount+=1
;    endwhile
;	zctime = tf - to ;time b/t first and last zero crossing
;
;	if zerocrossings eq 0 then quality=0 $
;	else begin
;		zc_pairs=zc_pairs[0:zerocrossings-1,*]
;		zc_index=(zc_pairs[*,0]+zc_pairs[*,1])/2
;		zc_index=[0,zc_index,npoints-1]
;		zc_spread=[zc_index[0:zerocrossings-1],npoints-1]-[0,zc_index[0:zerocrossings-1]]
;		if zerocrossings ge 3 then zc_spread=zc_spread[1:zerocrossings-1] ; ignore endpoints if we have enough ZCs
;		std_zc_spread=stddev(zc_spread)
;		if std_zc_spread eq 0 then std_zc_spread=0.1 ; prevent divide-by-zero
;		mean_zc_spread=mean(zc_spread)
;		quality=mean_zc_spread/std_zc_spread
;	endelse
;	frequency=float(zerocrossings/2)/evlength
;	frequency = float(zerocrossings/2)/zctime  ;more accurate frequencies

;	if frequency lt fmin or frequency gt fmax then quality=0.



;endif ;for zero crossings
;######################################################
;######################################################
;######################################################


;if keyword_set(startdate) then yyyymmdd=startdate else read,yyyymmdd,prompt="Enter date (yyyymmdd):  "
;if keyword_set(enddate) then yyyymmdd2 = enddate

;yyyymmdd = long(yyyymmdd)
;if keyword_set(enddate) then yyyymmdd2 = long(yyyymmdd2)

;if not keyword_set(starttime) then read,starttime,prompt="Enter start time (hhmmss): "
;if keyword_set(enddate) and not keyword_set(endtime) then read,endtime,prompt="Enter end time (hhmmss): "

;if strlen(starttime) eq 4 then starttime = starttime + '00'
;if strlen(starttime) eq 3 then starttime = '0' + starttime + '00'
;if strlen(starttime) eq 5 then starttime = '0' + starttime
;starttime = float(starttime)

;if keyword_set(endtime) then begin
;	if strlen(endtime) eq 4 then endtime = endtime + '00'
;	if strlen(endtime) eq 3 then endtime = '0' + endtime + '00'
;	if strlen(endtime) eq 5 then endtime = '0' + endtime
;	endtime = float(endtime)
;endif


;make_scet,yyyymmdd,starttime,UR8start
;if keyword_set(yyyymmdd2) then make_scet,yyyymmdd2,endtime,UR8end


;err=TM_UR8_to_String_Absolute(UR8start,SCETstring)
;if err ne 0 then begin
;    sw_errordump,sid,errorid="TIME CONVERSION"
;    return
;endif

;print,'start time conversion OK...'
;print,"  UR8:  ",UR8start
;print,"  STR:  ",SCETstring
;print,''


;-----------------------
;test time stuff
;-----------------------
