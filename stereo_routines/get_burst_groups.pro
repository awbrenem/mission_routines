;+
;*****************************************************************************************
;
;  FUNCTION :   get_burst_groups.pro
;  PURPOSE  :   Finds the burst wave groups for STEREO TDS based on the structure returned
;				from either get_TDS.pro (method1) or identify_waves.pro (method2). Probably
;				best to use method2 since identify_waves.pro does a better job accurately
;				identifying the waves
;
;  CALLED BY:   
;               
;
;  CALLS:
;               
;
;  REQUIRES:    
;               
;
;  INPUT:		sc -> 'STA' or 'STB'
;				tds -> structure from the program get_TDS.pro (method1) or identify_waves.pro (method2)
;					   that contains all of the TDS captures for a certain timerange. The only difference
;						b/t the two methods is that method2 already has the wave type identified. Therefore
;						this program won't try to identify the waves. Program automatically identifies
;						the correct method
;						
;				sepmax -> max time (minutes) allowed b/t events. Defaults to 10 min
;				freqmin ->
;				freqmax -> max frequency of TDS capture to be considered (Hz). Defaults to 500 Hz
;				bwmin ->
;				bwmax -> max bandwidth (Hz) value. Defaults to 70
;				dfmin ->
;				dfmax -> max value of df/f. Defaults to 5
;				ampmin ->
;				ampmax ->
;;;				bwstddevmin -> NO LONGER USED
;;;				bwstddevmax ->
;				qualmin -> min quality of burst to consider. Defaults to 100
;				qualmax -> 
;
;				densmin -> min "event density" (events/minute). Defaults to 1
;				noevents_min -> min number of events that a group must be comprised of. Defaults to 2
;				avgqual_min -> avg quality in group must be at least this value. Defaults to 0
;				avgqual_max -> avg quality in group must be less than this value. Defaults to 1d6
;				channel -> which channel (SWAVES coord) should we use to determine the groups?
;							'Ex','Ey','Ez' or 'Ex-Ey'
;				type -> the kind of wave you want to identify groups for. Defaults to 'whistler narrowband'
;						Other options include:  (for an updated list see identify_waves.pro)
;							'whistler narrowband'
;							'whistler broadband'
;							'ion acoustic'
;							'langmuir'
;							'dust'
;							'dust2'
;							'noise'
;				criteria -> select which criteria to use in selecting groups. An array of 0's and 1's
;							['freq','bw','deltaf_f','amp','type']
;							Ex. to select by only freq and bw use 
;								criteria = [1,1,0,0,0]
;
;
;  EXAMPLES:    
;                
;
;  KEYWORDS:    
;               
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:      
;
;   CREATED:  01/24/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-

function get_burst_groups,sc,tds,$
		sepmax=maxsep,$
		freqmin=min_freq,$
		freqmax=max_freq,$
		bwmin=min_bw,$
		bwmax=max_bw,$
		dfmin=min_df,$
		dfmax=max_df,$
		ampmin=min_amp,$
		ampmax=max_amp,$
		qualmin=min_qual,$
		qualmax=max_qual,$
		densmin=min_dens,$
		noevents_min=min_noevents,$
		avgqual_min=min_avgqual,$
		avgqual_max=max_avgqual,$
		channel=channel,$
		type=type,$
		criteria=crit
;		bwstddevmin=min_bwstddev,$
;		bwstddevmax=max_bwstddev,$



if ~keyword_set(crit) then crit = [1,1,1,1,1]
if total(crit) eq 0 then crit = [1,1,1,1,1]


channels = ['','Ex','Ey','Ez','Ex-Ey','','','']  ;TDS channels


if not keyword_set(sc) then begin
	sc = ''
	read,sc,prompt='Which spacecraft? (STA or STB): '
endif

if not keyword_set(type) then type = 'whistler narrowband'

;----------------
;INPUT PARAMETERS
;----------------

	;df_f < 1
	;qual > 100
	;freq < 500
	;bw < 100

;If not specified, then set up parameters to identify narrowband whistlers
if not keyword_set(sepmax) then maxsep = 5.  
if not keyword_set(freqmin) then min_freq = 0.
if not keyword_set(freqmax) then max_freq = 500.
if not keyword_set(bwmin) then min_bw = 0.
if not keyword_set(bwmax) then max_bw = 100.
if not keyword_set(dfmin) then min_df = 0.
if not keyword_set(dfmax) then max_df = 1.
if not keyword_set(ampmin) then min_amp = 0.
if not keyword_set(ampmax) then max_amp = 1d6
;if not keyword_set(bwstddevmin) then min_bwstddev = 0.
;if not keyword_set(bwstddevmax) then max_bwstddev = 1d6
if not keyword_set(qualmin) then min_qual = 100.
if not keyword_set(qualmax) then max_qual = 1d6
if not keyword_set(densmin) then min_dens = 1.      	
if not keyword_set(noevents_min) then min_noevents = 2.
if not keyword_set(avgqual_min) then min_avgqual = 0.
if not keyword_set(avgqual_max) then max_avgqual = 1d6
if not keyword_set(channel) then begin
	channel2 = 3 
	channel = 'Ez'
endif else channel2 = where(channels eq channel)

maxsep = maxsep * 60. 


t0 = tds[channel2].timeS[0]
t1 = tds[channel2].timeS[n_elements(tds[channel2].timeS)-1]


;---------------------------------------
;SEE IF WE'RE USING METHOD 1 OR METHOD 2
;---------------------------------------

if tag_exist(tds,'type') eq 0 then method = 'method1' else method = 'method2'



;----------------------------------------------------------------------------------
;IDENTIFY THE WAVES SATISFYING INPUT CONDITIONS, WHETHER THEY ARE IN GROUPS OR NOT
;----------------------------------------------------------------------------------


tmp = indgen(n_elements(tds[channel2]))

if crit[0] eq 1 then begin
	goodfreq1 = where(tds[channel2].freq ge min_freq)
	tmp = setintersection(tmp,goodfreq1)
	goodfreq2 = where(tds[channel2].freq le max_freq)
	tmp = setintersection(tmp,goodfreq2)
endif
if crit[1] eq 1 then begin
	goodbw1 = where(tds[channel2].bandwidth ge min_bw)
	tmp = setintersection(tmp,goodbw1)
	goodbw2 = where(tds[channel2].bandwidth le max_bw)
	tmp = setintersection(tmp,goodbw2)
endif
if crit[2] eq 1 then begin
	gooddf1   = where(tds[channel2].deltaf_f ge min_df)
	tmp = setintersection(tmp,gooddf1)
	gooddf2   = where(tds[channel2].deltaf_f le max_df)
	tmp = setintersection(tmp,gooddf2)
endif
if crit[3] eq 1 then begin
	goodamp1 = where(tds[channel2].amp ge min_amp)
	tmp = setintersection(tmp,goodamp1)
	goodamp2 = where(tds[channel2].amp le max_amp)
	tmp = setintersection(tmp,goodamp2)
endif
;goodbwstddev1 = where(tds[channel2].bandwidth_stddev ge min_bwstddev)
;tmp = setintersection(tmp,goodbwstddev1)
;goodbwstddev2 = where(tds[channel2].bandwidth_stddev le max_bwstddev)
;tmp = setintersection(tmp,goodbwstddev2)


if crit[4] eq 1 then begin

	if method eq 'method1' then begin
	;method 1
	goodqual1 = where(tds[channel2].quality ge min_qual)
	tmp = setintersection(tmp,goodqual1)
	goodqual2 = where(tds[channel2].quality le max_qual)
	tmp = setintersection(tmp,goodqual2)
	endif
	
	if method eq 'method2' then begin
	;method 2 
	goodqual1 = where(tds[channel2].type eq type)
	tmp = setintersection(tmp,goodqual1)
	
	endif
endif


goodwaves = tmp


;----------------------------------------------------------------------------
;TAKE THE goodwaves I'VE IDENTIFIED AND FIND OUT WHICH ONES ARE IN GROUPS
;----------------------------------------------------------------------------

if goodwaves[0] ne -1 then begin
	tdiff = shift(tds[channel2].timeD[goodwaves],-1) - tds[channel2].timeD[goodwaves]
	tdiff[n_elements(tdiff)-1] = !values.f_nan
	goodv = where(tdiff le maxsep)

	;----------------------------------------------------
	;RETURN DUMMY STRUCT IF NO GROUPS ARE FOUND IN SPECIFIED TIMERANGE
	;----------------------------------------------------
	
	if goodv[0] eq -1 then begin
		print,'#####  NO GROUPS FOUND IN SPECIFIED TIMERANGE #####
		
		grouplist = {begintimeS:'',$
			endtimeS:'',$
			begintimeD:!values.f_nan,$
			endtimeD:!values.f_nan,$
			eventdensity:!values.f_nan,$
			nevents:!values.f_nan,$
			minqual:!values.f_nan,$
			maxqual:!values.f_nan,$
			minamp:!values.f_nan,$
			maxamp:!values.f_nan,$
			minfreq:!values.f_nan,$
			maxfreq:!values.f_nan,$
			minbw:!values.f_nan,$
			maxbw:!values.f_nan,$
;			minbw_stddev:!values.f_nan,$
;			maxbw_stddev:!values.f_nan,$
			mindf_f:!values.f_nan,$
			maxdf_f:!values.f_nan,$
			avgamp:!values.f_nan,$
			avgqual:!values.f_nan,$
			avgfreq:!values.f_nan,$
			avgbw:!values.f_nan,$
;			avgbw_stddev:!values.f_nan,$
			avg_df_f:!values.f_nan,$
			notes:'',$
			notes2:'',$
			params:''}

		
		return,grouplist
	endif

	diff2 = fltarr(n_elements(tdiff))
	diff2[goodv] = 1.

	gooddiff = where(diff2 eq 1,complement=comp)
	beginv = lonarr(n_elements(gooddiff))
	endv = lonarr(n_elements(gooddiff))

	t3 = 0.
	beginv[0] = gooddiff[0]

	i=0L

	cond = ''

	CATCH, Error_status    
	IF Error_status NE 0 THEN BEGIN  
	   cond = 'STOP'  
	   CATCH, /CANCEL  
	ENDIF  


	while cond ne 'STOP' do begin
		if i ne 0 then beginv[i] = gooddiff[t3[0]]
		t2 = where(comp gt beginv[i])
		endv[i] = comp[t2[0]]
		t3 = where(gooddiff gt endv[i])

		i++
	endwhile



	goo = where(endv ne 0)
	beginv = beginv[goo]
	endv = endv[goo]


	;the original array elements in groups are now referenced by goodwaves[beginv] and goodwaves[endv]
	GbegintimeD = tds[channel2].timeD[goodwaves[beginv]]
	GendtimeD = tds[channel2].timeD[goodwaves[endv]]
	GbegintimeS = tds[channel2].timeS[goodwaves[beginv]]
	GendtimeS = tds[channel2].timeS[goodwaves[endv]]


	nelem = n_elements(beginv)
	avgfreq = fltarr(nelem)
	avgqual = fltarr(nelem)
	avgamp = fltarr(nelem)
	avgbw = fltarr(nelem)
;	avgbw_stddev = fltarr(nelem)
	avg_df_f = fltarr(nelem)
	nevents = fltarr(nelem)
	evdensity = fltarr(nelem)

	minfreq = fltarr(nelem)
	maxfreq = fltarr(nelem)
	minqual = fltarr(nelem)
	maxqual = fltarr(nelem)
	minamp = fltarr(nelem)
	maxamp = fltarr(nelem)
	minbw = fltarr(nelem)
	maxbw = fltarr(nelem)
;	minbw_stddev = fltarr(nelem)
;	maxbw_stddev = fltarr(nelem)
	mindf_f = fltarr(nelem)
	maxdf_f = fltarr(nelem)

	for i=0,nelem-1 do begin
		avgfreq[i] = total(tds[channel2].freq[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
		avgqual[i] = total(tds[channel2].quality[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
		avgamp[i] = total(tds[channel2].amp[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
		avgbw[i] = total(tds[channel2].bandwidth[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
;		avgbw_stddev[i] = total(tds[channel2].bandwidth_stddev[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
		avg_df_f[i] = total(tds[channel2].deltaf_f[goodwaves[beginv[i]:endv[i]]])/(endv[i]-beginv[i]+1)
		nevents[i] = (endv[i]-beginv[i]+1)
		evdensity[i] = 60.*nevents[i]/(tds[channel2].timeD[goodwaves[endv[i]]]-tds[channel2].timeD[goodwaves[beginv[i]]])

		minfreq[i] = min(tds[channel2].freq[goodwaves[beginv[i]:endv[i]]],/nan)
		minqual[i] = min(tds[channel2].quality[goodwaves[beginv[i]:endv[i]]],/nan)
		minamp[i] = min(tds[channel2].amp[goodwaves[beginv[i]:endv[i]]],/nan)
		minbw[i] = min(tds[channel2].bandwidth[goodwaves[beginv[i]:endv[i]]],/nan)
;		minbw_stddev[i] = min(tds[channel2].bandwidth_stddev[goodwaves[beginv[i]:endv[i]]],/nan)
		mindf_f[i] = min(tds[channel2].deltaf_f[goodwaves[beginv[i]:endv[i]]],/nan)

		maxfreq[i] = max(tds[channel2].freq[goodwaves[beginv[i]:endv[i]]],/nan)
		maxqual[i] = max(tds[channel2].quality[goodwaves[beginv[i]:endv[i]]],/nan)
		maxamp[i] = max(tds[channel2].amp[goodwaves[beginv[i]:endv[i]]],/nan)
		maxbw[i] = max(tds[channel2].bandwidth[goodwaves[beginv[i]:endv[i]]],/nan)
;		maxbw_stddev[i] = max(tds[channel2].bandwidth_stddev[goodwaves[beginv[i]:endv[i]]],/nan)
		maxdf_f[i] = max(tds[channel2].deltaf_f[goodwaves[beginv[i]:endv[i]]],/nan)
	endfor

;-----------------------------------------------------------
;Now filter by the variables that deal with group quantities
;-----------------------------------------------------------

	gooddens = where(evdensity ge min_dens)
	good_min_noevents = where(nevents ge min_noevents)
	good_min_avgqual = where(avgqual ge min_avgqual)
	good_max_avgqual = where(avgqual le max_avgqual)

	tmp = setintersection(gooddens,good_min_noevents)
	tmp = setintersection(tmp,good_min_avgqual)
	goodwaves = setintersection(tmp,good_max_avgqual)



	GbegintimeS = GbegintimeS[goodwaves]
	GendtimeS = GendtimeS[goodwaves]
	GbegintimeD = GbegintimeD[goodwaves]
	GendtimeD = GendtimeD[goodwaves]
	avgfreq = avgfreq[goodwaves]
	avgqual = avgqual[goodwaves]
	avgamp = avgamp[goodwaves]
	avgbw = avgbw[goodwaves]
;	avgbw_stddev = avgbw_stddev[goodwaves]
	avg_df_f = avg_df_f[goodwaves]
	nevents = nevents[goodwaves]
	evdensity = evdensity[goodwaves]
	minfreq = minfreq[goodwaves]
	minqual = minqual[goodwaves]
	minamp = minamp[goodwaves]
	minbw = minbw[goodwaves]
;	minbw_stddev = minbw_stddev[goodwaves]
	mindf_f = mindf_f[goodwaves]
	maxfreq = maxfreq[goodwaves]
	maxqual = maxqual[goodwaves]
	maxamp = maxamp[goodwaves]
	maxbw = maxbw[goodwaves]
;	maxbw_stddev = maxbw_stddev[goodwaves]
	maxdf_f = maxdf_f[goodwaves]


endif

;--------------------
;CREATE THE STRUCTURE
;--------------------

notes = ['max TDS capture separation time (min) for consecutive events in group = ' + strtrim(maxsep/60.,2),$
	'min allowed event density (events/min) of group = ' + strtrim(min_dens,2),$
	'min allowed number of events in group = ' + strtrim(min_noevents,2),$
	'min allowed freq (Hz) for waves in group = ' + strtrim(min_freq,2),$
	'max allowed freq (Hz) for waves in group = ' + strtrim(max_freq,2),$
	'min allowed bandwidth (Hz) for waves in group = ' + strtrim(min_bw,2),$
	'max allowed bandwidth (Hz) for waves in group = ' + strtrim(max_bw,2),$
	'min allowed df/f for waves in group = ' + strtrim(min_df,2),$
	'max allowed df/f for waves in group = ' + strtrim(max_df,2),$
	'min allowed amplitude (mV/m) for waves in group = ' + strtrim(min_amp,2),$
	'max allowed amplitude (mV/m) for waves in group = ' + strtrim(max_amp,2),$
;	'min allowed stddev of bandwidth for waves in group = ' + strtrim(min_bwstddev,2),$
;	'max allowed stddev of bandwidth for waves in group = ' + strtrim(max_bwstddev,2),$
	'min allowed quality of group = ' + strtrim(min_qual,2),$
	'max allowed quality of group = ' + strtrim(max_qual,2),$
	'min allowed average quality of group = ' + strtrim(min_avgqual,2),$
	'max allowed average quality of group = ' + strtrim(max_avgqual,2),$
	'Channel used to identify groups = ' + channel]



notes2 = ['begintimeS, endtimeS = string begin and end times of each group',$
		  'begintimeD, endtimeD = double begin and end times of each group',$
		  'eventdensity = # events/minute',$
		  'nevents = number of events in each group',$
		  'qual = quality measure from wave_quality.pro',$
		  'amp = amplitude (0-peak) in mV/m',$
		  'freq = frequency (Hz)',$
		  'bw = bandwidth (Hz)',$
;		  'bw_stddev = standard deviation of bandwidth',$
		  'df_f = deltaf/f']


struct = {t0:t0,$
		 t1:t1,$
		 maxsep:maxsep,$
		 min_dens:min_dens,$
		 min_noevents:min_noevents,$
		 min_freq:min_freq,$
		 max_freq:max_freq,$
		 min_bw:min_bw,$
		 max_bw:max_bw,$
		 min_df:min_df,$
		 max_df:max_df,$
		 min_amp:min_amp,$
		 max_amp:max_amp,$
;		 min_bwstddev:min_bwstddev,$
;		 max_bwstddev:max_bwstddev,$
		 min_qual:min_qual,$
		 max_qual:max_qual,$
		 min_avgqual:min_avgqual,$
		 max_avgqual:max_avgqual,$
		 channel:channel}


if goodwaves[0] ne -1 then begin
grouplist = {begintimeS:GbegintimeS,$
			endtimeS:GendtimeS,$
			begintimeD:GbegintimeD,$
			endtimeD:GendtimeD,$
			eventdensity:evdensity,$
			nevents:nevents,$
			minqual:minqual,$
			maxqual:maxqual,$
			minamp:minamp,$
			maxamp:maxamp,$
			minfreq:minfreq,$
			maxfreq:maxfreq,$
			minbw:minbw,$
			maxbw:maxbw,$
;			minbw_stddev:minbw_stddev,$
;			maxbw_stddev:maxbw_stddev,$
			mindf_f:mindf_f,$
			maxdf_f:maxdf_f,$
			avgamp:avgamp,$
			avgqual:avgqual,$
			avgfreq:avgfreq,$
			avgbw:avgbw,$
;			avgbw_stddev:avgbw_stddev,$
			avg_df_f:avg_df_f,$
			notes:notes,$
			notes2:notes2,$
			params:struct}
endif else begin
grouplist = {begintimeS:'',$
			endtimeS:'',$
			begintimeD:!values.f_nan,$
			endtimeD:!values.f_nan,$
			eventdensity:!values.f_nan,$
			nevents:!values.f_nan,$
			minqual:!values.f_nan,$
			maxqual:!values.f_nan,$
			minamp:!values.f_nan,$
			maxamp:!values.f_nan,$
			minfreq:!values.f_nan,$
			maxfreq:!values.f_nan,$
			minbw:!values.f_nan,$
			maxbw:!values.f_nan,$
;			minbw_stddev:!values.f_nan,$
;			maxbw_stddev:!values.f_nan,$
			mindf_f:!values.f_nan,$
			maxdf_f:!values.f_nan,$
			avgamp:!values.f_nan,$
			avgqual:!values.f_nan,$
			avgfreq:!values.f_nan,$
			avgbw:!values.f_nan,$
;			avgbw_stddev:!values.f_nan,$
			avg_df_f:!values.f_nan,$
			notes:notes,$
			notes2:notes2,$
			params:struct}
endelse


return,grouplist

end
