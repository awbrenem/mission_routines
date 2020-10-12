;+
;*****************************************************************************************
;
;  FUNCTION :   get_group_elements.pro
;  PURPOSE  :   Returns the "reducedlist" structure. This contains all of the TDS captures
;				that fall within the group burst times from get_group_elements.pro
;
;  CALLED BY:   
;               
;
;  CALLS:
;               
;
;  REQUIRES:    get_burst_groups.pro
;               
;
;  INPUT:		sc -> 'STA' or 'STB'
;				tds -> structure from the program get_TDS.pro that contains all of the TDS captures
;					   for a certain timerange
;				groups -> structure from program get_burst_groups.pro containing all of the TDS
;						groups for specified timerange and under various conditions
;				channel -> which channel (SWAVES coord) should we use to determine the groups?
;							'Ex','Ey','Ez' or 'Ex-Ey'
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
;
;   CREATED:  01/25/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-


function get_group_elements,sc,tds,groups,channel=channel


tst = size(groups)
tst1 = size(tds)

;----------------------------------------------------
;RETURN DUMMY STRUCT IF NO GROUPS ARE FOUND IN SPECIFIED TIMERANGE
;----------------------------------------------------

if tst[0] eq 0 or tst1[0] eq 0 then begin
	print,'###### NO GROUPS...RETURNING ######'
	
		reducedlist = {times:!values.f_nan,$
		   timed:!values.f_nan,$
		   source:'',$
		   freq:!values.f_nan,$
		   quality:!values.f_nan,$
		   amp:!values.f_nan,$
		   bandwidth:!values.f_nan,$
		   bandwidth_stddev:!values.f_nan,$
		   deltaf_f:!values.f_nan,$
		   sample_rate:'',$
		   npoints:!values.f_nan,$
		   sample_length:!values.f_nan,$
		   filter:!values.f_nan,$
		   eventid:!values.f_nan,$
		   header:!values.f_nan}

	
	return,reducedlist
endif


if not keyword_set(sc) then begin
	sc = ''
	read,sc,prompt='Which spacecraft? (STA or STB): '
endif
if not keyword_set(channel) then begin
	channel = ''
	read,channel,prompt='Which data channel? (Ex,Ey,Ez,Ex-Ey): '
endif


channels = ['','Ex','Ey','Ez','Ex-Ey','','','']
channel2 = where(channels eq channel)


t0 = tds[channel2].timeS[0]
t1 = tds[channel2].timeS[n_elements(tds[channel2].timeS)-1]


goodb = !values.f_nan
q=0.
for i=0,n_elements(groups.begintimes)-1 do begin

	tmp = where((tds[channel2].timed ge groups.begintimed[i]) and (tds[channel2].timed le groups.endtimed[i]))
	if tmp[0] ne -1 then goodb = [goodb,tmp]  
	q = q + n_elements(tmp)

endfor


tst = size(goodb) ;if there's only a single matched element, make sure it's not a NaN
if tst[0] eq 0 then tst2 = reform(goodb[0]) else tst2 = 1.

if finite(tst2) ne 0 then begin
	goodb = goodb[1:n_elements(goodb)-1]

	reducedlist = {times:tds[channel2].times[goodb],$
		   timed:tds[channel2].timed[goodb],$
		   source:tds[channel2].source[goodb],$
		   freq:tds[channel2].freq[goodb],$
		   quality:tds[channel2].quality[goodb],$
		   amp:tds[channel2].amp[goodb],$
		   bandwidth:tds[channel2].bandwidth[goodb],$
		   bandwidth_stddev:tds[channel2].bandwidth_stddev[goodb],$
		   deltaf_f:tds[channel2].deltaf_f[goodb],$
		   sample_rate:tds[channel2].sample_rate[goodb],$
		   npoints:tds[channel2].npoints[goodb],$
		   sample_length:tds[channel2].sample_length[goodb],$
		   filter:tds[channel2].filter[goodb],$
		   eventid:tds[channel2].eventid[goodb],$
		   header:tds[channel2].header[goodb]}
endif else begin
	reducedlist = {times:!values.f_nan,$
		   timed:!values.f_nan,$
		   source:'',$
		   freq:!values.f_nan,$
		   quality:!values.f_nan,$
		   amp:!values.f_nan,$
		   bandwidth:!values.f_nan,$
		   bandwidth_stddev:!values.f_nan,$
		   deltaf_f:!values.f_nan,$
		   sample_rate:'',$
		   npoints:!values.f_nan,$
		   sample_length:!values.f_nan,$
		   filter:!values.f_nan,$
		   eventid:!values.f_nan,$
		   header:!values.f_nan}


endelse

return,reducedlist
end