;+
;*****************************************************************************************
;
;  FUNCTION :   get_TDS.pro
;  PURPOSE  :   returns a structure of all the STEREO TDS burst captures b/t t0 and t1
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
;  INPUT:		t0 -> initial time (ex '2007-04-08/22:09:44.793')
;				t1 -> final time (ex '2007-04-08/22:59:20.832')
;				file -> name of the file output from wave_quality.pro
;
;
;  EXAMPLES:    This program is meant to be used in conjunction with
;				get_burst_groups.pro and get_group_elements.pro
;
;				sc = 'STA'
;				t0 = '2007-01-01/00:00:00'
;				t1 = '2007-12-31/24:00:00'
;				file2 = 'STEREO_TDS_STA_q_waveclass_2007_complete.txt'
;				tds = identify_waves(t0,t1,file='~/Desktop/code/Aaron/datafiles/stereo/wave_quality/'+file2)
;				groups = get_burst_groups('STA',tds,channel='Ex')
;				elements = get_group_elements('STA',tds,groups,channel='Ex')
;
;               
;
;  KEYWORDS:    
;               
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:	    structure has form (STRUCT    = -> <Anonymous> Array[8]), one embedded
;					structure for each time (8 possible channels)
;				
;
;   CREATED:  01/22/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-


function get_TDS,t0,t1,file=file

init_path

;####TEST INPUT#######
;initial and final times
;t0 = '2007-04-08/22:09:44.793'
;t1 = '2007-04-08/22:59:20.832'

;test file
;file = '~/Desktop/STEREO_TDS_STA_q_waveclass_20070408_to_20070409.txt'
;#######################

if not keyword_set(t0) then begin
	t0 = ''
	read,t0,prompt='Enter begin time (yyyy-mm-dd/hh:mm:ss.msc): '
endif
if not keyword_set(t1) then begin
	t1 = ''
	read,t1,prompt='Enter end time (yyyy-mm-dd/hh:mm:ss.msc): '
endif
if not keyword_set(file) then begin
	file = ''
	file = dialog_pickfile(path=!data.stereo.wavequality)
endif


header = strarr(20)


openr,lun,file,/get_lun

;--------------------------------------------------
;Read in the header file and get to the actual data
;--------------------------------------------------

junk = ''
for i=0,18 do begin
	readf,lun,junk
	header[i] = junk
endfor


format = ''
readf,lun,format
format = strmid(format,8,1000)

header[19] = format

readf,lun,junk
hoo = ''
readf,lun,hoo



goo = ''
times = strarr(1000000)
timed = dblarr(1000000)
numb = dblarr(1000000)

cond = ''


CATCH, Error_status    
IF Error_status NE 0 THEN BEGIN  
   cond = 'STOP'  
   CATCH, /CANCEL  
ENDIF  

;------------------------------
;Get the times for each channel
;------------------------------

jj=0L
while cond ne 'STOP' do begin
	readf,lun,goo
	times[jj] = strtrim(strmid(goo,0,25),2)	
	timed[jj] = time_double(times[jj])
	for b=0,2 do readf,lun,junk
	jj++
endwhile


goo = where(times eq '')
times = times[0:goo[0]-2]
timed = timed[0:goo[0]-2]


close,lun
free_lun,lun

;-------------------------------------------------------------------------------
;reopen file so I can read the actual data from all the channels (for each time)
;-------------------------------------------------------------------------------

openr,lun2,file,/get_lun
for i=0,21 do readf,lun2,junk

goo = where((timed ge time_double(t0)) and (timed le time_double(t1)))

if goo[0] eq -1 then print,'NO CAPTURES IN SELECTED TIMERANGE. RETURNING DUMMY STRUCTURE'
;line numbers of the good captures
;numb = goo*8L
numb = goo*4L


;Reduced arrays with only the elements we are interested in
if goo[0] ne -1 then timeSf = times[goo]
if goo[0] ne -1 then timeDf = timed[goo]

;----------------------------
;Now get the rest of the data
;----------------------------

;ch = strarr(n_elements(numb),15,8)  ;;array with all the data for each channel
ch = strarr(n_elements(numb),14,4)  ;;array with all the data for each channel

;14 elements including time

cond = ''

CATCH, Error_status    
IF Error_status NE 0 THEN BEGIN  
   cond = 'STOP'  
   CATCH, /CANCEL  
ENDIF  

while cond ne 'STOP' do begin
	;get to the first element
	if numb[0] ne 0 then begin
		for qq=0D,numb[0]-1 do begin
			junk = ''
			readf,lun2,junk
		endfor
	endif
    
;read in the 4 channels for each burst time
	for bb=0D,n_elements(numb)-1 do begin
		for uu=0,3 do begin
			junk = strarr(14)
			readf,lun2,junk,format=format
			;print,junk
			ch[bb,*,uu] = junk
		endfor
	endfor
	cond ='STOP'
endwhile


close,lun2
free_lun,lun2

nelem = n_elements(numb)

;create a structure with all the data
dummy = {timeS:replicate('',nelem),$
	  timeD:replicate(0D,nelem),$
	  source:replicate('',nelem),$
	  freq:replicate(0L,nelem),$
	  quality:replicate(0L,nelem),$
	  amp:replicate(0L,nelem),$
	  bandwidth:replicate(0L,nelem),$
	  deltaf_f:replicate(0.,nelem),$
	  sample_rate:replicate(0L,nelem),$
	  npoints:replicate(0L,nelem),$
	  sample_length:replicate(0.,nelem),$
	  filter:replicate(0L,nelem),$
	  eventID:replicate(0L,nelem),$
	  triggersource:replicate('',nelem),$
	  tdsQ:replicate(0L,nelem),$
	  header:replicate('',n_elements(header))}

 
struct = replicate(dummy,4)  ;one for each data channel


if goo[0] eq -1 then return,struct

for i=0,3 do begin
	struct[i].timeS = timeSf
	struct[i].timeD = timeDf
	struct[i].source = ch[*,1,i]
	struct[i].freq = long(ch[*,2,i])
	struct[i].quality = float(ch[*,3,i])
	struct[i].amp = long(ch[*,4,i])
	struct[i].bandwidth = long(ch[*,5,i])
	struct[i].deltaf_f = float(ch[*,6,i])
	struct[i].sample_rate = long(ch[*,7,i])
	struct[i].npoints = long(ch[*,8,i])
	struct[i].sample_length = float(ch[*,9,i])
	struct[i].filter = long(ch[*,10,i])
	struct[i].eventID = long(ch[*,11,i])
	struct[i].triggersource = ch[*,12,i]
	struct[i].tdsQ = float(ch[*,13,i])
	struct[i].header = header
endfor

return,struct

end
