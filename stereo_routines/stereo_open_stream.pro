;Open a STEREO TMLIB data stream.
;Returns evid and UR8start.
;endtime is optional.


pro stereo_open_stream,sid,mission,spcr,inst,ev_type,starttime,nchannel,endtime=endtime,$
                       evid,ev_err,UR8start,UR8end

;-----------------------
;-try to open the stream
;-----------------------

print,"Selecting Domain:  "+mission+"/"+spcr+"/"+inst+"/"+ev_type
err=TM_Select_Domain(sid,mission,spcr,inst,ev_type)
if err ne 0 then begin
    sw_errordump,sid,errorid="DOMAIN SELECTION"
    return
endif

;-----------------------
;-get a start time
;-----------------------


;if not keyword_set(starttime) then read,starttime,prompt="Enter start time (yyyy-mm-dd hh:mm:ss.msc): "

err=TM_UR8_from_string_comparison(UR8start,starttime)
if err ne 0 then begin
    sw_errordump,sid,errorid="TIME CONVERSION"
    return
endif

print,'start time conversion OK...'
print,"  UR8:  ",UR8start
print,''

if not keyword_set(endtime) then endtime = '2020-01-01 00:00:01'   ;large endtime to make sure all events are gotten

err=TM_UR8_from_string_comparison(UR8end,endtime)
if err ne 0 then begin
    sw_errordump,sid,errorid="TIME CONVERSION"
    return
endif

print,'end time conversion OK...'
print,"  UR8:  ",UR8end
print,''

if (UR8end - UR8start) le 0 then begin
	sw_errordump,sid,errorid="BAD TIMES"
	RETURN
endif

;---------------------------
;-set up the stream timerange (from start time to EOI)
;---------------------------

err=TM_Select_Stream_TimeRange(sid,UR8start,-1)
if err ne 0 then begin
    sw_errordump,sid,errorid="TIME RANGE SELECTION"
    return
endif

;-------------------------------------------------------------
;Get an event to start things off and make sure it's channel 1
;-------------------------------------------------------------

while nchannel ne 1 do begin
    ev_err=TM_Find_Event(sid)
    if ev_err ne 0 then begin
        sw_errordump,sid,errorid='COULD NOT FIND FIRST EVENT'
        return
    endif
    err=tm_get_item_i4(sid,"Channel",nchannel,1l,n_items)
    err=err+tm_get_item_i4(sid,"EventID",evid,1l,n_items)
	err=err+tm_get_item_i4(sid,"Nsamples",npoints,1l,n_items)
    if err ne 0 then begin
        sw_errordump,sid,errorid='FIRST EVENT: DATA ACQUISITION'
        return
    endif
endwhile


end
