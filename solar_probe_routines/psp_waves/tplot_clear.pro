;+
;PROCEDURE: tplot_clear.pro
;
;PURPOSE: to delete tplot variables
;
;ARGUMENTS:
;	ID    ->  ID(s) of tplot quanties to remov
;
;KEYWORDS: 
;       ALL   /   Clear ALL tplot variables
;       STR   ->  String for TPLOT variables to remove
;
;CALLING SEQUENCE: tplot_clear,/all
;                  tplot_clear,[1,2,3,4]
;
;NOTES: 
;      None
;
;CREATED BY: John Dombeck Oct 24 2001
;
;LAST MODIFICATION: 
;       10/24/2001-J.Dombeck    Creation
;       07/07/2004-J.Dombeck    Added STR keyword
;       11/21/2007-J.Dombeck    Fixed crash when no tplot qtys
;-
;INCLUDED MODULES:
;   tplot_clear
;
;LIBRARIES USED:
;   tplot_com
;
;DEPENDANCIES
;   data_type
;   store_data
;
;-



;*** MAIN *** : * TPLOT_CLEAR *

pro tplot_clear,id,all=all,str=str
@tplot_com.pro

  num=n_elements(id)
  if (num ne 0 and keyword_set(all)) then begin
    message,"Passed in ID(s) with keyword ALL set",/cont
    return
  endif

  if (num eq 0 and not keyword_set(all) and not keyword_set(str)) then begin
    message,"ID(s) required",/cont
    return
  endif

  if data_type(data_quants) ne 8 then return  ; No Tplot Quantities

  names=data_quants.name
  tnum=n_elements(names)

  if keyword_set(all) then begin
    for xx=1,tnum-1 do store_data,names[xx],/delete
    return
  endif

  if keyword_set(str) then begin
    if data_type(str) ne 7 then begin
      message,'STR requires string',/cont
      return
    endif
    tplot_names,names=tnames
    match=strpos(tnames,str)
    ids=where(match ne -1,cnt)
    if cnt ne 0 then tplot_clear,tnames[ids]
    return
  endif

  if data_type(id) eq 7 then begin
    for xx=0,num-1 do store_data,id[xx],/delete
    return
  endif

  bogus=where(id ne fix(id),cnt)
  if cnt ne 0 then begin
    message,"ID(s) require integers",/cont
    return
  endif

  bogus=where(id lt 1,cnt)
  if cnt ne 0 then begin
    message,"ID(s) > 0 required",/cont
    return
  endif

  bogus=where(id gt tnum-1,cnt)
  if cnt ne 0 then begin
    message,"ID(s) out of range",/cont
    return
  endif

  for xx=0,num-1 do store_data,names[id[xx]],/delete

return
end        ;*** MAIN *** : * TPLOT_CLEAR *
