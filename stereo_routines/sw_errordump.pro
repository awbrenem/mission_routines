;+
;PROGRAM: sw_errordump.pro
;
;PURPOSE: Uses TMlib error reporting routines to generate useful error dumps.
;
;ARGUMENTS:
;       STREAMID    ->  TMlib stream ID (integer) of the broken stream
;
;KEYWORDS:
;       ERRORID     ->  optional string identifier
;                       (useful for bug tracing, ie. 'program.pro, line ###'
;
;RETURNS: N/A
;
;CALLING SEQUENCE:
;       IDL> sw_errordump, streamid
;       IDL> sw_errordump, streamid, errorid='Some useful keyword, etc.'
;
;NOTES:
;       (none)
;-
;CREATED BY:    Kris Kersten, Jan 2007
;
;MODIFICATION HISTORY:
;       01/01/2007  created, KK
;
;INCLUDED MODULES:
;       sw_errordump
;
;LIBS USED:
;
;
;DEPENDENCIES:
;       *REQUIRES TMlib SIDL <http://homepage.mac.com/swaves>
;-

pro sw_errordump,streamid,errorid=errorid

    e_err=0
    e_buff=800
    e_dump=''

    print,""
    print,"----- BEGIN ERROR DUMP -----"
    if keyword_set(errorid) then print,"ERROR IN:  ",errorid
    e_err=TM_Error_Stack(streamid,"code",e_dump,e_buff,size)
    print,"    CODE:  ",e_dump
    e_err=TM_Error_Stack(streamid,"name",e_dump,e_buff,size)
    print,"    NAME:  ",e_dump,2
    e_err=TM_Error_Stack(streamid,"description",e_dump,e_buff,size)
    print,"    DESC:  "+strtrim(e_dump,2)
    e_err=TM_Error_Stack(streamid,"message",e_dump,e_buff,size)
    print,"    MSG :  "+strtrim(e_dump,2)
    print,"------ END ERROR DUMP ------"
    print,""

end