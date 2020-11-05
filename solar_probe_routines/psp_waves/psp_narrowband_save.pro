;+
;PROGRAM:   psp_narrowband_save.pro
;
;ARGUMENTS:
;
;
;KEYWORDS:
;
;
;RETURNS:
;   WRITES OUTPUT FILES:
;     'STEREO_TDS_STB_q_waveclass.txt' --> for TDS quality data on STB
;
;
;
;INCLUDED MODULES:
;       swtds_lff
;
;LIBS USED:
;       (none)
;
;DEPENDENCIES:
;
;
;CREATED BY:    Zac A. Cohen
;
;
;HISTORY:
;   Created: 06/29/2017 - ZAC
;
;-

function psp_narrowband_save,savelun,path,mission,inst,version,ev,count_a,t0,tf,$
  bgn=bgn,mid=mid,ed=ed
  format1 =    '(a24,5x,a3,7x,i4,10x,f4.2,3x,f7.2,2x,f7.2,3x,f5.3,5x,f6.3,8x, f6.3, 8x,f7.2,8x,f7.2, 8x, i2)'
                ;date  srce  freq   f/fce    Amp     BW    d/df   duration                                quality
                ;date  srce  freq   f/fce    Amp     BW    d/df   duration   waveangle   coev   countev   quality
  if keyword_set(bgn) then begin


    openu,savelun,path

    printf,savelun,'Run started: ' + systime()
    ;printf,savelun,'TDS waveform capture data for: ' + strupcase(strmid(ev_type,0,3))
    printf,savelun,'Spacecraft: Parker Solar Probe'
    ;printf,savelun,'Quality or honesty: ' + strupcase(hq)
    printf,savelun,'Instrument: ' + inst
    printf,savelun,'Mission: ' + mission
    printf,savelun,'Version: ' + version
    printf,savelun,'Channels = [Ex, Ey, Ez]'
    printf,savelun,'Freq = Waveform Average Frequency in Hz'
    printf,savelun,'Amp = Wave maximum amplitude in mV/m (peak-peak)'
    printf,savelun,'Bw = bandwidth of peak power in Hz'
    printf,savelun,'df/f = normalized bandwidth'
    ;printf,savelun,'drate = data rate'
    ;printf,savelun,'npts = number of samples in burst'
    ;printf,savelun,'blen = burst length in sec'
    ;printf,savelun,'maxf = data removed above this frequency '
    ;printf,savelun,'ID = official event ID number'
    ;printf,savelun,'ts = channel TDS triggered off of'
    ;printf,savelun,'TDSQ = quality of event as defined onboard the sc'
    printf,savelun,'format: ' + format1

    printf,savelun,''
    print,'Opened output files: '+path
    print,''


    printf,savelun,'        Time                  Src  Frequency(Hz)   f/fce   Amp(pk-pk)   BW    df/f   duration (sec) waveangle(deg)      coeV          counteV    quality'

    close,savelun

    return,path
  endif
  
  if keyword_set(mid) then begin
    openu,savelun,path
    nullstr = string(replicate(32B,strlen(ev.scetstr)))
    header = strarr(file_lines(path))

    readf,savelun,header    ;(a24,5x,a3,7x,i4,10x,f4.2,3x,f7.2,2x,f7.2,3x,f5.3,5x,f6.3,8x, f6.3, 8x,f7.2,8x,f7.2, 8x, i2)
    ;psp_waves_ev_check,ev
    
    timestr = strarr(4)
    timestr[0] = ev.scetstr
    timestr[1:3] = nullstr
    
    printf,savelun,format=format1,timestr[0],ev.source[0],ev.freq[0],ev.f_fce[0],(0.5)*round(2*ev.max_amp[0]),ev.bw[0],ev.df_f[0],ev.duration[0],$
       ev.norm_angle[0],ev.coev[0],ev.countev[0], min(ev.quality[0,*]) ;chn 0
    
    printf,savelun,format=format1,timestr[0],ev.source[1],ev.freq[1],ev.f_fce[1],(0.5)*round(2*ev.max_amp[1]),ev.bw[1],ev.df_f[1],ev.duration[1],$
       ev.norm_angle[0],ev.coev[0],ev.countev[0], min(ev.quality[1,*]) ;chn 1
    
    printf,savelun,format=format1,timestr[0],ev.source[2],ev.freq[2],ev.f_fce[2],(0.5)*round(2*ev.max_amp[2]),ev.bw[2],ev.df_f[2],ev.duration[2],$
       ev.norm_angle[0],ev.coev[0],ev.countev[0], min(ev.quality[2,*]) ;chn 2
    
    ;printf,savelun,format=format1,nullstr,ev.source[3],ev.freq[3],ev.f_fce[3],(0.5)*round(2*ev.max_amp[3]),ev.bw[3],ev.df_f[3],ev.duration[3], min(ev.quality[3,*]) ;chn 3

    
    print,'CAPTURE NUMBER ' + strtrim(count_a,2) + ' : ' + ev.scetstr
    close,savelun
    return,path
  endif
  
  if keyword_set(ed) then begin
    openu,savelun,path

    header = strarr(file_lines(path))
    readf,savelun,header
    printf,savelun,'Run terminated: ' + systime()
    print,""
    close,savelun,exit_status=fstatus1
    print,'Wrote ' + strcompress(string(count_a),/remove_all) + ' events to file: ' + path
    print,'Output files closed with status: ',fstatus1 ;,fstatus2
    ;goo = TM_UR8_to_string_comparison(ev.scetUR8,tf)
    ;print, goo
    print,'Final time saved: ' + tf
    ;err = tm_close(StreamId)

    sttm = strmid(t0,0,4) + strmid(t0,5,2) + strmid(t0,8,2)
    edtm = strmid(tf,0,4) + strmid(tf,5,2) + strmid(tf,8,2)

    goo = strpos(path,'.txt')
    patht = strmid(path,0,goo)
    ;print, patht
    path2 = patht + '_' + sttm + '_to_' + edtm +'.txt'

    file_move,path,path2,/overwrite

    close,savelun
    free_lun,savelun

    return,path2
  endif

end   ;procedure stereo_narrowband_save