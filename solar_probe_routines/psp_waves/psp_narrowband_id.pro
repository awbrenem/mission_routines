;
;
;
;
;
;
;
;
;
;
;

PRO psp_narrowband_id,t0,tf,file

  tmp0=t0
  ;15 lines before data
  header=strarr(15)

  openr,lun,file,/get_lun
  ;the function file_lines(file) returns an array of length 1
  rows = file_lines(file)
  ;make rows just the value of file_lines, so a scalar
  rows = rows[0]

  ;--------------------------------------------------
  ;Read in the header file and get to the actual data
  ;--------------------------------------------------

  junk = ''
  for i=0,n_elements(header)-4 do begin
    readf,lun,junk
    header[i] = junk
  endfor

  format = ''
  readf,lun,format
  format = strmid(format,8,1000)

  ;the last entry in header should be the format
  header[-1] = format

  ;the titles of each column
  for i=0,1 do readf,lun,junk

  cond = ''

  ;CATCH, Error_status
  ;IF Error_status NE 0 THEN BEGIN
  ;  print,!Error_state.msg
  ;  cond = 'STOP'
  ;  CATCH, /CANCEL
  ;ENDIF

  ;------------------------------;
  ;Get the times for each channel;
  ;------------------------------;

  array1 = ''
  line = ''
  while ~EOF(lun) do begin
    READF, lun, line
    array1 = [array1, line]
  endwhile

  array1 = array1[0:n_elements(array1)-2]
  nTDS = n_elements(array1)

  x = {duration: fltarr(nTDS),$
    quality: intarr(nTDS) $
  }

  x.duration[*] = float(strtrim(strmid(array1[*],89,6)))
  x.quality[*] = long(strtrim(strmid(array1[*],147,2)))
  ;print, x.duration[6]
  array2 = strarr(nTDS)
  j=0
  for i=0,nTDS-1 do begin
    if (x.quality[i] ne -1) then begin ;(x.duration[i] ne -1) and (x.duration[i] ne 0) and 
      array2[j] = array1[i]
      j+=1
    endif
  endfor
  
  bwhere = where(array2 eq '')
  array2 = array2[0:bwhere[0]-1]

  nTDS = n_elements(array2)

  y = {date: strarr(nTDS),$
    duration: fltarr(nTDS),$
    amp: fltarr(nTDS) $
  }

  y.date[*] = strmid(array2[*], 1,23)
  y.duration[*] = float(strtrim(strmid(array2[*],89,6)))
  y.amp[*] =float(strtrim(strmid(array2[*],60,7)))
  array3 = strarr(nTDS)

  j=0
  for i=0,nTDS-3 do begin
    if (y.date[i] eq y.date[i+1]) and (y.date[i] eq y.date[i+2]) then begin
      max = max(y.amp[i:i+2], sub)
      array3[j] = array2[i + sub]
      j+=1
      i+=3
    endif else if (y.date[i] eq y.date[i+1]) and (y.date[i] ne y.date[i+2]) then begin
      max = max(y.amp[i:i+1], sub)
      array3[j] = array2[i + sub]
      j+=1
      i+=2
    endif else begin
      array3[j] = array2[i]
      j+=1
    endelse
  endfor

  cwhere = where(array3 eq '')
  array3 = array3[0:cwhere[0]-1]
  nTDS = n_elements(array3)

  z = {date: strarr(nTDS), $
    source: strarr(nTDS), $
    freq: intarr(nTDS),$
    f_fce: fltarr(nTDS),$
    amp: fltarr(nTDS),$
    bw: fltarr(nTDS),$
    df_f: fltarr(nTDS),$
    duration: fltarr(nTDS),$
    wave_angle: fltarr(nTDS),$
    fce: intarr(nTDS),$
    tper_tpar: fltarr(nTDS),$
    heatflux: fltarr(nTDS),$
    coev:fltarr(nTDS),$
    countev:fltarr(nTDS) $
  }

  z.date[*] = strmid(array3[*], 1,23)
  z.source[*] = strmid(array3[*], 31, 1)
  z.freq[*] = long(strtrim(strmid(array3[*], 39, 4)))
  z.f_fce[*] = float(strtrim(strmid(array3[*],53, 4)))
  z.amp[*] = float(strtrim(strmid(array3[*],60,7)))
  z.bw[*] = float(strtrim(strmid(array3[*],69,7)))
  z.df_f[*] = float(strtrim(strmid(array3[*],79,5)))
  z.duration[*] = float(strtrim(strmid(array3[*],89,6)))
  z.wave_angle[*] = float(strtrim(strmid(array3[*],103,6)))
  z.coev[*] = float(strtrim(strmid(array3[*],117,7)))
  z.countev[*] = float(strtrim(strmid(array3[*],132,7)))

  for i=0,nTDS-1 do begin
    z.fce[i] = z.freq[i]/z.f_fce[i]
  endfor


  format = '(a24,5x,a3,7x,i4,10x,f4.2,3x,f7.2,2x,f7.2,3x,f5.3,5x,f6.3,8x,f6.3,10x,i4, 8x, f7.2, 8x, f7.2)'
  openw,outlun,'~/Desktop/psp_identify_narrowband/'+strmid(tmp0,0,10)+'.txt',/get_lun
  for i=0,12 do printf,outlun, header[i]
  printf, outlun, format
  printf,outlun,'        Time                  Src  Frequency(Hz)   f/fce   Amp(pk-pk)   BW    df/f   duration (sec)   wave_angle(deg)   fce(Hz)    co(eV)   count(eV)'
  for j=0,nTDS-1 do begin
    if z.date[j] ne '' then begin
      printf, outlun, format=format, z.date[j], z.source[j], z.freq[j], z.f_fce[j], z.amp[j], z.bw[j], z.df_f[j], z.duration[j], z.wave_angle[j], z.fce[j], z.coev[j], z.countev[j]
    endif

  endfor

  close,outlun
  free_lun,outlun, /FORCE

end ;program stereo_narrowband_id