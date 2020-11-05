;+
    ; :Author: benshort
    ;-
    ;+
    ; :Description:
    ;    A procedure created in order to scroll through
    ;    Parker Solar Probe 3.5 second captures.
    ;     
    ;    the CHN parameter allows you to choose which antenna to look at.
    ;    the LVL parameter specifies what level data you want to look at.
    ;
    ; :Params:
    ;    t0 - date for waveforms
    ;    lvl - choose what level
    ;    chn -
    ;    
    ; :Dependencies:
    ;    
    ;    SPEDAS bleeding edge
    ;    tplot_clear.pro
    ;    colorbar.pro
    ;    
    ; :Author: benshort
    ; 
    ; LAST MODIFIED: 2020-05-21
    ; 
    ; ID:short186
    ; PW:poLuSHISheoc
    ; 
    ;-
pro psp_tds, t0, lvl=lvl, chn=chn
  
 ; setenv, 'PSP_STAGING_ID=short186'
  ;setenv, 'PSP_STAGING_PW=poLuSHISheoc'
  
  setenv, 'PSP_STAGING_ID=cattell'
  setenv, 'PSP_STAGING_PW=flds@psp'
  
  if getenv('PSP_STAGING_ID') EQ '' then begin 
    id = ''
    READ, id, PROMPT='Enter Staging ID: '
    setenv, 'PSP_STAGING_ID='+id
  endif  
  
  if getenv('PSP_STAGING_PW') EQ '' then begin 
    pw = ''
    READ, pw, PROMPT='Enter Staging PW: '
    setenv, 'PSP_STAGING_PW='+pw
  endif 
  
  if ~keyword_set(lvl) then begin
     lvl = 2
     print, ''
     print, "No level specified, using level 2 data."
     print, ''
     wait, 3
  endif
  
  if ~keyword_set(chn) then begin
     chn = 'e'
     print, ''
     print, "No channel specified, using electric field data."
     print, ''
     wait, 3
  endif
  
  if strlen(t0) eq 10 then t0=t0+'/00:00:00'
  
  t_tmp = t0
  
  date_dbl = time_double(t0)
  
  hr_sum = long(strmid(t0,11,2))+long(strmid(t0,14,2))+long(strmid(t0,17,2))
  
  if hr_sum gt 0 then begin
     
     if long(strmid(t0,11,2)) ge 0 and long(strmid(t0,11,2)) lt 6 then hr_tmp = '00:00:00'
     if long(strmid(t0,11,2)) ge 6 and long(strmid(t0,11,2)) lt 12 then hr_tmp = '06:00:00'
     if long(strmid(t0,11,2)) ge 12 and long(strmid(t0,11,2)) lt 18 then hr_tmp = '12:00:00'
     if long(strmid(t0,11,2)) ge 18 and long(strmid(t0,11,2)) lt 24 then hr_tmp = '18:00:00'
     
  endif else hr_tmp = '00:00:00'
  
  t_tmp = strmid(t_tmp, 0, 11)+hr_tmp
  
  timespan, t_tmp, 0.25
  t = t_tmp
  
  Bt = 3.49525 ;burst length (sec)
  sr = 150000. ;samples/sec
  pi = !PI
  q = 1.60217662*10.^(-19) ; coulombs
  mass_e = 9.10938356*10.^(-31) ;kg
  
  e=0
  jhnycash = 0
  
  datesave = -1
  ;b=''
  
  move_day:
  
  if lvl eq 1 then begin
  
    if chn eq 'e' then begin
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_1',/load
      get_data, 'spp_fld_dfb_dbm_1_dbm_data', data = dbm_dat_1
  
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_2',/load
      get_data, 'spp_fld_dfb_dbm_2_dbm_data', data = dbm_dat_2
  
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_3',/load
      get_data, 'spp_fld_dfb_dbm_3_dbm_data', data = dbm_dat_3
      
      fieldtype = 'Electric Field (adc counts) '
      color = 'blue'
      channel = 'E'
      length1 = 1
      lengthz = 1
    endif else if chn eq 'b' then begin
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_4',/load
      get_data, 'spp_fld_dfb_dbm_4_dbm_data', data = dbm_dat_1
  
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_5',/load
      get_data, 'spp_fld_dfb_dbm_5_dbm_data', data = dbm_dat_2
  
      spp_fld_make_or_retrieve_cdf,'dfb_dbm_6',/load
      get_data, 'spp_fld_dfb_dbm_6_dbm_data', data = dbm_dat_3
      
      fieldtype = 'Magnetic Field (adc counts) '
      color = 'goldenrod'
      channel = 'B'
      length1 = 1
      lengthz = 1
      
    endif
  
    spp_fld_make_or_retrieve_cdf,'mago_survey', /load
    get_data, 'spp_fld_mago_survey_nT_mag', data=mag
    get_data, 'spp_fld_mago_survey_mag_bx_nT', data=magx
    get_data, 'spp_fld_mago_survey_mag_by_nT', data=magy
    get_data, 'spp_fld_mago_survey_mag_bz_nT', data=magz
    

    
    magx = magx.y
    magy = magy.y
    magz = magz.y
  
  endif else if lvl eq 2 then begin
                                                                ;THIS IS A CRITICAL LINE THAT NEEDS TO BE CHANGED PER USER, AUTHOR PREFERED DATA BE WRITTEN TO AND READ FROM
    ;setenv, 'ROOT_DATA_DIR=/Volumes/500GB/Users/benshort/data/' ;EXTERNAL HARD DRIVE, IF YOU WANT DATA WRITTEN ELSEWHERE YOU MUST CHANGE THIS DIRECTORY.
                                                                ;IF YOU WANT DATA SENT TO AND FROM YOUR DEFAULT DIRECTORY, COMMENT OUT THIS LINE.
    setenv, 'PSP_STAGING_DIR=/Users/aaronbreneman/Desktop/data/PSP/
    
    if chn eq 'e' then begin
      
      get_timespan, tget                            
      spp_fld_load, trange=tget, type='dfb_dbm_dvac',/no_staging
      
      get_data, 'psp_fld_l2_dfb_dbm_dvac12', data = dbm_dat_1
      get_data, 'psp_fld_l2_dfb_dbm_dvac34', data = dbm_dat_2
      get_data, 'psp_fld_l2_dfb_dbm_dvacz', data = dbm_dat_3
      
      fieldtype = 'Electric Field (mV/m)'
      color = 'blue'
      channel = 'E'
      length1 = 3.5
      lengthz = 1
      
    endif else if chn eq 'b' then begin
      
      get_timespan, tget                            
      spp_fld_load, trange=tget, type='dfb_dbm_scm',/no_staging

      get_data, 'psp_fld_l2_dfb_dbm_scmlgu', data = scm_xlg
      get_data, 'psp_fld_l2_dfb_dbm_scmlgv', data = scm_ylg
      get_data, 'psp_fld_l2_dfb_dbm_scmlgw', data = scm_zlg

      get_data, 'psp_fld_l2_dfb_dbm_scmhgu', data = scm_xhg
      get_data, 'psp_fld_l2_dfb_dbm_scmhgv', data = scm_yhg
      get_data, 'psp_fld_l2_dfb_dbm_scmhgw', data = scm_zhg

      checklg = ISA(scm_xlg, 'INT')
      checkhg = ISA(scm_xhg, 'INT')

      if checklg eq 0 and checkhg eq 0 then begin

        dbm_dat_1 = { x: [scm_xlg.x,scm_xhg.x],$
          y: [scm_xlg.y,scm_xhg.y],$
          v: [scm_xlg.v,scm_xhg.v]}

        dbm_dat_2 = { x: [scm_ylg.x,scm_yhg.x],$
          y: [scm_ylg.y,scm_yhg.y],$
          v: [scm_ylg.v,scm_yhg.v]}

        dbm_dat_3 = { x: [scm_zlg.x,scm_zhg.x],$
          y: [scm_zlg.y,scm_zhg.y],$
          v: [scm_zlg.v,scm_zhg.v]}

      endif else if checklg eq 0 and checkhg eq 1 then begin

        dbm_dat_1 = scm_xlg
        dbm_dat_2 = scm_ylg
        dbm_dat_3 = scm_zlg

      endif else if checklg eq 1 and checkhg eq 0 then begin

        dbm_dat_1 = scm_xhg
        dbm_dat_2 = scm_yhg
        dbm_dat_3 = scm_zhg

      endif else begin

        get_data, 'psp_fld_l2_dfb_dbm_scmu', data = dbm_dat_1
        get_data, 'psp_fld_l2_dfb_dbm_scmv', data = dbm_dat_2
        get_data, 'psp_fld_l2_dfb_dbm_scmw', data = dbm_dat_3

      endelse

;      get_data, 'psp_fld_l2_dfb_dbm_scmu', data = dbm_dat_1
;      get_data, 'psp_fld_l2_dfb_dbm_scmv', data = dbm_dat_2
;      get_data, 'psp_fld_l2_dfb_dbm_scmw', data = dbm_dat_3
      
      fieldtype = 'Magnetic Field (nT) '
      color = 'goldenrod'
      channel = 'B'
      length1 = 1
      lengthz = 1

      ;length
    endif  
    
    spp_fld_load,trange=tget,type= 'mag_SC',/no_staging
    
    get_data, 'psp_fld_l2_mag_SC', data=mago
    
    ;help, mago
    
    magx = reform(mago.y[*,0])
    magy = reform(mago.y[*,1])
    magz = reform(mago.y[*,2])
    
    mag = { x: mago.x,$
            y: sqrt(magx^2+magy^2+magz^2)} 
    
  endif
  

  
  tplot_clear, /all

  if lvl eq 1 then begin
    
    xwavecheck = dbm_dat_1.x
    ywavecheck = dbm_dat_2.x
    zwavecheck = dbm_dat_3.x

    nwavex = n_elements(xwavecheck)
    nwavey = n_elements(ywavecheck)
    nwavez = n_elements(zwavecheck)

    wavecheck = [xwavecheck,ywavecheck,zwavecheck]

    uniqwavecheck = wavecheck[uniq(wavecheck, SORT(wavecheck))]

    n_uniqwaves = n_elements(uniqwavecheck)


    newx = fltarr(n_uniqwaves,524288)
    newy = fltarr(n_uniqwaves,524288)
    newz = fltarr(n_uniqwaves,524288)

    for lol=0,n_uniqwaves-1 do begin
      
      tmpwherex = where(dbm_dat_1.x eq uniqwavecheck[lol])
      tmpwherey = where(dbm_dat_2.x eq uniqwavecheck[lol])
      tmpwherez = where(dbm_dat_3.x eq uniqwavecheck[lol])
      
      if tmpwherex[0] ne -1 then begin
        newx[lol,*] = dbm_dat_1.y[tmpwherex[0],*]
      endif else newx[lol,*] = 0.
      
      if tmpwherey[0] ne -1 then begin
        newy[lol,*] = dbm_dat_2.y[tmpwherey[0],*]
      endif else newy[lol,*] = 0.
      
      if tmpwherez[0] ne -1 then begin
        newz[lol,*] = dbm_dat_3.y[tmpwherez[0],*]
      endif else newz[lol,*] =0
      


    endfor

;    help, uniqwavecheck
;    help, dbm_dat_1.x

    datetimes = uniqwavecheck

    truex = newx
    truey = newy
    truez = newz
    
    
  endif else begin
    datetimes = dbm_dat_1.x
    
    truex = dbm_dat_1.y
    truey = dbm_dat_2.y
    truez = dbm_dat_3.y
  endelse

  
  n_bursts = n_elements(datetimes)
  wavetime = dblarr(n_bursts, 524288)

  for n=0,n_bursts-1 do begin
    t1 = datetimes[n]
    times = dindgen(524288)/sr + t1
    wavetime[n,*] = times[*]
  endfor
  
  
  
  dbm_data = { starttimes: datetimes,$
    times: wavetime,$
    dataV12: truex,$
    dataV34: truey,$
    dataV5: truez,$
    magtime: mag.x,$
    mag: mag.y,$
    magx: magx, $
    magy: magy, $
    magz: magz}
    
  tplot_clear, /all
  
  
  tmp_where = where(dbm_data.starttimes ge date_dbl)
  
  
  i=tmp_where[0]
  
  
  ;-------------------------create plots-------------------------------;
  
  
  
  datewhere = where(dbm_data.starttimes eq datesave)
    
  if jhnycash gt 0 then if e eq -1 then i=n_bursts-1 else i=0
  if datewhere[0] ne -1 then i=datewhere[0]
  
;               ; Library of Indexes
;               ;
;  i=0          ; i = master index, this number controls flipping through waveforms
;               ; line is commented because i is handled elsewhere
  check = ''
  datesave = -1 ;
  j=0           ; j = menu index, this number controls when the menu is displayed
  ;k=0          ; k = doesnt appear to do anything anymore
  s=0           ; s = fft display, this number toggles 
  d=0           ; d = toggles decibel/linear scale for colorbar
  e=0           ; e = handles continuous transition between days
  l=0           ; l = toggles interpolation of waveform data
  p12=0         ; p12 = toggles waveform zoom.
  p13=0         ; p13 = handles waveform zoom window memory. 
  num=1         ; num = changes the scale of the FFT y axis.
  kpop=0        ; kpop = bad music, but also iterates between instrument coordinate system and spacecraft coordinate system.
  
  tag_int = 0
  tag = ''
  
  repeat begin
    
    if l mod 2 eq 0 then time = dindgen(524288)/sr ;else time = dindgen(52428*2)/(sr/5)  ; time for one capture
    
    
    if i lt 0 then begin
      t = time_string(time_double(t)-21600)
      timespan, t, 0.25
      e = -1
      goto, move_day
    endif else if i gt n_bursts-1 then begin
      t = time_string(time_double(t)+21600)
      timespan, t, 0.25
      e = 1
      goto, move_day
    endif
    
    
    
    ;print, i
    ;a = string(i)
    magdif = abs(dbm_data.magtime - dbm_data.starttimes[i])
    ;print, magdif
    a = where(magdif eq min(magdif))
    
    bfield = dbm_data.mag[a] ; in nT
    bfield_tesla = bfield*10.^(-9) ; in Tesla
    
    fce = (1./(2*pi))*(q*bfield_tesla)/mass_e ;should be in Hz
    fce = fce[0]
    
    fcestr = string(fce)
    fcestr = strtrim(fcestr,1)
    
    if num eq 0 then begin
      fce = fce
    endif else if num eq 1 then begin
      fce = 2.*fce
    endif else if num eq 2 then begin
      fce = 4.*fce
    endif else if num eq 3 then begin
      fce = 6.*fce
    endif else if num eq 4 then begin
      fce = 8.*fce
    endif else if num eq 5 then begin
      
      
      if check eq '' then begin
        print, ""
        READ, check, PROMPT='Enter custom fce multiple: '
        print, "'SHIFT'+'0' to reset custom scale."
        print, ""
        check = float(check)
      endif

      
      
      
      fce = check*fce
    endif
    
    ;if chn eq 1 or chn eq 2 then length = 5.
    ;if chn eq 3 then length = 6. else length = 1.
    
    ;length1 = 5.
    ;length2 = 6.
    
    waveduration = reform(dbm_data.times[i,*]-dbm_data.starttimes[i])
    
    if chn eq 'e' then begin
      wavedatax = dbm_data.dataV12[i,*]
      wavedatay = dbm_data.dataV34[i,*]
      wavedataz = dbm_data.dataV5[i,*]
      
      wavedatax = reform(wavedatax)
      wavedatay = reform(wavedatay)
      wavedataz = reform(wavedataz)
      
      wavedatax = wavedatax/length1  ;
      wavedatay = wavedatay/length1  ; Effective Length calculations
      wavedataz = wavedataz/lengthz  ;
      
      wavedatax = wavedatax - mean(wavedatax)
      wavedatay = wavedatay - mean(wavedatay)
      wavedataz = wavedataz - mean(wavedataz)
      
      ;PSP L2 Ew files need the -1 sign correction. 
      ;This should be fixed in L3.
      if lvl eq 2 then begin
        wavedatax = -1*wavedatax*1000.
        wavedatay = -1*wavedatay*1000. ; convert from V/m to mV/m?
        wavedataz = -1*wavedataz*1000.
      endif

      
      ;rotate to spacecraft coordinates;
      if kpop mod 2 eq 0 then begin

        ;Rotating V12, V34 into SC coord (x,y)
        rot_mat = [[0.64524,-0.82228,0.],$
                  [0.76897,0.57577,0.],$
                  [0,0,1]]

        wavedataSENSOR = [[wavedatax],[wavedatay],[wavedataz]]
        wavedataSC = reform(rot_mat ## wavedataSENSOR)

        wavedatax = wavedataSC[*,0]
        wavedatay = wavedataSC[*,1]
        wavedataz = wavedataSC[*,2]

        coordsys = 'SC'
      endif else coordsys = 'IC'

      if p12 mod 2 ne 0 then begin
        
        try_againE:
        
        if p13 eq 0 then begin
           zoomstart = ''
           zoomend = ''
        endif
        
        ;if b eq 'p' or b eq 'o' then zoomstart = zoomstarttmp & zoomend = zoomendtmp      

        if zoomstart eq '' then READ, zoomstart, PROMPT='Enter zoom window start time: '
        if zoomend eq '' then begin
          READ, zoomend, PROMPT='Enter zoom window end time: '
          print, ""
        endif
        
        
        zoomstarttmp = zoomstart
        zoomendtmp = zoomend
        
        startflt = float(zoomstart)
        endflt = float(zoomend)
        
        if endflt le startflt then begin
          print, ""
          print, "Endtime cannot be before starttime, try different window."
          print, ""
          zoomstart = ''
          zoomend = ''
          wait,2
          goto,try_againE
        endif
        
        zoomwhere = where(waveduration ge startflt and waveduration le endflt)
        
        if zoomwhere[0] eq -1 then begin
          print, ""
          print, "Invalid window, try different window."
          print, ""
          zoomstart = ''
          zoomend = ''
          wait, 2
          goto, try_againE
        endif
        
        wavedatax = wavedatax[zoomwhere]
        wavedatay = wavedatay[zoomwhere]
        wavedataz = wavedataz[zoomwhere]
        waveduration = waveduration[zoomwhere]
        
        time = dindgen(n_elements(wavedatax))/sr; + startflt
        
        
      endif

      
      ;--------------------------------;
      
    endif else if chn eq 'b' then begin
      
      
      wavedatax = dbm_data.dataV12[i,*]
      wavedatay = dbm_data.dataV34[i,*]
      wavedataz = dbm_data.dataV5[i,*]
      ;-------------------------------------------------------;
      wavedatax = reform(wavedatax)
      wavedatay = reform(wavedatay)
      wavedataz = reform(wavedataz)
      
      wavedatax = wavedatax - mean(wavedatax)
      wavedatay = wavedatay - mean(wavedatay)
      wavedataz = wavedataz - mean(wavedataz)
      
      if kpop mod 2 eq 0 then begin

        rot_mat_scm = [[0.81654,-0.40827,-0.40827],$
                        [0.,-0.70715,0.70715],$
                        [-0.57729,-0.57729,-0.57729]]

        wavedataSENSOR = [[wavedatax],[wavedatay],[wavedataz]]
        wavedataSC = reform(rot_mat_scm ## wavedataSENSOR)

        wavedatax = wavedataSC[*,0]
        wavedatay = wavedataSC[*,1]
        wavedataz = wavedataSC[*,2]
        
        coordsys = 'SC'
        
      endif else coordsys = 'IC'
      
      if p12 mod 2 ne 0 then begin

        try_againB:

        if p13 eq 0 then begin
           zoomstart = ''
           zoomend = ''
        endif
        ;if b eq 'p' or b eq 'o' then zoomstart = zoomstarttmp & zoomend = zoomendtmp

        if zoomstart eq '' then READ, zoomstart, PROMPT='Enter zoom window start time: '
        if zoomend eq '' then begin
          READ, zoomend, PROMPT='Enter zoom window end time: '
          print, ""
        endif

        zoomstarttmp = zoomstart
        zoomendtmp = zoomend

        startflt = float(zoomstart)
        endflt = float(zoomend)

        if endflt le startflt then begin
          print, ""
          print, "Endtime cannot be before starttime, try different window."
          print, ""
          zoomstart = ''
          zoomend = ''
          wait,2
          goto,try_againB
        endif

        zoomwhere = where(waveduration ge startflt and waveduration le endflt)

        if zoomwhere[0] eq -1 then begin
          print, ""
          print, "Invalid window, try different window."
          print, ""
          zoomstart = ''
          zoomend = ''
          wait, 2
          goto, try_againB
        endif

        wavedatax = wavedatax[zoomwhere]
        wavedatay = wavedatay[zoomwhere]
        wavedataz = wavedataz[zoomwhere]
        waveduration = waveduration[zoomwhere]

        time = dindgen(n_elements(wavedatax))/sr; + startflt

      endif

      
    endif
    
;    help, wavedatax, wavedatay, wavedataz

    starttime = dbm_data.starttimes[i]
    
    srt = 150000.
    lf = 20.    ;Hz
    hf = 1000.  ;Hz
    
    if l mod 2 ne 0 then begin
      ;help, wavedatax
      
      x_tmp = [[wavedatax],[wavedatay],[wavedataz]]
      
      ;help, x_tmp
      
      x = float(rbsp_vector_bandpass(x_tmp,srt,lf,hf))
      
      ;help,x
      
      wavedatax = x[*,0]
      wavedatay = x[*,1]
      wavedataz = x[*,2]
      
      time = dindgen(n_elements(wavedatax))/sr
      
      ;print, n_elements(wavedatax), ' TEST'
      
    endif
    
    
    if s mod 2 eq 0 then begin
      acol = 2
      bcolx = 1
      bcoly = 3
      bcolz = 5
    endif else begin
      acol = 1
      bcolx = 1
      bcoly = 2
      bcolz = 3
    endelse
    
    n_dur = n_elements(waveduration)
    
    px = plot(waveduration, wavedatax,$
       TITLE=time_string(starttime)+' '+coordsys,DIM=[1500,750], XRANGE=[waveduration[0],waveduration[n_dur-1]],$
       YTITLE=fieldtype+'    '+channel+'x', color=color, XTITLE='Time (Seconds)', LAYOUT=[acol,3,bcolx])
       
    py = plot(waveduration, wavedatay,$
       TITLE=time_string(starttime)+' '+coordsys,DIM=[1500,750], XRANGE=[waveduration[0],waveduration[n_dur-1]],$
       YTITLE=fieldtype+'    '+channel+'y', color=color, XTITLE='Time (Seconds)', LAYOUT=[acol,3,bcoly],/CURRENT)
       
    pz = plot(waveduration, wavedataz,$
       TITLE=time_string(starttime)+' '+coordsys,DIM=[1500,750], XRANGE=[waveduration[0],waveduration[n_dur-1]],$
       YTITLE=fieldtype+'    '+channel+'z', color=color, XTITLE='Time (Seconds)', LAYOUT=[acol,3,bcolz],/CURRENT)
       
    ;-----------------------------create sliding FFT/Regular FFT----------------------------------;   
       
    if s mod 2 eq 0 then begin
      ct = colortable(74, /REVERSE)
    
      ;w = window(DIMENSIONS=[950,400])
      zoomlength = (waveduration[n_dur-1] - waveduration[0]) 
       
      if d mod 2 eq 0 then begin
        powerspecx = slide_spec(time,wavedatax,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs,/db)
        powerspecy = slide_spec(time,wavedatay,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs,/db)
        powerspecz = slide_spec(time,wavedataz,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs,/db)
        
        powerspecx=powerspecx[2:40,*]
        powerspecy=powerspecy[2:40,*]
        powerspecz=powerspecz[2:40,*]
        
        xfft = 10*alog10(2*(abs(fft(wavedatax)))^2) ;10*alog10(2*padded_step_length*(abs(fft(wavedatax)))^2)
        yfft = 10*alog10(2*(abs(fft(wavedatay)))^2)
        zfft = 10*alog10(2*(abs(fft(wavedataz)))^2)
        
        title_tmp = 'dB'
        
      endif else begin
        powerspecx = slide_spec(time,wavedatax,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs)
        powerspecy = slide_spec(time,wavedatay,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs)
        powerspecz = slide_spec(time,wavedataz,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs)
        
        powerspecx=powerspecx[2:40,*]
        powerspecy=powerspecy[2:40,*]
        powerspecz=powerspecz[2:40,*]
        
        xfft = 2*(abs(fft(wavedatax)))^2
        yfft = 2*(abs(fft(wavedatay)))^2
        zfft = 2*(abs(fft(wavedataz)))^2
        
        title_tmp = 'Linear Units'
      endelse
      
      
      
      if zoomlength ge 0.05 then begin
        
        freqs = freqs*1000. ; in Hz
        
        a1 = where((freqs gt 15) and (freqs le fce))
        nanwhere = where(powerspecx eq min(powerspecx))
        powerspecx[nanwhere]=!values.f_nan

        b1 = where((freqs gt 15) and (freqs le fce))
        nanwhere = where(powerspecy eq min(powerspecy))
        powerspecy[nanwhere]=!values.f_nan

        c1 = where((freqs gt 15) and (freqs le fce))
        nanwhere = where(powerspecz eq min(powerspecz))
        powerspecz[nanwhere]=!values.f_nan

        ;print, a
        ;help, powerspec[*,a]
        ;print, freqs[a]

        gx = image(powerspecx[*,a1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
          AXIS_STYLE=4, LAYOUT=[2,3,2], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]

        gy = image(powerspecy[*,b1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
          AXIS_STYLE=4, LAYOUT=[2,3,4], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]

        gz = image(powerspecz[*,c1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
          AXIS_STYLE=4, LAYOUT=[2,3,6], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]

        clrbar1 = [total(wavedatax),total(wavedatay),total(wavedataz)]
        ;print, clrbar1

        clrbar2 = [gx,gy,gz]

        clrbarwhere = where(abs(clrbar1) gt 0)


;.compile /Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/mission_routines/solar_probe_routines/psp_waves/colorbar.pro

        if d mod 2 eq 0 then begin
          c = colorbar_psp_waves(TARGET=clrbar2[clrbarwhere[0]], ORIENTATION=1, POSITION=[.945,.1,.96,.9], TITLE='Decibels (dB)',TEXTPOS=1, BORDER=1, TICKDIR=1)
        endif else c = colorbar_psp_waves(TARGET=clrbar2[clrbarwhere[0]], ORIENTATION=1, POSITION=[.945,.1,.96,.9], TITLE='Linear Units',TEXTPOS=1, BORDER=1, TICKDIR=1)

        yaxisx = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gx)
        xaxisx = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gx)

        yaxisy = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gy)
        xaxisy = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gy)

        yaxisz = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gz)
        xaxisz = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gz)
        
      endif else begin
;        
;        help, xfft
;        help, wavedatax
        
        nyq_freq = 150000./2.; Hz
        
        freqs = (dindgen(n_dur)/n_dur)*nyq_freq
        
        fft_where = where(freqs lt fce)
        
        gx = plot(freqs,xfft,LAYOUT=[2,3,2],TITLE='FFT'+'   fce='+fcestr, /current, XRANGE=[15,fce],YTITLE=title_tmp)
        gy = plot(freqs,yfft,LAYOUT=[2,3,4],TITLE='FFT'+'   fce='+fcestr, /current, XRANGE=[15,fce],YTITLE=title_tmp)
        gz = plot(freqs,zfft,LAYOUT=[2,3,6],TITLE='FFT'+'   fce='+fcestr, /current, XRANGE=[15,fce],YTITLE=title_tmp)
        
      endelse
      
;      a1 = where((freqs gt 15) and (freqs le fce))
;      nanwhere = where(powerspecx eq min(powerspecx))
;      powerspecx[nanwhere]=!values.f_nan
;      
;      b1 = where((freqs gt 15) and (freqs le fce))
;      nanwhere = where(powerspecy eq min(powerspecy))
;      powerspecy[nanwhere]=!values.f_nan
;      
;      c1 = where((freqs gt 15) and (freqs le fce))
;      nanwhere = where(powerspecz eq min(powerspecz))
;      powerspecz[nanwhere]=!values.f_nan
;      
;      ;print, a
;      ;help, powerspec[*,a]
;      ;print, freqs[a]
;      
;      gx = image(powerspecx[*,a1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
;         AXIS_STYLE=4, LAYOUT=[2,3,2], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
;         
;      gy = image(powerspecy[*,b1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
;         AXIS_STYLE=4, LAYOUT=[2,3,4], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
;         
;      gz = image(powerspecz[*,c1], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
;         AXIS_STYLE=4, LAYOUT=[2,3,6], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
;       
;      clrbar1 = [total(wavedatax),total(wavedatay),total(wavedataz)]
;      ;print, clrbar1
;      
;      clrbar2 = [gx,gy,gz]
;      
;      clrbarwhere = where(abs(clrbar1) gt 0)
;       
;      if d mod 2 eq 0 then begin 
;       c = colorbar(TARGET=clrbar2[clrbarwhere[0]], ORIENTATION=1, POSITION=[.945,.1,.96,.9], TITLE='Decibels (dB)',TEXTPOS=1, BORDER=1, TICKDIR=1)
;      endif else c = colorbar(TARGET=clrbar2[clrbarwhere[0]], ORIENTATION=1, POSITION=[.945,.1,.96,.9], TITLE='Linear Units',TEXTPOS=1, BORDER=1, TICKDIR=1)
;       
;      
;      
;      yaxisx = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gx)
;      xaxisx = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gx)
;      
;      yaxisy = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gy)
;      xaxisy = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gy)
;      
;      yaxisz = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gz)
;      xaxisz = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[waveduration[0],zoomlength/1750.], TARGET=gz)
;       
       
       
    endif

    ;-------------------------------------------------------------------------------------------;
       
    ;endif
    ;if j eq 0 then tlimit,0,0
    ;if j mod 10 eq 0 then print, "'p' for next plot. 'o' for previous. 'r' to reset plot. 'Enter' for ironic exit."
    if j mod 10 eq 0 then begin
       print, ""
       print, "'p' for next plot. 'o' for previous."
       print, "'s' to save plots. 'f' to hide/unhide FFT."
       print, "'d' to toggle dB color scale."
       print, "'1','2','4','6', or '8' to change freq scale in multiples of fce."
       print, "'0' to enter custom freq scale. (In multiples of fce.)"
       print, "'SHIFT'+'1','2' to change field type to E, B respectively."
       print, "'g' to hide/unhide zoom. 'h' to reset zoom window."
       print, "'k' to toggle SC/Instrument coordinates."
       print, "'l' to apply bandpass filter to data."
       print, ""
       print, "Press 'q' to exit."
       print, ""
    endif
    
    b = get_kbrd()
    ;if b eq 'r' then tlimit,0,0
    
    ;------------scroll thru captures---------;
    
    if b eq 'p' then begin 
      i+=1
      tag_int = 0
      tag=''
    endif  
    
    if b eq 'o' then begin 
      i-=1
      tag_int = 0
      tag=''
    endif
    
    j+=1
    
    ;-----------adjust freq scale-------------; 
    
    if b eq '1' then num=0 ; set maximum to fce
    if b eq '2' then num=1 ; max to 2*fce
    if b eq '4' then num=2 ; max to 4*fce
    if b eq '6' then num=3 ; max to 6*fce
    if b eq '8' then num=4 ; max to 8*fce
    if b eq '0' then num=5 ; max to 20*fce
    
    ;-----------------save--------------------;
    

    
    if b eq 's' then begin
      
      date = strmid(time_string(dbm_data.starttimes[i]),0,4)+$
        strmid(time_string(dbm_data.starttimes[i]),5,2)+$
        strmid(time_string(dbm_data.starttimes[i]),8,2)
      
      time0 = strmid(time_string(dbm_data.starttimes[i]),11,2)+$
        strmid(time_string(dbm_data.starttimes[i]),14,2)+$
        strmid(time_string(dbm_data.starttimes[i]),17,2)
      

      
      datetime = date+'_'+time0  
      
      if s mod 2 eq 0 then datetime = datetime+'_fft'
      
      ;print, datetime
      if ~file_test('~/Desktop/PSP_Whist_Plots/') then begin
        
        file_mkdir,'~/Desktop/PSP_Whist_Plots/'
        print, 'Created new directory ~/Desktop/PSP_Whist_Plots/'
        print, 'Plots are saved there.'
        
        
      endif
      
      recheck:
      
      if file_test('~/Desktop/PSP_Whist_Plots/'+channel+'_'+datetime+tag+'.png') then begin
        tag_int+=1
        tag = '_'+strtrim(string(tag_int),1)
        goto, recheck
      endif
      
      print,""
      print, "Saving..."
      
      px.save, '~/Desktop/PSP_Whist_Plots/'+channel+'_'+datetime+tag+'.png'
      
      print,""
      print,'Saved as: '+channel+'_'+datetime+tag+'.png'
      print,""
      ;if k mod 2 ne 0 then g.save,'~/Desktop/PSP_Whist_Plots/'+channel+'_'+datetime+'_sfft.png'
      
    endif
    
    px.close
    
    ;if k mod 2 ne 0 then gx.close
    if p12 mod 2 ne 0 then p13+=1
    ;------------change channel---------------;

    if b eq '!' then begin
      chn='e'
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '@' then begin
      chn='b'
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    
    ;---------------FFT options menu---------------;
    ;if b eq 'f' then k+=1 ;toggle fft menu
    if b eq 'f' then s+=1 ;toggle sliding fft
    if b eq 'd' then d+=1 ;toggle dB colorbar scale
    
    ;------------------less data-------------------;   
    
    if b eq 'g' then p12+=1
    
    if p12 mod 2 ne 0 and b eq 'h' then p13=0
    
    if b eq 'l' then l+=1
    
    if b eq 'k' then kpop+=1
    
    if b eq ')' then check=''
    
    jhnycash+=1
    
  endrep until b eq 'q' ;press q to quit
  
  
  
  
end