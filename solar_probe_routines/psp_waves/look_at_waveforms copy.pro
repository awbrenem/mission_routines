;+
    ; :Author: benshort
    ;-
    ;+
    ; :Description:
    ;    A procedure created in order to scroll through
    ;    Parker Solar Probe 3.5 second captures.
    ;     
    ;    the CHN parameter allows you to choose which antenna to look at.
    ;    
    ;
    ; :Params:
    ;    t0
    ;    chn
    ;
    ;    
    ; :Dependencies:
    ;    tplot_clear.pro
    ;    colorbar.pro
    ;
    ; :Author: benshort
    ;-
pro look_at_waveforms, t0, chn
  
  if getenv('PSP_STAGING_ID') EQ '' then setenv, 'PSP_STAGING_ID=short186'
  if getenv('PSP_STAGING_PW') EQ '' then setenv, 'PSP_STAGING_PW=poLuSHISheoc'

  timespan, t0
  t = t0
  
  Bt = 3.49525 ;burst length (sec)
  sr = 150000. ;samples/sec
  pi = 3.14159
  q = 1.60217662*10.^(-19) ; coulombs
  mass_e = 9.10938356*10.^(-31) ;kg
  
  e=0
  i=0
  datesave = -1
  
  move_day:
  
  
  spp_fld_make_or_retrieve_cdf,'mago_survey', /load
  get_data, 'spp_fld_mago_survey_nT_mag', data=mag
  
  tplot_clear, /all
  
  if chn eq 1 then begin
    channel = 'V12'
    fieldtype = 'Electric Field (V/m) '
    color = 'blue'
  endif else if chn eq 2 then begin
    channel = 'V34' 
    fieldtype = 'Electric Field (V/m) ' 
    color = 'blue'
  endif else if chn eq 3 then begin
    channel = 'V5' 
    fieldtype = 'Electric Field (V/m) ' 
    color = 'blue'
  endif else if chn eq 4 then begin
    channel = 'SCM Bx'
    fieldtype = 'Magnetic Field (nT) ' 
    color = 'goldenrod'
  endif else if chn eq 5 then begin
    channel = 'SCM By'
    fieldtype = 'Magnetic Field (nT) '
    color = 'goldenrod'
  endif else if chn eq 6 then begin
    channel = 'SCM Bz'
    fieldtype = 'Magnetic Field (nT)'
    color = 'goldenrod'
  endif
  
  dbm_str = 'dfb_dbm_'+strtrim(string(chn),1)
  dbm_chn_str = 'spp_fld_dfb_dbm_'+strtrim(string(chn),1)+'_dbm_data'
  
  ;print, channel, fieldtype, color
  
  spp_fld_make_or_retrieve_cdf,dbm_str,/load   ;waveform bursts
  get_data, dbm_chn_str, data = dbm_dat
  
  n_bursts = n_elements(dbm_dat.x)
  wavetime = dblarr(n_bursts, 524288)
  
  
  
  for n=0,n_bursts-1 do begin
    t1 = dbm_dat.x[n]
    times = dindgen(524288)/sr + t1
    wavetime[n,*] = times[*]
  endfor
  
  dbm_data = { starttimes: dbm_dat.x,$
                     times: wavetime,$
                     data: dbm_dat.y,$
                     magtime: mag.x,$
                     mag: mag.y}
  
  tplot_clear, /all
  
  datewhere = where(dbm_data.starttimes eq datesave)
  
  
  ;-------------------------create plots-------------------------------;
  
  if e eq -1 then i=n_bursts-1 else i=0
  if datewhere[0] ne -1 then i=datewhere[0]
  
  datesave = -1
  j=0
  k=0
  s=0
  d=0
  e=0
  l=0
  num=1
  repeat begin
    
    if l mod 2 eq 0 then time = dindgen(524288)/sr else time = dindgen(52428*2)/(sr/5)  ; time for one capture
    
    
    if i lt 0 then begin
      t = time_string(time_double(t)-86400)
      timespan, t
      e = -1
      goto, move_day
    endif else if i gt n_bursts-1 then begin
      t = time_string(time_double(t)+86400)
      timespan, t
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
      fce = 20.*fce
    endif
    
    if chn eq 1 or chn eq 2 then length = 5.
    if chn eq 3 then length = 6. else length = 1.
    
    waveduration = dbm_data.times[i,*]-dbm_data.starttimes[i]
    wavedata = dbm_data.data[i,*]/length
    wavedata = wavedata - mean(dbm_data.data[i,*]/length)
    starttime = dbm_data.starttimes[i]
    
    if l mod 2 ne 0 then wavedata = interpol(wavedata, 52428*2)
    
    px = plot(waveduration, wavedata,$
       TITLE=time_string(starttime),DIM=[1500,750], XRANGE=[0,3.495],$
       YTITLE=fieldtype+'    '+channel, color=color, XTITLE='Time (Seconds)', LAYOUT=[2,3,1])
       
    py = plot(waveduration, wavedata,$
       TITLE=time_string(starttime),DIM=[1500,750], XRANGE=[0,3.495],$
       YTITLE=fieldtype+'    '+channel, color=color, XTITLE='Time (Seconds)', LAYOUT=[2,3,3],/CURRENT)
       
    pz = plot(waveduration, wavedata,$
       TITLE=time_string(starttime),DIM=[1500,750], XRANGE=[0,3.495],$
       YTITLE=fieldtype+'    '+channel, color=color, XTITLE='Time (Seconds)', LAYOUT=[2,3,5],/CURRENT)
       
       
    ;if k mod 2 ne 0 then begin
       
    ;-----------------------------create sliding FFT/Regular FFT----------------------------------;   
       
    if s mod 2 eq 0 then begin
      ct = colortable(74, /REVERSE)
    
      ;w = window(DIMENSIONS=[950,400])
       
      if d mod 2 eq 0 then begin
        powerspec = slide_spec(time,wavedata,0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs,/db)
      endif else powerspec = slide_spec(time,dbm_data.data[i,*],0.05,0.5,iwindow=1,time_index=time_index,freq_bins=freqs)
      
      ;help, powerspec
      
      powerspec=powerspec[2:40,*]
      
      freqs = freqs*1000. ; in Hz
      
      a = where((freqs gt 15) and (freqs le fce))
      nanwhere = where(powerspec eq min(powerspec))
      powerspec[nanwhere]=!values.f_nan
       ;print, a
      ;help, powerspec[*,a]
       
       
      ;print, freqs[a]
      
      gx = image(powerspec[*,a], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
         AXIS_STYLE=4, LAYOUT=[2,3,2], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
         
         
      gy = image(powerspec[*,a], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
         AXIS_STYLE=4, LAYOUT=[2,3,4], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
         
      gz = image(powerspec[*,a], RGB_TABLE=ct,IMAGE_DIMENSIONS=[1750,400],TITLE='Sliding FFT'+'   fce='+fcestr,$
         AXIS_STYLE=4, LAYOUT=[2,3,6], MARGIN=0.13,/CURRENT) ;time_index,freqs[1:2000], ASPECT_RATIO=[200,900]
       
      if d mod 2 eq 0 then begin 
       c = colorbar(TARGET=gy, ORIENTATION=1, POSITION='right', TITLE='Decibels (dB)',TEXTPOS=1, BORDER=1, TICKDIR=1)
      endif else c = colorbar(TARGET=gx, ORIENTATION=1, POSITION='right', TITLE='Linear Units',TEXTPOS=1, BORDER=1, TICKDIR=1)
       
       
      yaxisx = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gx)
      xaxisx = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[0,3.495/1750.], TARGET=gx)
      
      yaxisy = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gy)
      xaxisy = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[0,3.495/1750.], TARGET=gy)
      
      yaxisz = AXIS('Y', LOCATION='left', TITLE='Freq (Hz)', COORD_TRANSFORM=[0, fce/400.], TARGET=gz)
      xaxisz = AXIS('X',LOCATION='bottom', TITLE='Time (Seconds)', COORD_TRANSFORM=[0,3.495/1750.], TARGET=gz)
       
       
       
    endif else begin
      
      h = fft(dbm_data.data[i,*])
      g = plot(h, YRANGE=[0,0.5], LAYOUT=[1,2,2],TITLE='FFT'+'   fce='+fcestr+'(Hz)');,/CURRENT)
    endelse
    
    ;-------------------------------------------------------------------------------------------;
       
    ;endif
    ;if j eq 0 then tlimit,0,0
    ;if j mod 10 eq 0 then print, "'p' for next plot. 'o' for previous. 'r' to reset plot. 'Enter' for ironic exit."
    if j eq 0 then begin
       print, "'p' for next plot. 'o' for previous. 'Enter' to exit."
       print, "'s' to save plots. 'f' to toggle sliding fft/normal fft."
       print, "'d' to toggle dB color scale."
       print, "'1','2','4','6', or '8' to change freq scale in multiples of fce."
       print, "'0' to change freq scale to 20*fce"
       print, "'SHIFT'+'1','2','3','4','5', or '6' to change channel."
       print, "'l' to interpolate data"
    endif
    
    b = get_kbrd()
    ;if b eq 'r' then tlimit,0,0
    
    ;------------scroll thru captures---------;
    
    if b eq 'p' then i+=1
    if b eq 'o' then i-=1
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
      ;print, datetime
      if ~file_test('~/Desktop/PSP_Whist_Plots/') then begin
        
        file_mkdir,'~/Desktop/PSP_Whist_Plots/'
        print, 'Created new directory ~/Desktop/PSP_Whist_Plots/'
        print, 'Plots are saved there.'
      endif
      
      p.save, '~/Desktop/PSP_Whist_Plots/'+channel+'_'+datetime+'.png'
      ;if k mod 2 ne 0 then g.save,'~/Desktop/PSP_Whist_Plots/'+channel+'_'+datetime+'_sfft.png'
      
    endif
    px.close
    if k mod 2 ne 0 then gx.close
    
    ;------------change channel---------------;

    if b eq '!' then begin
      chn=1
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '@' then begin
      chn=2
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '#' then begin
      chn=3
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '$' then begin
      chn=4
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '%' then begin
      chn=5
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    if b eq '^' then begin
      chn=6
      datesave = dbm_data.starttimes[i]
      goto, move_day
    endif
    
    ;---------------FFT options menu---------------;
    ;if b eq 'f' then k+=1 ;toggle fft menu
    if b eq 'f' then s+=1 ;toggle sliding fft
    if b eq 'd' then d+=1 ;toggle dB colorbar scale
    
    ;------------------less data-------------------;   
    
    if b eq 'l' then l+=1
    
  endrep until byte(b) eq 10 ;press enter
  
  
  
  
end