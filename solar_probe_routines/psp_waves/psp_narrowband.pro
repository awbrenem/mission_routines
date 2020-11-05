;
;
;
;
;
;


FUNCTION new_ev_struct,keeptime=keeptime
  if keyword_set(keeptime) then begin
    ev =  {raw_data:dblarr(524288,3),$       ;time series data
      mv_data:dblarr(52428*2,3),$
      quality:fltarr(3,39),$         ;quality of event at each frequency
      freq:fltarr(3),$          ;highest power frequency
      max_amp:fltarr(3),$  ;maximum mV amplitude
      bw:fltarr(3),$   ;bandwidth
      df_f:fltarr(3),$   ;normalized bandwidth
      source:['--','--','--'],$     ;source being recorded
      drate:['--','--','--'],$      ;data rate
      npoints:replicate(!values.f_nan,3),$ ;number of samples
      stime:replicate(!values.f_nan,3),$   ;time value
      filter:['--','--','--'],$     ;filter
      id:0L,$                ;event id
      triggersource:['--','--','--'],$  ;source that triggered the recording
      tds_q:replicate(!values.f_nan,3),$   ;onboard quality
      deltat:0D,$              ;sample time
      scetstr:'',$             ;string event time 'YYYY-MM-DD/hh:mm:ss.fff'
      scetur8:keeptime,$           ;ur8 event time
      thrust:0.,$              ;thruster state
      cal_state:0,$              ;calibration state
      ev_err:0,$               ;event error
      prog_err:0,$              ;program error
      duration:fltarr(3),$      ;time length of maximum coherency
      f_fce: fltarr(3), $       ;f_fce for each channel
      norm_angle: fltarr(1), $   ; wave propogation angle
      coev: fltarr(1), $ ;costreaming electron energy resonances
      countev:fltarr(1) $ ;counterstreaming electron energy resonances
    }
  endif else begin

    ev =  {raw_data:dblarr(524288,3),$       ;time series data
      mv_data:dblarr(52428*2,3),$
      quality:fltarr(3,39),$       ;quality of event at each frequency
      freq:fltarr(3),$            ;highest power frequency
      max_amp:fltarr(3),$  ;maximum mV amplitude
      bw: fltarr(3),$   ;bandwidth
      df_f:fltarr(3),$ ;normalized bandwidth
      source:['--','--','--'],$     ;source being recorded
      drate:['--','--','--'],$      ;data rate
      npoints:replicate(!values.f_nan,3),$ ;number of samples
      stime:replicate(!values.f_nan,3),$   ;time value
      filter:['--','--','--'],$     ;filter
      id:0L,$                ;event id
      triggersource:['--','--','--'],$  ;source that triggered the recording
      tds_q:replicate(!values.f_nan,3),$   ;onboard quality
      deltat:0D,$              ;sample time
      scetstr:'',$             ;string event time 'YYYY-MM-DD/hh:mm:ss.fff'
      scetur8:0D,$             ;ur8 event time
      thrust:0.,$              ;thruster state
      cal_state:0,$              ;calibration state
      ev_err:0,$               ;event error
      prog_err:0,$              ;program error
      duration:fltarr(3),$      ;time length of maximum coherency
      f_fce: fltarr(3), $       ;f_fce for each channel
      norm_angle: fltarr(1), $   ; wave propogation angle
      coev: fltarr(1), $ ;costreaming electron energy resonances
      countev:fltarr(1) $ ;counterstreaming electron energy resonances
    }
  endelse
  return,ev
END   ;function new_ev_struct

pro psp_narrowband, t0, tf, type=type, savepath=savepath,$
  rotsc=rotsc,rotminvar=rotminvar,rotfa=rotfa

  if getenv('PSP_STAGING_DIR') EQ '' then setenv, 'PSP_STAGING_DIR=~/Users/benshort/data/PSP/'
  
  
  ;---------Sets username and password if not set already----------:
  if getenv('USER') EQ '' then begin
    
    a = ''
    print, 'Please enter your User ID:'
    repeat begin

      b = get_kbrd()
      a = a + b
    endrep until byte(b) eq 10 ;press enter
    g = strmid(a,0,strlen(a)-1)
    setenv, 'USER='+g  
  endif

  if getenv('PSP_STAGING_ID') EQ '' then setenv, 'PSP_STAGING_ID=short186'

  if getenv('PSP_STAGING_PW') EQ '' then setenv, 'PSP_STAGING_PW=poLuSHISheoc'

;  if getenv('PSP_STAGING_ID') EQ '' then begin
;    
;    a = ''
;    print, 'Please enter your Staging ID (probably your User ID):'
;    repeat begin
;      b = get_kbrd(/escape)
;      a = a + b
;    endrep until byte(b) eq 10 ;press enter
;    g = strmid(a,0,strlen(a)-1)
;    print, 'Thanks!' 
;    setenv, 'PSP_STAGING_ID='+g  
;  endif
;
;  if getenv('PSP_STAGING_PW') EQ '' then begin
;    
;    a = ''
;    print, 'Please enter your Staging Password:'
;    repeat begin
;      b = get_kbrd(/escape)
;      a = a + b
;    endrep until byte(b) eq 10 ;press enter
;    g = strmid(a,0,strlen(a)-1)
;    print, 'Thanks!'
;    setenv, 'PSP_STAGING_PW='+g
;  endif
  ;--------------------------end password setting----------------------------;
  
  ;---------------sets savepath for output file-------------;

  if ~keyword_set(savepath) then begin
    case !VERSION.OS of
      'Windows': $
        begin
        ;This is not tested, I'm only assuming it works...
        spawn,'echo %HOMEPATH%',home
      end
      ELSE: $;basically, UNIX
        begin
        spawn,'echo ~/',home
        savepath = home+'Desktop/PSP_Raw_Narrowband/'
        if ~file_test(savepath,/directory) then file_mkdir,savepath,/noexpand_path
      end
    endcase
  endif
;--------------------end savepath--------------------;

  if ~keyword_set(t0) then begin
  
    a = ''
    print, 'Please enter start date:'
    repeat begin
      b = get_kbrd()
      a = a + b
    endrep until byte(b) eq 10 ;press enter
    t0 = strmid(a,0,strlen(a)-1)
    
  endif
  
  if ~keyword_set(tf) then begin
  
    a = ''
    print, 'Please enter end date:'
    repeat begin
      b = get_kbrd()
      a = a + b
    endrep until byte(b) eq 10 ;press enter
    tf = strmid(a,0,strlen(a)-1)
  
  endif
  
  version = 'PSP Narrowband ID v0.0.1'
  print, version
  ;----------------------load data---------------------;
  
  if ~keyword_set(type) then type='spec'
  
  
  
  daynum = CEIL(abs(time_double(t0)-time_double(tf))/86400)
  
  init_crib_colors
  t = t0
  
  
  days=0
  read_quit=''
  counter=0
  count_a=0L
  while days lt daynum do begin
    
    ev = new_ev_struct()
    ev.npoints[*] = 524288
    
    timespan, t
    
    t1 = time_string(time_double(t)+86400)
    
    spp_fld_make_or_retrieve_cdf,'dfb_dbm_1',/load ;V12
    get_data, 'spp_fld_dfb_dbm_1_dbm_data', data = dbm_1_dat
    
    spp_fld_make_or_retrieve_cdf,'dfb_dbm_2',/load ;V34
    get_data, 'spp_fld_dfb_dbm_2_dbm_data', data = dbm_2_dat
    
    spp_fld_make_or_retrieve_cdf,'dfb_dbm_3',/load ;V5
    get_data, 'spp_fld_dfb_dbm_3_dbm_data', data = dbm_3_dat
    
    spp_fld_make_or_retrieve_cdf, 'mago_survey', /load
    get_data, 'spp_fld_mago_survey_mag_bx_nT', data=magx
    get_data, 'spp_fld_mago_survey_mag_by_nT', data=magy
    get_data, 'spp_fld_mago_survey_mag_bz_nT', data=magz
    
    datewhere = where(dbm_1_dat.x eq dbm_3_dat.x)
    datewhere2 = where(dbm_3_dat.x eq dbm_1_dat.x)
    ;help, datewhere
    
    n_bursts = n_elements(datewhere)
    tplot_clear, /all
    
    data = {date: dblarr(n_bursts,3), $
            datestr: strarr(n_bursts,3), $
            efield: dblarr(n_bursts,3,524288), $
            bx: magx.y, $
            by: magy.y, $
            bz: magz.y, $
            bfieldtime: magx.x $
            }
    
    
    
    data.date[*,0] = dbm_1_dat.x[datewhere]
    data.date[*,1] = dbm_2_dat.x[datewhere]
    
;    help,dbm_1_dat.y
;    help,dbm_3_dat.y
;    
;    print, dbm_1_dat.x[0]
;    print, dbm_3_dat.x[0]
;    
;    print, dbm_1_dat.x[79]
;    print, dbm_3_dat.x[79]
    
    data.date[*,2] = dbm_3_dat.x[datewhere2]
    
    data.datestr = time_string(data.date)
    
    data.efield[*,0,*] = dbm_1_dat.y[datewhere,*]
    data.efield[*,1,*] = dbm_2_dat.y[datewhere,*]
    data.efield[*,2,*] = dbm_3_dat.y[datewhere2,*]
    
    ;data.efield[*,0,*] = E12*cos((!pi/180.)*55.)-E34*cos((!pi/180.)*40.)
    ;data.efield[*,1,*] = E12*sin((!pi/180.)*55.)+E34*sin((!pi/180.)*40.)

    
    Bt = 3.49525 ;burst length (sec)
    Sr = 150000./5 ;samples/sec
    pi = !DPI
    q = 1.60217662*10.^(-19) ; coulombs
    mass_e = 9.10938356*10.^(-31) ;kg
    
    time_of_one_capture = indgen(52428*2, /DOUBLE)/Sr 
    
    mission = 'PSP'
    
    inst = 'DFB'
    
    case 1 of
      keyword_set(rotsc):   co = 'SC'
      keyword_set(rotminvar): co = 'MV'
      keyword_set(rotfa):   co = 'FA'
      else:         co = 'SW'
    endcase
    
    
    savefile = savepath + mission + '_' + 'ST' +$
      '_Narrow_'+co+'COORD'+'_'+$
      strmid(t,5,2)+'_'+strmid(t,8,2)+'.txt
      
    openw,savelun,savefile,/get_lun
    close,savelun

    saved1 = psp_narrowband_save(savelun,savefile,mission,inst,version,ev,count_a,t,t1,/bgn)
    
    
    ;print, time
  
    
    for k=0, n_bursts-1 do begin
      
      ev.source[0] = 'Ex'
      ev.source[1] = 'Ey'
      ev.source[2] = 'Ez'
      
      ;--------------------------------------
      ;ROTATE DATA TO SPACECRAFT COORD SYSTEM
      ;--------------------------------------
      
      data.efield[k,0,*] = data.efield[k,0,*] - mean(data.efield[k,0,*])
      data.efield[k,1,*] = data.efield[k,1,*] - mean(data.efield[k,1,*])
      data.efield[k,2,*] = data.efield[k,2,*] - mean(data.efield[k,2,*])
      
      ;effective antenna length: 7m for V12 and V34, 6m for V5
      
      data.efield[k,0,*] = data.efield[k,0,*]/3.5 ;V/m
      data.efield[k,1,*] = data.efield[k,1,*]/3.5 ;V/m
      data.efield[k,2,*] = data.efield[k,2,*]/4. ;V/m
      
      data.efield[k,0,*] = data.efield[k,0,*]*cos((!pi/180.)*55.)-data.efield[k,1,*]*cos((!pi/180.)*40.)
      data.efield[k,1,*] = data.efield[k,0,*]*sin((!pi/180.)*55.)+data.efield[k,1,*]*sin((!pi/180.)*40.)
     
      Ex = reform(data.efield[k,0,*])
      Ey = reform(data.efield[k,1,*])
      Ez = reform(data.efield[k,2,*])
      
      Ex = interpol(Ex,52428*2)
      Ey = interpol(Ey,52428*2)
      Ez = interpol(Ez,52428*2)
     
      ev.mv_data[*,0] = Ex
      ev.mv_data[*,1] = Ey
      ev.mv_data[*,2] = Ez
      
      ;---------------
      ;OTHER ROTATIONS
      ;---------------
      
      
;      if keyword_set(rotsc) then begin
;        for rotcount=0,ev.npoints[1]-1 do begin
;          tempvec = swrot_Vant2Esc([ev.mv_data[rotcount,0],ev.mv_data[rotcount,1],ev.mv_data[rotcount,2]])
;          ev.mv_data[rotcount,0] = tempvec[0]
;          ev.mv_data[rotcount,1] = tempvec[1]
;          ev.mv_data[rotcount,2] = tempvec[2]
;        endfor
;      endif
;
;      ;--------------------------------------------
;      ;ROTATE DATA TO MINIMUM VARIANCE COORD SYSTEM
;      ;--------------------------------------------
;      ;adapted from swtds
;
;      if keyword_set(rotminvar) then begin
;        ;first rotate SWAVES into SC coords with M2 matrix
;        for rotcount=0,ev.npoints[1]-1 do begin
;          tempvec = swrot_Vant2Esc([ev.mv_data[rotcount,0],ev.mv_data[rotcount,1],ev.mv_data[rotcount,2]])
;          ev.mv_data[rotcount,0] = tempvec[0]
;          ev.mv_data[rotcount,1] = tempvec[1]
;          ev.mv_data[rotcount,2] = tempvec[2]
;        endfor
;
;        rotmat=min_var(REFORM(ev.mv_data[*,0]), REFORM(ev.mv_data[*,1]), REFORM(ev.mv_data[*,2]), eig_vals=eigs)
;        
;        ;perform the rotation
;        rotated=reform([[[ev.mv_data[*,0]]],[[ev.mv_data[*,1]]],[[ev.mv_data[*,2]]]]) # rotmat
;        for rotcount=0,ev.npoints[1]-1 do begin
;          ev.mv_data[rotcount,0]=rotated[rotcount,0]
;          ev.mv_data[rotcount,1]=rotated[rotcount,1]
;          ev.mv_data[rotcount,2]=rotated[rotcount,2]
;        endfor
;      endif
;
;      ;-----------------------------------------
;      ;ROTATE DATA TO FIELD ALIGNED COORD SYSTEM
;      ;-----------------------------------------
;      ;adapted from swtds
;      if keyword_set(rotfa) then begin
;
;        ;first rotate SWAVES into SC coords with M2 matrix
;        FOR rotcount=0,npoints-1 DO BEGIN
;          tempvec=swrot_Vant2Esc([ev.mv_data[rotcount,0],ev.mv_data[rotcount,1],ev.mv_data[rotcount,2]])
;          ev.mv_data[rotcount,0]=tempvec[0]
;          ev.mv_data[rotcount,1]=tempvec[1]
;          ev.mv_data[rotcount,2]=tempvec[2]
;        ENDFOR
;
;        ;rotate ev_data into field aligned coordinates
;        ev = stereo_waves_rotfa(ev,sc,plasma_struct)
;
;      endif

      
      for chn=0,2 do begin
        
        timedif = abs(data.bfieldtime - data.date[k,chn])
        time = time_of_one_capture ;+ data.date[k,chn]
        timewhere = where(timedif eq min(timedif))
        
        b_field = sqrt(data.bx[timewhere]^2+data.by[timewhere]^2+data.bz[timewhere]^2) ; in nT
        b_field_tesla = b_field*10.^(-9) ; in Tesla
        
        channel = ev.source[chn]
        
        fce = (1./(2*pi))*(q*b_field_tesla)/mass_e
        ;p = plot(time,reform(data.efield[k,chn,*]), TITLE=channel+'  '+data.datestr[k,chn], XTITLE="Time (Seconds)")
        wa = psp_narrowband_analyze(data.datestr[k,chn], time, ev.mv_data[*,chn], fce[0], type=type, fft_step=0.05, step_overlap=0.5)
        
        ev.scetstr = data.datestr[k,chn]
        ev.quality[chn,*] = wa.quality[*]
        ev.freq[chn] = wa.freq_maxima
        ev.max_amp[chn] = wa.max_amp
        ev.f_fce[chn] = wa.f_fce
        ev.duration[chn] = wa.timelength
        ev.bw[chn] = wa.bw_maxima
        ev.df_f[chn] = wa.df_f_maxima
          
       ; lag2 = 2*indgen(ev.npoints[chn], /LONG) - (ev.npoints[chn] - 1)
        
       ; ac = a_correlate(ev.mv_data[*,chn]/max(ev.mv_data[*,chn],/nan),lag2)
        dust = 'no'
        if (chn eq 3) then begin
          max01 = max(abs(ev.mv_data[*,0]))/max(abs(ev.mv_data[*,1]))
          max12= max(abs(ev.mv_data[*,1]))/max(abs(ev.mv_data[*,2]))
          max23 = max(abs(ev.mv_data[*,2]))/max(abs(ev.mv_data[*,3]))

          if (((max01 gt 100) || (max01 lt 1./100.)) || $
            ((max12 gt 100) || (max12 lt 1./100.)) || $
            ((max23 gt 100) || (max23 lt 1./100.))) then dust='yes'

        endif


        ;------------------------------------------------------------------------------
        ;RESET WAVE QUALITY IF SIGNAL IS DUST OR FREQUENCY IS OUTSIDE OF ALLOWED RANGE
        ;------------------------------------------------------------------------------
      endfor ; this for loop is for channels
      
      ;---------------------------------------------
      ;wave angle and electron resonance calculation
      ;---------------------------------------------
      
      wf = [[Ex],[Ey],[Ez]]
      
      wf = {X: fltarr(n_elements(wf[*,1])), $
        Y: wf $
      }

      store_data, 'waveform', data=wf
      
      timedif = abs(data.bfieldtime - data.date[k,0])
      aaa = where(timedif eq min(timedif))
      
      Bx = data.bx[aaa]
      By = data.by[aaa]
      Bz = data.bz[aaa]
      
      bvectmp = [Bx,By,Bz]

      bmagtmp = sqrt(bvectmp[0]^2+bvectmp[1]^2+bvectmp[2]^2)

      Eminvar = rbsp_rotate_field_2_vec('waveform',bvectmp)

      get_data, 'waveform_FA_minvar', data=waveform
      
      Emax_Emin = mean(abs(waveform.y[*,0]))/mean(abs(waveform.y[*,1]))
      
      dispstruct = cold_plasma_dispersion(epol=Emax_Emin, dens=150, freq=mean(ev.freq), Bo=bmagtmp)
      ev.norm_angle = dispstruct.theta_kb
      ev.countev = dispstruct.cyclo_counterstream_res.Etots[0]*1000. ;convert to eV
      ev.coev = dispstruct.cyclo_costream_res.Etots[0]*1000. ;convert to eV
      
      ;----------
      ;DUST CHECK
      ;----------
      
      qwhere = where(ev.quality eq -2)
      ;stop
      if ((qwhere[0] eq -1) && (dust eq 'yes')) then ev.quality[*,*] = -1  ;identify dust that isn't nyquist

      if counter eq 10. then print,'PRESS "q" TO QUIT AND SAVE CURRENT RESULTS'

      saved2 = psp_narrowband_save(savelun,savefile,mission,inst,version,ev,count_a,t,t1,/mid)

      
      ev = new_ev_struct(keeptime=ev.scetur8)
      count_a+=1

      qtemp = 0.
      
      

      if counter eq 10. then counter = 0.
      counter = counter + 1

      read_quit=get_kbrd(0)
      if read_quit eq 'q' then goto, endpro

    endfor ; this for goes to bursts in one day
    
    
    savedfile = psp_narrowband_save(savelun,savefile,mission,inst,version,ev,count_a,t,t1,/ed)

    psp_narrowband_id,t,t1,savedfile

    free_lun, savelun, /FORCE
    
    t = time_string(time_double(t)+86400.)
    
    days+=1
    
  endwhile
  
  

  ;Catches errors and ensures that the username and password get reset even if a failure occurs in the program
  CATCH, Error_status
  
  if Error_status ne 0 then begin
    print, 'Error Index:', Error_status
    print, 'Error message: ', !ERROR_STATE.MSG
  endif
  
  endpro:

end