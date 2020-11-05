;+
  ; :Author: benshort
  ;-
  ;+
  ; :Description:
  ;    Describe the procedure.
  ;
  ; :Params:
  ;    t0
  ;
  ;
  ;
  ; :Author: benshort
  ;-
pro look_at_polar, t0
  
  if getenv('PSP_STAGING_ID') EQ '' then setenv, 'PSP_STAGING_ID=short186'
  if getenv('PSP_STAGING_PW') EQ '' then setenv, 'PSP_STAGING_PW=poLuSHISheoc'
  
  timespan, t0
  t = t0
  
  Bt = 3.49525 ;burst length (sec)
  sr = 150000. ;samples/sec
  pi = 3.14159
  q = 1.60217662*10.^(-19) ; coulombs
  mass_e = 9.10938356*10.^(-31) ;kg
  
  chn = 1
  
  change_channel:
  
  spp_fld_make_or_retrieve_cdf,'mago_survey', /load
  get_data, 'spp_fld_mago_survey_mag_bx_nT', data=magx
  get_data, 'spp_fld_mago_survey_mag_by_nT', data=magy
  get_data, 'spp_fld_mago_survey_mag_bz_nT', data=magz
  
  if chn eq 1 then begin
    spp_fld_make_or_retrieve_cdf,'dfb_dbm_1',/load
    get_data, 'spp_fld_dfb_dbm_1_dbm_data', data = dbm_dat_1

    spp_fld_make_or_retrieve_cdf,'dfb_dbm_2',/load
    get_data, 'spp_fld_dfb_dbm_2_dbm_data', data = dbm_dat_2

    spp_fld_make_or_retrieve_cdf,'dfb_dbm_3',/load
    get_data, 'spp_fld_dfb_dbm_3_dbm_data', data = dbm_dat_3
    
  endif else if chn eq 2 then begin
    
    spp_fld_make_or_retrieve_cdf,'dfb_dbm_4',/load
    get_data, 'spp_fld_dfb_dbm_4_dbm_data', data = dbm_dat_1

    spp_fld_make_or_retrieve_cdf,'dfb_dbm_5',/load
    get_data, 'spp_fld_dfb_dbm_5_dbm_data', data = dbm_dat_2

    spp_fld_make_or_retrieve_cdf,'dfb_dbm_6',/load
    get_data, 'spp_fld_dfb_dbm_6_dbm_data', data = dbm_dat_3
    
  endif
  
;  spp_fld_make_or_retrieve_cdf,'dfb_dbm_4',/load
;  get_data, 'spp_fld_dfb_dbm_4_dbm_data', data = dbm_dat_1
;  
;  spp_fld_make_or_retrieve_cdf,'dfb_dbm_5',/load
;  get_data, 'spp_fld_dfb_dbm_5_dbm_data', data = dbm_dat_2
;  
;  spp_fld_make_or_retrieve_cdf,'dfb_dbm_6',/load
;  get_data, 'spp_fld_dfb_dbm_6_dbm_data', data = dbm_dat_3
  
  tplot_clear, /all
  
  n_bursts = n_elements(dbm_dat_1.x)
  wavetime = dblarr(n_bursts, 524288)
  
  for n=0,n_bursts-1 do begin
    t1 = dbm_dat_1.x[n]
    times = dindgen(524288)/sr + t1
    wavetime[n,*] = times[*]
  endfor
  
  dbm_data = { starttimes: dbm_dat_1.x,$
    times: wavetime,$
    dataV12: dbm_dat_1.y,$
    dataV34: dbm_dat_2.y,$
    dataV5: dbm_dat_3.y,$
    magtime: magx.x,$
    magx: magx.y, $
    magy: magy.y, $
    magz: magz.y }
  
  ;help, dbm_data.dataV12
  i=0
  j=0
  k=1
  l=1
  repeat begin
    
    w1 = window(DIMENSIONS=[1200,800])
    
    E12 = dbm_data.dataV12[i,*] - mean(dbm_data.dataV12[i,*]) ; convert to only the perturbation of the electric field
    E34 = dbm_data.dataV34[i,*] - mean(dbm_data.dataV34[i,*]) ;
    E5 = dbm_data.dataV5[i,*] - mean(dbm_data.dataV5[i,*])    ;
    
    if chn eq 1 then begin
      E12 = reform(E12)/5.
      E34 = reform(E34)/5.
      E5 = reform(E5)/6.

      Ex = E12*cos((!pi/180.)*55.)-E34*cos((!pi/180.)*40.)
      Ey = E12*sin((!pi/180.)*55.)+E34*sin((!pi/180.)*40.)
      Ez = E5

      Ex = interpol(Ex,52428*2)
      Ey = interpol(Ey,52428*2)
      Ez = interpol(Ez,52428*2)
    endif else begin
      
      Ex = interpol(E12,52428*2)
      Ey = interpol(E34,52428*2)
      Ez = interpol(E5,52428*2)
      
    endelse
    
    
    
    print, chn

;    enew = fltarr(104857)
;
;    ; 524288/104857 = 
;
;    for l=0,104856 do begin
;
;      enew[l] = mean(Ey[floor(5*l):floor(5*(l+1))-1])
;
;    endfor
 
    ;Ex = interpol(Ex,52428)
    
    j1 = 40000
    
    
    if k mod 2 eq 0 then begin
      Ex = Ex[0+j1:10000+j1]
      Ey = Ey[0+j1:10000+j1]
      Ez = Ez[0+j1:10000+j1]
    endif
    
    wf = [[Ex],[Ey],[Ez]]
   
    wf = {X: fltarr(n_elements(wf[*,1])), $
      Y: wf $
    }
    
    store_data, 'waveform', data=wf
    
    magdif = abs(dbm_data.magtime - dbm_data.starttimes[i])
    aaa = where(magdif eq min(magdif))
    
    Bx = dbm_data.magx[aaa]
    By = dbm_data.magy[aaa]
    Bz = dbm_data.magz[aaa]
    
    bvectmp = [Bx,By,Bz]
    ;help, bvectmp
    bmagtmp = sqrt(bvectmp[0]^2+bvectmp[1]^2+bvectmp[2]^2)
    ;help, bmagtmp
    Eminvar = rbsp_rotate_field_2_vec('waveform',bvectmp)
    ;help, eminvar
    get_data, 'waveform_FA_minvar', data=waveform
    ;help,waveform
    if l mod 2 eq 0 then begin
      Ex = waveform.y[*,0]
      Ey = waveform.y[*,1]
      Ez = waveform.y[*,2]
    endif
    

    
    Emax_Emin = mean(abs(waveform.y[*,0]))/mean(abs(waveform.y[*,1]))
    ;print, emax_emin
    dispstruct = cold_plasma_dispersion(epol=Emax_Emin, dens=100, freq=1400, Bo=bmagtmp)
    ;help,dispstruct
    ;print, dispstruct.theta_kb, dispstruct.cyclo_counterstream_res.Etots[0]*1000.
    
    
    p11 = plot(Ex, color='blue',YTITLE='Ex', LAYOUT=[1,3,1],$
      TITLE=time_string(dbm_data.starttimes[i])+' Waveforms', /CURRENT)
      
    p12 = plot(Ey, color='blue',YTITLE='Ey', LAYOUT=[1,3,2],$
      /CURRENT)
       
    p13 = plot(Ez, color='blue',YTITLE='Ez', LAYOUT=[1,3,3],$
      /CURRENT)
    
    
    
    w2 = window(DIMENSIONS=[1200,500])
    p21 = plot(Ex,Ey,color='blue',XTITLE='Ex',YTITLE='Ey',$
       LAYOUT=[3,1,1],/CURRENT)
       
    p22 = plot(Ex,Ez,color='blue', XTITLE='Ex',YTITLE='Ez',$
       TITLE=time_string(dbm_data.starttimes[i])+' Hodograms', LAYOUT=[3,1,2],/CURRENT)
       
    p23 = plot(Ey,Ez,color='blue', XTITLE='Ey',YTITLE='Ez',$
       LAYOUT=[3,1,3],/CURRENT)
    
    if j eq 0 then begin
      print, "'p' for next plot. 'o' for previous. 'Enter' to exit."
      print, "'k' for zoom in. 'l' for minvar coordinates"
      print, "'SHIFT'+'1','2' to change field type."
    endif
    
    ;p = plot(enew)
    
    b = get_kbrd()
    
    if b eq 'p' then i+=1
    if b eq 'o' then i-=1
    
    if b eq 'k' then k+=1
    if b eq 'l' then l+=1
    
    if b eq '!' then begin
      chn = 1
      goto, change_channel
    endif
    
    if b eq '@' then begin
      chn = 2
      goto, change_channel
    endif
    w1.close
    w2.close
    j+=1
  endrep until byte(b) eq 10 ;press enter
end