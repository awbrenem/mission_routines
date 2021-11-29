;Creates the master list of bias sweep times by combining
;Sheng's and Gabe's lists.
;Gabe's used the HSK variable SDRPWRCTL (SDRAM protection bit) while Sheng's used the Usher and
;Guard voltages.
  ;NOTE: from Bonnell email on June 23, 2020:
  ;Here’s one reason why the 2019 “events” may look different in the HSK data from the earlier events:
  ;Most of the 2019 sweeps were designed to just look at the changes in saturation
  ;photocurrent, rather than to determine optimal current and voltage bias settings.
  ;Because of that, all or most of the sweeping was confined to IBIAS;
  ;the Usher and Guard voltages in a 2019 sweeps may not have changed at all,
  ;but the SDRAM protection bit was still being toggled (had to happen so that
  ;the SDT program could run at all).
  ;I’ve not looked closely at al, the HSK from the 2019 sweeps, but I remember
  ;that while I wanted to not sweep the guards and ushers at all,
  ;I may have had to because of the details of how the loops in SDT were
  ;implemented in assembly language.


;Aaron's list:
;   My list divides up the sweeps into too many chunks.
;   Doesn't work well near end of mission.
;Gabe's list:
;   has unusually long sweeps that are artificial (many days)
;   has sweeps that last < 1 min that aren't in Sheng's or Aaron's list
;Sheng's list:
;   Sweeps tend to last a bit longer (few min)
;   Doesn't pick up 2019 sweeps (except for one)  SEE above BONNELL note
;   Shouldn't be any false positives


;I simply combine the two lists to get the final list. If one list has the sweep
;lasting longer then I use the longer time.

;Final lists are called rbsp?_bias_sweep_times.txt



  rbsp_efw_init

  path = '/Users/aaronbreneman/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/calibration_files/'

  sc = 'b'


  fns = ['rbsp'+sc+'_bias_sweep_times_sheng.txt',$
         'rbsp'+sc+'_bias_sweep_times_gabe.txt',$
         'rbsp'+sc+'_bias_sweep_times_aaron.txt']


  jnk = ''
  t0s = '' & t1s = t0s
  t0g = '' & t1g = t0g
  t0a = '' & t1a = t0a

  for i=0,2 do begin
    openr,lun,path+fns[i],/get_lun
    while not eof(lun) do begin
      readf,lun,jnk
      goo = strsplit(jnk,' ',/extract)

      if i eq 0 then begin
        t0s = [t0s,goo[0]]
        t1s = [t1s,goo[1]]
      endif
      if i eq 1 then begin
        t0g = [t0g,goo[0]]
        t1g = [t1g,goo[1]]
      endif
      if i eq 2 then begin
        t0a = [t0a,goo[0]]
        t1a = [t1a,goo[1]]
      endif

    endwhile


    close,lun & free_lun,lun

  endfor

  t0ad = time_double(t0a) & t1ad = time_double(t1a)
  t0sd = time_double(t0s) & t1sd = time_double(t1s)
  t0gd = time_double(t0g) & t1gd = time_double(t1g)

  ;Create tplot variables of the sweep times
  ndays = 365 * 8.
  ntimes = ndays * 1440.
  tst = time_double('2012-01-01')
  tbase = 60.*dindgen(ntimes) + tst

  va = fltarr(n_elements(tbase))
  for i=0,n_elements(t0a)-1 do begin $
    goo = where((tbase ge t0ad[i]) and (tbase le t1ad[i])) & $
    if goo[0] ne -1 then begin
      dur = tbase[goo[n_elements(goo)-1]] - tbase[goo[0]]
      if (dur gt 60.) and (dur lt 5*60.*60.) then va[goo] = 1
    endif
  endfor
  store_data,'sweep_aaron_rbsp'+sc,tbase,va


  vs = fltarr(n_elements(tbase))
  for i=0,n_elements(t0s)-1 do begin $
    goo = where((tbase ge t0sd[i]) and (tbase le t1sd[i])) & $
    if goo[0] ne -1 then begin
      dur = tbase[goo[n_elements(goo)-1]] - tbase[goo[0]]
      if (dur gt 60.) and (dur lt 5*60.*60.) then vs[goo] = 1
    endif
  endfor
  store_data,'sweep_sheng_rbsp'+sc,tbase,vs


  vg = fltarr(n_elements(tbase))
  for i=0,n_elements(t0g)-1 do begin $
    goo = where((tbase ge t0gd[i]) and (tbase le t1gd[i])) & $
    if goo[0] ne -1 then begin
      dur = tbase[goo[n_elements(goo)-1]] - tbase[goo[0]]
      if (dur gt 60.) and (dur lt 5*60.*60.) then vg[goo] = 1
    endif
  endfor
  store_data,'sweep_gabe_rbsp'+sc,tbase,vg


  ;Take difference b/t Gabe's and Sheng's lists.
  ;+/- values mean that Sheng's/Gabe's list is missing a "sweep"
  diffv = vg - vs
  store_data,'diffv',tbase,diffv
  ylim,'diffv*',-2,2

  rbsp_detrend,'diffv',10*60. & ylim,'diffv_smoothed',-2,2

  ylim,'sweep*',0,2
  tplot,['sweep_sheng_rbsp'+sc,'sweep_gabe_rbsp'+sc,'diffv','diffv_smoothed'] & stop

  get_data,'diffv_smoothed',data=dd
  goon = where((dd.y gt -0.7) and (dd.y lt 0.))
  goop = where((dd.y lt 0.7) and (dd.y gt 0.))
  dd.y[goon] = 0.
  dd.y[goop] = 0.

  store_data,'diffv_smoothed_adj',data=dd & ylim,'diffv_smoothed_adj',-2,2
  tplot,['sweep_sheng_rbsp'+sc,'sweep_gabe_rbsp'+sc,'diffv','diffv_smoothed_adj'] & stop


  ;construct final SDT test list that's a combination of Sheng's and Gabe's lists.
  ;1) use Sheng's list since his events are almost always longer than Gabe's.
  ; There are also no false positives in this list that I can tell.
  ;2) Add in the times from Gabe's list that aren't in Sheng's list

  finallist = vs
  get_data,'diffv_smoothed_adj',data=dd
  goop2 = where(dd.y gt 0.)

  finallist[goop2] = 1.
  store_data,'finalsweeps',dd.x,finallist & ylim,'finalsweeps',0,2
  tplot,['sweep_sheng_rbsp'+sc,'sweep_gabe_rbsp'+sc,'diffv_smoothed_adj','finalsweeps'] & stop

  ;smooth the final list somewhat to remove blips
  rbsp_detrend,'finalsweeps',60.*10. & ylim,'finalsweeps_smoothed',0,2
  tplot,['diffv_smoothed_adj','finalsweeps','finalsweeps_smoothed'] & stop

  get_data,'finalsweeps_smoothed',data=dd
  goo = where(dd.y lt 0.2)
  dd.y[goo] = 0.
  boo = where(dd.y gt 0.2)
  dd.y[boo] = 1.
  doo = where(dd.y lt 1.)
  dd.y[doo] = 0.
  store_data,'finalsweeps_smoothed_adj',data=dd & ylim,'finalsweeps_smoothed_adj',0,2
  tplot,['sweep_sheng_rbsp'+sc,'sweep_gabe_rbsp'+sc,'diffv_smoothed_adj','finalsweeps','finalsweeps_smoothed','finalsweeps_smoothed_adj']



;stop
;  tplot,['fina','finb']
;
;  tplot,['fina','sweep_aaron_rbspa']
;  tplot,['finb','sweep_aaron_rbspb']


  ;Now turn the final times into a list of start/stop times.
  get_data,'finalsweeps',data=d
  vshift = d.y - shift(d.y,1)
  store_data,'vshifta',d.x,vshift & ylim,'vshifta',-2,2
  tplot,['finalsweeps','vshifta']

  startt = where(vshift eq 1)
  stopt = where(vshift eq -1)


  svals = fltarr(n_elements(tbase))
  svals[startt] = 1.
  store_data,'startt',tbase,svals

  stvals = fltarr(n_elements(tbase))
  stvals[stopt] = 1.
  store_data,'stopt',tbase,stvals


  store_data,'comb',data=['finalsweeps','startt','stopt'] & ylim,'comb',0,2
  options,'startt','colors',50
  options,'stopt','colors',254
  tplot,'comb'



  for i=0,1000 do print,time_string(d.x[startt[i]])+' '+time_string(d.x[stopt[i]])


stop

end
