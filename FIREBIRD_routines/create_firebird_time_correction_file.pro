;Create file that has the time correction so that I can time-correct data without having to run 
;daily load routines (e.g. Mike's microburst list). 
;
;Creates files FU3_time_corrections.txt and FU4_time_corrections.txt
;
;This code only needs to be run once. 
;

fb = '4' 



t0 = time_double('2015-02-01')
t1 = time_double('2019-05-16')

ndays = (t1 - t0)/86400 


paths = get_project_paths()


tcorr = []
times = []

file_fail = 0.


timebase = dindgen(1440)*60.   ;time base to interpolate to (once/min)

for i=0.,ndays-1 do begin

  timespan,t0 + 86400.*i  

  tb = t0 + 86400.*i  + timebase

  firebird_load_context_data_cdf_file,fb,file_fail=file_fail
  
  tinterpol_mxn,'Count_Time_Correction',tb,/overwrite
  get_data,'Count_Time_Correction',data=dd

  if is_struct(dd) then begin 
    tcorr = [tcorr,dd.y]
    times = [times,dd.x]
  endif
  
endfor

;tmax = time_double('2019-05-16')
;goo = where(times lt tmax)
;times = times[goo]
;tcorr = tcorr[goo]

store_data,'tcorr',times,tcorr



openw,lun,paths.root+'FU'+fb+'_time_corrections.txt',/get_lun
printf,lun,'Time correction values (sec) for the FIREBIRD data'
for i=0.,n_elements(times)-1 do printf,lun,times[i],tcorr[i],format='(f11.0,5x,f9.3)'
close,lun
free_lun,lun

stop

end

