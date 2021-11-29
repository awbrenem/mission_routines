;Creates a .txt file that has all the EFW burst 2 start/stop times

;NOTE: this really only needs to be run once to create the .txt files.
;After this code has been run, to read in this list and create an IDL save file
;use create_burst2_times_idl_save_file.pro
;


pro create_rbsp_burst2_times_list,probe,date


  ;Output file with burst times and rates
  fn = '~/Desktop/burst2_times_RBSP'+probe+'.txt'
  ;Create file with header if it doesn't already exist
  ftst = file_test(fn)
  if not ftst then begin
    openw,lun,fn,/get_lun
    printf,lun,'Burst 2 times,                          duration (sec) for RBSP'+probe+', Rate=16,384 Samples/sec'
    close,lun & free_lun,lun
  endif


  timespan,date
  tr = timerange()



  ;Get burst start/stop times
  rbsp_load_efw_burst_times,probe=probe,/force_download,b2_times=b2t



  if KEYWORD_SET(b2t) then begin

    for i=0,n_elements(b2t[*,0])-1 do begin

      dur = b2t[i,1] - b2t[i,0]
      ;Sometimes burst duration is negative. Skip these...
      if dur gt 0. then begin
        burst_duration = strtrim(dur,1)
        burst_duration = string(burst_duration,format='(A4)')

        s1 = time_string(b2t[i,0])
        s2 = time_string(b2t[i,1])
        s3 = burst_duration

        openw,lun,fn,/get_lun,/append
        printf,lun,s1+' - '+s2+'  '+s3+'  '
        close,lun & free_lun,lun
      endif
    endfor
  endif


  store_data,tnames(),/del
end
