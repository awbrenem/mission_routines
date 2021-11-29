;Create a list of bias sweep values using HSK data.
;NOTE: this is intended to be called in a loop of dates for
;each spacecraft. After the file is finished creating, this program
;no longer needs to be run.

;************************
;NOTE: THIS CODE IS OBSOLETE B/C I'M USING SHENG'S AND GABE'S LISTS FOR THE
;FINAL AUTOBIAS TIMES.
;************************


;----------------------------------------------
;I now have 3 bias sweep lists for each sc.
;1) Aaron's list from this code (rbsp?_bias_sweep_times_aaron.txt)
;2) Sheng's list (rbsp?_bias_sweep_times_sheng.txt)
;3) Gabe's list (rbsp?_bias_sweep_times_gabe.txt)
;----------------------------------------------


pro create_rbsp_bias_sweep_list,sc,date

;  date = '2015-11-14'
;  sc = 'a'


  rbx = 'rbsp'+sc

  timespan,date
  tr = timerange()


  rbsp_load_efw_hsk,probe=sc,/get_support_data




;--------------------------------------------------
;Determine times of bias sweeps
;--------------------------------------------------


  get_data, rbx+'_efw_hsk_beb_analog_CONFIG0', data = BEB_config
  if is_struct(BEB_config) then begin


      bias_sweep = intarr(n_elements(BEB_config.x))
      boo = where(BEB_config.y eq 64)
      if boo[0] ne -1 then bias_sweep[boo] = 1
      store_data,'bias_sweep',data={x:BEB_config.x,y:bias_sweep}

      diffv = bias_sweep - shift(bias_sweep,1)
      store_data,'bias_switch',BEB_config.x,diffv

      goo = where(diffv eq 1)
      if goo[0] ne -1 then startT = time_string(BEB_config.x[goo])
      goo = where(diffv eq -1)
      if goo[0] ne -1 then endT = time_string(BEB_config.x[goo])

      if goo[0] ne -1 then begin
        openw,lun,'~/Desktop/rbsp'+sc+'_bias_sweep_times.txt',/append,/get_lun
        for i=0,n_elements(startT)-1 do printf,lun,startT[i]+'  '+endT[i]

        close,lun & free_lun,lun

      endif

  endif
  store_data,tnames(),/delete
end
