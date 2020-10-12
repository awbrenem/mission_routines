
;Create tplot variables from LANL sat data sent to me by Geoff Reeves.
;Saves in filenames like lanl_sat_1994-084_20130103.tplot
;****TAKES PROBABLY 10 MIN FOR EACH FILE (single day, single LANL sat)

pro read_lanl_sat


;campaign = '2'
;date = ['20140101','20140102','20140103','20140104','20140105','20140106','20140107','20140108','20140109','20140110','20140111','20140112','20140113','20140114']
;folder = ['Jan1','Jan2','Jan3','Jan4','Jan5','Jan6','Jan7','Jan8','Jan9','Jan10','Jan11','Jan12','Jan13','Jan14','Jan15']
;satname = ['1991-080','1994-084','LANL-01A','LANL-02A','LANL-04A','LANL-97A']

campaign = '1'
;date = ['20130114','20130116','20130117','20130118','20130119','20130120','20130125','20130126','20130127','20130128','20130130','20130131','20130201','20130202','20130203','20130204','20130205','20130206']
;date = ['20130130','20130131','20130201','20130202','20130203','20130204','20130205','20130206']
;folder = ['Jan30','Jan31','Feb1','Feb2','Feb3','Feb4','Feb5','Feb6']

date = ['20130122']
folder = ['Jan22']
satname = ['1991-080','1994-084','LANL-01A','LANL-02A','LANL-04A','LANL-97A']




;For each date...
for b=0,n_elements(date)-1 do begin
;  path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events/'+folder[b]+'/'
  path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events_campaign'+campaign+'/'+folder[b]+'/'

  fn = [date[b]+'_1991-080_SOPA_ESP_v2.1.0.txt',$
  date[b]+'_1994-084_SOPA_ESP_v2.1.0.txt',$
  date[b]+'_LANL-01A_SOPA_ESP_v2.1.0.txt',$
  date[b]+'_LANL-02A_SOPA_ESP_v2.1.0.txt',$
  date[b]+'_LANL-04A_SOPA_ESP_v2.1.0.txt',$
  date[b]+'_LANL-97A_SOPA_ESP_v2.1.0.txt']

  ;For each LANL sat for a single date
  for l=0,n_elements(fn)-1 do begin



    openr,lun,path+fn[l],/get_lun
    jnk = ''
    for i=0,222 do readf,lun,jnk


    datetime = ''
    mlt = '' & mlat = '' & re = ''
    avg_flux_sopa_e = strarr(1,10)
    avg_flux_esp_e = strarr(1,9)
    avg_flux_sopa_esp_e = strarr(1,19)
    avg_flux_i = strarr(1,12)


    while not eof(lun) do begin

      readf,lun,jnk
      vals = strsplit(jnk,' ',/extract)

      datetime = [datetime,vals[0]]
      mlat = [mlat,vals[9]]
      mlt = [mlt,vals[11]]
      re = [re,vals[8]]

      avg_flux_sopa_e = [avg_flux_sopa_e,reform(vals[12:21],1,10)] ;(keV) [ 60.64, 86.78, 128.71, 185.23, 278.68, 395.38, 567.46, 821.18, 1175.14, 1494.84 ]
      avg_flux_esp_e =  [avg_flux_esp_e,reform(vals[22:30],1,9)] ;(keV) [ 1122.50, 1989.97, 2437.21, 3074.09, 3968.63, 5196.15, 6841.05, 9178.24, 16692.51 ]
      avg_flux_sopa_esp_e = [avg_flux_sopa_esp_e,reform(vals[31:49],1,19)] ;(keV) [ 60.64, 86.78, 128.71, 185.23, 278.68, 395.38, 567.46, 821.18, 1122.50, 1175.14, 1494.84, 1989.97, 2437.21, 3074.09, 3968.63, 5196.15, 6841.05, 9178.24, 16692.51 ]
      avg_flux_i = [avg_flux_i,reform(vals[50:61],1,12)] ;(keV) [ 61.24, 92.06, 138.6, 206.16, 316.23, 517.69, 896.66, 1509.97, 2426.93, 3937, 6204.84, 19621.4 ]

    endwhile

    close,lun & free_lun,lun

    nelem = n_elements(datetime)-1

    datetime = time_double(datetime[1:nelem])
    mlt = float(mlt[1:nelem]) & mlat = float(mlat[1:nelem]) & re = float(re[1:nelem])
    avg_flux_sopa_e = float(avg_flux_sopa_e[1:nelem,*])
    avg_flux_esp_e = float(avg_flux_esp_e[1:nelem,*])
    avg_flux_sopa_esp_e = float(avg_flux_sopa_esp_e[1:nelem,*])
    avg_flux_i = float(avg_flux_i[1:nelem,*])


    Lshell = re/(cos(!dtor*mlat)^2)  ;L-shell in centered dipole


    sopa_e_v = [ 60.64, 86.78, 128.71, 185.23, 278.68, 395.38, 567.46, 821.18, 1175.14, 1494.84 ]
    esp_e_v = [ 1122.50, 1989.97, 2437.21, 3074.09, 3968.63, 5196.15, 6841.05, 9178.24, 16692.51 ]
    sopa_esp_e_v = [ 60.64, 86.78, 128.71, 185.23, 278.68, 395.38, 567.46, 821.18, 1122.50, 1175.14, 1494.84, 1989.97, 2437.21, 3074.09, 3968.63, 5196.15, 6841.05, 9178.24, 16692.51 ]
    i_v = [ 61.24, 92.06, 138.6, 206.16, 316.23, 517.69, 896.66, 1509.97, 2426.93, 3937, 6204.84, 19621.4 ]

    store_data,satname[l]+'_'+'lshell',datetime,lshell
    store_data,satname[l]+'_'+'mlat',datetime,mlat
    store_data,satname[l]+'_'+'mlt',datetime,mlt
    store_data,satname[l]+'_'+'RE',datetime,re
    store_data,satname[l]+'_'+'avg_flux_sopa_e',datetime,avg_flux_sopa_e,sopa_e_v
    store_data,satname[l]+'_'+'avg_flux_esp_e',datetime,avg_flux_esp_e,esp_e_v
    store_data,satname[l]+'_'+'avg_flux_sopa_esp_e',datetime,avg_flux_sopa_esp_e,sopa_esp_v
    store_data,satname[l]+'_'+'avg_flux_i',datetime,avg_flux_i,i_v

    options,satname[l]+'_'+['avg_flux_sopa_e','avg_flux_esp_e','avg_flux_sopa_esp_e','avg_flux_i'],'spec',0
    ylim,satname[l]+'_'+['avg_flux_sopa_e','avg_flux_esp_e','avg_flux_sopa_esp_e','avg_flux_i'],1,1d5,1
    ylim,satname[l]+'_'+'mlat',-30,30
    ylim,satname[l]+'_'+'mlt',0,24
    ylim,satname[l]+'_'+'RE',5,7
;    tplot,satname[l]+'_'+['lshell','mlat','mlt','RE','avg_flux_sopa_e','avg_flux_esp_e','avg_flux_sopa_esp_e','avg_flux_i']


    tplot_save,'*',filename='~/Desktop/lanl_sat_'+satname[l]+'_'+date[b]

    store_data,tnames(),/delete

  endfor ;each LANL sat for a single date
endfor  ;each date

end
