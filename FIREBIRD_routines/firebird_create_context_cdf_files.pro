;Take the FIREBIRD context files that span each campaign and divide them up into daily CDF files 
;for easier use.  
;Note that the survey files have only counts and not flux values.


;Written by Aaron W Breneman, Summer 2020
;-----------------------------------------------------------------------------------


path = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/research_projects/RBSP_Firebird_microburst_conjunctions_all/context_files/'
;fn = ['FU3_Context_camp01_2015-02-01_2015-02-22_L2.txt',$
;'FU3_Context_camp02_2015-03-20_2015-04-20_L2.txt',$
;'FU3_Context_camp03_2015-05-15_2015-06-15_L2.txt',$
;'FU3_Context_camp04_2015-07-03_2015-08-05_L2.txt',$
;'FU3_Context_camp05_2015-08-07_2015-09-05_L2.txt',$
;'FU3_Context_camp06_2015-11-15_2015-12-16_L2.txt',$
;'FU3_Context_camp07_2016-01-14_2016-02-04_L2.txt',$
;'FU3_Context_camp08_2016-05-20_2016-06-14_L2.txt',$
;'FU3_Context_camp09_2016-08-12_2016-09-07_L2.txt',$
;'FU3_Context_camp10_2016-12-21_2017-01-03_L2.txt',$
;'FU3_Context_camp11_2017-05-01_2017-05-25_L2.txt',$
;'FU3_Context_camp12_2017-07-01_2017-07-25_L2.txt',$
;'FU3_Context_camp13_2017-11-19_2017-12-14_L2.txt',$
;'FU3_Context_camp14_2018-02-26_2018-03-31_L2.txt',$
;'FU3_Context_camp15_2018-04-19_2018-05-13_L2.txt',$
;'FU3_Context_camp16_2018-06-24_2018-07-19_L2.txt',$
;'FU3_Context_camp17_2018-07-30_2018-08-24_L2.txt',$
;'FU3_Context_camp18_2018-09-17_2018-10-13_L2.txt',$
;'FU3_Context_camp19_2018-12-15_2019-01-10_L2.txt',$
;'FU3_Context_camp20_2019-01-23_2019-02-20_L2.txt',$
;'FU3_Context_camp21_2019-03-15_2019-04-10_L2.txt',$
;'FU3_Context_camp22_2019-05-04_2019-05-17_L2.txt']



fn = ['FU4_Context_camp01_2015-02-01_2015-02-22_L2.txt',$
'FU4_Context_camp02_2015-03-22_2015-04-20_L2.txt',$
'FU4_Context_camp03_2015-05-15_2015-06-15_L2.txt',$
'FU4_Context_camp04_2015-07-03_2015-08-02_L2.txt',$
'FU4_Context_camp05_2015-08-07_2015-09-05_L2.txt',$
'FU4_Context_camp06_2015-11-14_2015-12-14_L2.txt',$
'FU4_Context_camp07_2016-01-14_2016-02-04_L2.txt',$
'FU4_Context_camp08_2016-06-09_2016-06-20_L2.txt',$
'FU4_Context_camp09_2016-08-12_2016-09-07_L2.txt',$
'FU4_Context_camp10_2016-12-21_2017-01-02_L2.txt',$
'FU4_Context_camp11_2017-05-01_2017-05-26_L2.txt',$
'FU4_Context_camp12_2017-07-01_2017-07-25_L2.txt',$
'FU4_Context_camp13_2017-11-19_2017-12-12_L2.txt',$
'FU4_Context_camp14_2018-02-25_2018-04-04_L2.txt',$
'FU4_Context_camp15_2018-04-19_2018-05-14_L2.txt',$
'FU4_Context_camp16_2018-06-25_2018-07-19_L2.txt',$
'FU4_Context_camp17_2018-07-30_2018-08-25_L2.txt',$
'FU4_Context_camp18_2018-09-17_2018-10-13_L2.txt',$
'FU4_Context_camp19_2018-12-15_2019-01-11_L2.txt',$
'FU4_Context_camp20_2019-01-23_2019-02-23_L2.txt',$
'FU4_Context_camp21_2019-03-15_2019-04-11_L2.txt',$
'FU4_Context_camp22_2019-05-04_2019-05-17_L2.txt',$
'FU4_Context_camp23_2019-07-03_2019-07-30_L2.txt',$
'FU4_Context_camp24_2019-09-10_2019-10-08_L2.txt',$
'FU4_Context_camp25_2020-01-02_2020-01-27_L2.txt',$
'FU4_Context_camp26_2020-02-12_2020-03-09_L2.txt',$
'FU4_Context_camp27_2020-05-03_2020-05-28_L2.txt',$
'FU4_Context_camp28_2020-07-07_2020-07-08_L2.txt']


fb = strmid(fn[0],0,3)

;Final list of variables to NOT delete
varsave_general = ['epoch',$
                  'D0',$
                  'D1',$
                  'Count_Time_Correction',$
                  'McIlwainL',$
                  'Lat',$
                  'Lon',$
                  'Alt',$
                  'MLT',$
                  'kp',$
                  'Flag',$
                  'Loss_cone_type']




;----------------------------------------------------------------
;Grab the skeleton file.
skfile = fb + '_context_00000000_v01.cdf'
pathsk = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/research_projects/RBSP_Firebird_microburst_conjunctions_all/'
skeleton = pathsk + skfile
;make sure we have the skeleton CDF
found = 1
skeletonFile=file_search(skeleton,count=found)
if ~found then begin
  dprint,'Could not find skeleton CDF, returning.'
  stop
endif
skeletonFile = skeletonFile[0]
;----------------------------------------------------------------




;Define the structure of the FIREBIRD text files that we load in
ftype = [7,4,4,4,4,4,4,4,4,4,4,3]
fns = ['time','d0','d1','count_time_correction','L','lat','lon','alt','mlt','Kp','flag','lct']
floc = [0,27,33,38,56,74,92,112,130,149,154,158]
fg = indgen(12)

template = {version:1.,$
            datastart:98L,$
            delimiter:32B,$
            missingvalue:!values.f_nan,$
            commentsymbol:'',$
            fieldcount:12,$
            fieldtypes:ftype,$
            fieldnames:fns,$
            fieldlocations:floc,$
            fieldgroups:fg}






;for each campaign
for b=0,n_elements(fn)-1 do begin




  ;Load in all the data for current campaign
  x = read_ascii(path+fn[b],template=template)




  ;Divide up times by days
  nelem = n_elements(x.time)

  d0 = strmid(x.time[0],0,10)
  d1 = strmid(x.time[nelem-1],0,10)

  tunix = time_double(x.time)

  ndays = 1+(time_double(d1) - time_double(d0))/86400


  ;Loop through each day and save data in CDF file
  dcurr = time_double(d0)
  for i=0,ndays-1 do begin

    difft = tunix - dcurr
    goo = where((difft gt 0.) and (difft le 86400.))

    ;If data exists for current day
    if goo[0] ne -1 then begin

      ;Epoch times for the CDF file
      epoch = tplot_time_to_epoch(x.time[goo],/epoch16)




      ;--------------------------------------------------
      ;Rename skeleton file
      dtmp = time_string(dcurr)
      datestr = strmid(dtmp,0,4)+strmid(dtmp,5,2)+strmid(dtmp,8,2)
      fnfin = fb + '_context_'+datestr+'_v01.cdf'
      datafile = path + fnfin
      file_copy, skeletonFile, datafile, /overwrite ; Force to replace old file.
      cdfid = cdf_open(datafile)


      ;----------------------------------------------------


      ;Now that we have renamed some of the variables to our liking,
      ;get list of all the variable names in the CDF file.
      inq = cdf_inquire(cdfid)
      CDFvarnames = ''
      for varNum = 0, inq.nzvars-1 do begin $
        stmp = cdf_varinq(cdfid,varnum,/zvariable) & $
        if stmp.recvar eq 'VARY' then CDFvarnames = [CDFvarnames,stmp.name]
      endfor
      CDFvarnames = CDFvarnames[1:n_elements(CDFvarnames)-1]



      ;Delete all variables we don't want to save.
      for qq=0,n_elements(CDFvarnames)-1 do begin $
        tstt = array_contains(varsave_general,CDFvarnames[qq]) & $
        if not tstt then print,'Deleting var:  ', CDFvarnames[qq]
        if not tstt then cdf_vardelete,cdfid,CDFvarnames[qq]
      endfor


      cdf_varput,cdfid,'epoch',epoch
      cdf_varput,cdfid,'D0',transpose(x.d0[goo])
      cdf_varput,cdfid,'D1',transpose(x.d1[goo])
      cdf_varput,cdfid,'Count_Time_Correction',transpose(x.count_time_correction[goo])
      cdf_varput,cdfid,'McIlwainL',transpose(x.l[goo])
      cdf_varput,cdfid,'Lat',transpose(x.lat[goo])
      cdf_varput,cdfid,'Lon',transpose(x.lon[goo])
      cdf_varput,cdfid,'Alt',transpose(x.alt[goo])
      cdf_varput,cdfid,'MLT',transpose(x.mlt[goo])
      cdf_varput,cdfid,'kp',transpose(x.kp[goo])
      cdf_varput,cdfid,'Flag',transpose(x.flag[goo])
      cdf_varput,cdfid,'Loss_cone_type',transpose(x.lct[goo])


      cdf_close, cdfid

    endif

    dcurr += 86400d


  endfor ;each day
endfor  ;each campaign



end
