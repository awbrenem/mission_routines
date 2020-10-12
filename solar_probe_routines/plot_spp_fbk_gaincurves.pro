;Plot the PSP filterbank gain curves sent to me by Dave Malaspina on 2019-05-03


rbsp_efw_init
path = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/solar_probe_routines/'


;AC values go from 104 - 89817 Hz
fn1 = 'PSP_FIELDS_DFB_AC_FilterBank_BandPass_Response_60dB_and_Above_20190502_DMM.sav'
restore,path+fn1
plot,ac_f_in_use,db_0,/nodata,xrange=[1,30000],/xlog
oplot,ac_f_in_use,db_0 & oplot,ac_f_in_use,db_1 & oplot,ac_f_in_use,db_2
oplot,ac_f_in_use,db_3 & oplot,ac_f_in_use,db_4 & oplot,ac_f_in_use,db_5
oplot,ac_f_in_use,db_6




;DC values go from 0.1 - 10000 Hz
fn2 = 'PSP_FIELDS_DFB_DC_FilterBank_BandPass_Response_60dB_and_Above_20190502_DMM.sav'
restore,path+fn2
plot,f_in_use,db_0,/nodata,xrange=[1,30000],/xlog
oplot,f_in_use,db_0 & oplot,f_in_use,db_1 & oplot,f_in_use,db_2
oplot,f_in_use,db_3 & oplot,f_in_use,db_4 & oplot,f_in_use,db_5
oplot,f_in_use,db_6 & oplot,f_in_use,db_7 & oplot,f_in_use,db_8
oplot,f_in_use,db_9 & oplot,f_in_use,db_10 & oplot,f_in_use,db_11
oplot,f_in_use,db_12 & oplot,f_in_use,db_13 & oplot,f_in_use,db_14
