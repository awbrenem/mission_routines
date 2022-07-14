;Test time corrections on FIREBIRD context and hires data 
;See July 8th, 2021 email I sent to Arlo. 
;This has an example event (23:00:08.799) that Arlo confirms has the correct time-correction. 
;The event itself is a false microburst followed by a data dropout. 
;This day is a good one to compare time-corrections b/c the time-correction~30 sec.


;NOTES:
;(1) Note that for this date (2017-12-05 on FU3 - campaign 13 - see Arlo's paper) the second energy channel is used to calibrate the survey data.
;(see tplot variable created below called channel_used_for_survey_calibration3 from firebird_get_calibration_counts2flux.pro)
;(2) The time-corrections and proper flux calibrations are done within the load routines 
;firebird_load_data.pro and firebird_load_context_data_cdf_file.pro 



;RESULTS
;--The survey and hires data match up very closely. Both the flux values and the time-correcitons
;--The ephemeris values match up very closely. 




;--------------------------------------------------------------------------------------------------














rbsp_efw_init


;Time-correction about 30.5 sec for this day
timespan,'2017-12-05'

fb = '3'




firebird_load_data,'3'    ;time=corrected and calibrated 

firebird_load_context_data_cdf_file,'3'   ;time-corrected and calibrated

split_vec,'fu3_fb_col_hires_flux'

;Downsample hires data to 6-sec cadence of survey data
;Use second energy channel for comparison. 
rbsp_detrend,'fu3_fb_col_hires_flux_1',6.
options,'fu3_fb_col_hires_flux_1_smoothed','color',250


store_data,'fluxcomb2',data=['fu3_fb_col_hires_flux_1_smoothed','flux_context_FU3']
loadct,39

ylim,'fu3_fb_col_hires_flux_?',0,0,0

tplot,['fluxcomb2','fu3_fb_col_hires_flux_?']
stop

;See if the ephemeris context data need to be time-shifted
options,['fu3_fb_mcilwainL_from_hiresfile','McIlwainL'],'panel_size',1
ylim,['fu3_fb_mcilwainL_from_hiresfile','McIlwainL'],0,10
tplot,['fu3_fb_mcilwainL_from_hiresfile','McIlwainL']


options,['fu3_fb_mlt_from_hiresfile','MLT'],'panel_size',1
ylim,['fu3_fb_mlt_from_hiresfile','MLT'],0,24
tplot,['fu3_fb_mlt_from_hiresfile','MLT']

end


