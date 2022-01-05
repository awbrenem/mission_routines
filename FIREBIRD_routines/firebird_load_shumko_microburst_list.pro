;Load the list that Mike Shumko's FIREBIRD microburst id code ("microburst_detection" package)
;spits out. These are individual microbursts pulled from the FIREBIRD flux data. 
;These microbursts have been time-corrected, as comparison to my detrended data from firebird_subtract_tumble...pro shows 

;Listed microburst amplitudes ("sig0" through "sig5" in the .csv file) are actually signal/background ratio for each energy channel.
;To turn these into counts use:
;   counts = sig*sqrt(A+1) + A, 
;where A is the running average.
;Once in counts, these can be turned into fluxes using calibrations from 
;  firebird_get_calibration_counts2flux.pro. 
;Since the above conversions require loading of daily data, they aren't performed in this code. 


;fb = '3' or '4'  corresponding to FIREBIRD FU3 or FU4

function load_firebird_microburst_list,fb,filename=fntmp

    pathtmp = '/Users/abrenema/Desktop/code/Aaron/github/research_projects/RBSP_Firebird_Colpitts_Chen/microburst_detection/data/'

    if not keyword_set(fntmp) then fntmp = 'FU'+fb+'_microbursts.csv'


    ft = [7,4,4,4,4,4,4,4,4,4,4,4,4]
    fn = ['time','lat','lon','alt','mcilwainL','MLT','kp','flux_ch1','flux_ch2','flux_ch3','flux_ch4','flux_ch5','flux_ch6']
    floc = [ 0,27,45,64,82,101,119,124,143,161,179,198,217]
    fg = indgen(13)

    template = {version:1.,$
                datastart:1L,$
                delimiter:44B,$
                missingvalue:!values.f_nan,$
                commentsymbol:'',$
                fieldcount:13L,$
                fieldtypes:ft,$
                fieldnames:fn,$
                fieldlocations:floc,$
                fieldgroups:fg}



    vals = read_ascii(pathtmp+fntmp,template=template)

    vals.time = time_string(vals.time,prec=6)


    return,vals 

end