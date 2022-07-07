;Returns FIREBIRD II (FU3 and FU4) geometric factors, energy channel ranges, and cadence 
;for the collimated detectors and surface detector (FU3 only) based on the date.
;Note that the date has to fall within the start/stop times of a campaign. See list below or in Table IV
;in Johnson2020. 


;This comes from Arlo Johnson's 2020 paper (Tables 1,2,3)
;The FIREBIRD-II CubeSat mission: Focused investigations of relativistic electron burst intensity, range, and dynamics
;doi: 10.1063/1.5137905


;********************************************************************
;Code uses the correct (GEANT-determined) geometric factors. See this following note from the hires files:

;#"CAVEAT": "Flux has been calculated using the analytic geometric factor of 9 cm^2 sr. 
;#Modeling in GEANT4 suggests the real geometric factor is closer to 6, so actual flux is about 50% larger than reported.


;Because of the above fact, it's best to use the "counts" (not flux) from the FIREBIRD files and then calibrate those into flux 
;using this code rather than directly using the flux values. 
;********************************************************************


;*********************************************************************
;NOTE ON HOW TO CALIBRATE FROM COUNTS TO FLUX (differential channels only)
;Cadence = msec 
;energy_width = keV 
;geometric_factor = cm2 - sr 

;***Counts as input (e.g. from FIREBIRD hires files gotten from firebird_load_data.pro)     
;   flux = counts/(cadence/1000.)/energy_width/geometric_factor

;***Counts/sec as input (e.g. from Mike's microburst auto-id code, which outputs "counts_s").
;   flux = counts_sec/energy_width/geometric_factor
;*********************************************************************



;Example usage (see firebird_load_context_data_cdf_file.pro):
;   x = firebird_get_calibration_counts2flux('2017-12-05','3')


function firebird_get_calibration_counts2flux,date,fb


    if time_double(date) lt time_double('2015-02-01') then begin
      print,'CHOSEN DATE OCCURS BEFORE START OF FIREBIRD MISSION ON 2015-02-01'
      return,-1
    endif



    ;----------------------------------------------------
    ;Campaign start and end dates for FU3 and FU4 
    ;----------------------------------------------------

    case fb of 
        
        '3': begin
            campaign_start = time_double(['2015-02-01','2015-03-20','2015-05-15',$
            '2015-07-02','2015-08-07','2015-11-14','2016-01-14',$
            '2016-05-19','2016-08-11','2016-12-20','2017-04-31',$
            '2017-06-30','2017-11-18','2018-02-26','2018-04-19',$
            '2018-06-24','2018-07-30','2018-09-16','2018-12-15',$
            '2019-01-23','2019-03-15','2019-05-04','2019-07-04',$
            '2019-09-09','2020-01-01'])
            campaign_end = time_double(['2015-02-22','2015-04-20','2015-06-16',$
            '2015-08-05','2015-09-05','2015-12-16','2016-02-04',$
            '2016-06-15','2016-09-08','2017-01-05','2017-05-25',$
            '2017-07-25','2017-12-15','2018-03-30','2018-05-14',$
            '2018-07-19','2018-08-24','2018-10-14','2019-01-11',$
            '2019-02-21','2019-04-11','2019-05-18','2019-07-30',$
            '2019-10-09','2020-01-28'])
        end

        '4': begin
            campaign_start = time_double(['2015-02-01','2015-03-20','2015-05-15',$
            '2015-07-03','2015-08-07','2015-11-14','2016-01-13',$
            '2016-06-08','2016-08-11','2016-12-20','2017-04-31',$
            '2017-06-30','2017-11-18','2018-02-26','2018-04-19',$
            '2018-06-24','2018-07-30','2018-09-16','2018-12-15',$
            '2019-01-23','2019-03-15','2019-05-04','2019-07-04',$
            '2019-09-09','2020-01-01'])
            campaign_end = time_double(['2015-02-22','2015-04-20','2015-06-17',$
            '2015-08-05','2015-09-05','2015-12-16','2016-02-04',$
            '2016-06-21','2016-09-08','2017-01-05','2017-05-26',$
            '2017-07-25','2017-12-15','2018-04-04','2018-05-14',$
            '2018-07-19','2018-08-25','2018-10-14','2019-01-11',$
            '2019-02-22','2019-04-11','2019-05-18','2019-07-30',$
            '2019-10-09','2020-01-28'])
        end

    else: begin
        print, 'ERROR: NEED TO SELECT FIREBIRD 3 OR 4'
        return,-1
    end 
    endcase



    ;----------------------------------------------------
    ;Sampling cadence for the different campaigns.
    ;Values are the same for each detector
    ;----------------------------------------------------

    cadence = [18.75,18.75,18.75,18.75,18.75,18.75,12.5,50,50,12.5,50,50,50,$
                50,50,50,50,50,50,50,50,12.5,50,50,50]



    ;----------------------------------------------------
    ;Determine which campaign corresponds to the requested date


    goo = where(time_double(date) ge campaign_start)
    campstart = campaign_start[goo[n_elements(goo)-1]]
    campstop =  campaign_end[goo[n_elements(goo)-1]]

    ;Check to see if the requested data falls b/t the start and stop times of 
    ;nearest campaign 
    tmp = ((time_double(date) ge campstart) and (time_double(date) le campstop))
    if tmp then begin
        campaign_index = goo[n_elements(goo)-1]
        campaign = campaign_index + 1   
    endif else begin 
        campaign_index = -1
        campaign = -1
        print,'ERROR: DATE OUTSIDE OF DESIGNATED CAMPAIGN TIMES'
        return,-1
    endelse



    cadence_fin = cadence[campaign_index]

    ;----------------------------------------------------
    ;Energy ranges for the different detectors
    ;I'll reference these with respect to FU3 since they're the same for FU4
    ;----------------------------------------------------

    if fb eq '3' then begin 
        if campaign le 20 then begin
            goo = fltarr(6,2)
            goo[0,0] = 231.0 & goo[0,1] = 299.7
            goo[1,0] = 299.7 & goo[1,1] = 407.6
            goo[2,0] = 407.6 & goo[2,1] = 554.8
            goo[3,0] = 554.8 & goo[3,1] = 770.7
            goo[4,0] = 770.7 & goo[4,1] = 1055.2
            goo[5,0] = 1055.2 & goo[5,1] = 1055.2
            energy_range_collimated = goo

            goo = fltarr(6,2)
            goo[0,0] = 176.4 & goo[0,1] = 240.9
            goo[1,0] = 240.9 & goo[1,1] = 342.2
            goo[2,0] = 342.2 & goo[2,1] = 480.3
            goo[3,0] = 480.3 & goo[3,1] = 683.0
            goo[4,0] = 683.0 & goo[4,1] = 950.1
            goo[5,0] = 950.1 & goo[5,1] = 950.1
            energy_range_surface = goo

            g_factor_collimated = [5.4,5.8,6.0,5.8,5.9,3.8]
            g_factor_surface = [13.1,13.0,13.8,13.9,14.1,11.2]
        endif 
        if campaign ge 21 then begin 
            goo = fltarr(6,2)
            goo[0,0] = 201.2 & goo[0,1] = 250.2
            goo[1,0] = 250.2 & goo[1,1] = 299.2
            goo[2,0] = 299.2 & goo[2,1] = 348.2
            goo[3,0] = 348.2 & goo[3,1] = 446.2
            goo[4,0] = 446.2 & goo[4,1] = 1055.2
            goo[5,0] = 1055.2 & goo[5,1] = 1055.2
            energy_range_collimated = goo

            energy_range_surface = goo 
            energy_range_surface[*] = !values.f_nan

            g_factor_collimated = [4.9,5.1,5.0,5.9,6.6,3.6]
            g_factor_surface = replicate(!values.f_nan,6)
        endif 
    endif  ;for FU-3



    if fb eq '4' then begin
        if campaign le 20 then begin

            goo = fltarr(6,2)
            goo[0,0] = 219.7 & goo[0,1] = 283.4
            goo[1,0] = 283.4 & goo[1,1] = 383.6
            goo[2,0] = 383.6 & goo[2,1] = 520.3
            goo[3,0] = 520.3 & goo[3,1] = 720.7
            goo[4,0] = 720.7 & goo[4,1] = 985.0
            goo[5,0] = 985.0 & goo[5,1] = 985.0
            energy_range_collimated = goo

            energy_range_surface = goo 
            energy_range_surface[*] = !values.f_nan

            g_factor_collimated = [5.5,5.4,5.7,5.9,5.8,4.2]
            g_factor_surface = replicate(!values.f_nan,6)

        endif 

        if campaign ge 21 then begin 

            goo = fltarr(6,2)
            goo[0,0] = 201.0 & goo[0,1] = 246.5
            goo[1,0] = 246.5 & goo[1,1] = 301.1
            goo[2,0] = 301.1 & goo[2,1] = 346.2
            goo[3,0] = 346.2 & goo[3,1] = 446.7
            goo[4,0] = 446.7 & goo[4,1] = 983.6
            goo[5,0] = 983.6 & goo[5,1] = 983.6
            energy_range_collimated = goo

            energy_range_surface = goo 
            energy_range_surface[*] = !values.f_nan

            g_factor_collimated = [4.9,5.3,4.9,5.9,6.5,4.2]
            g_factor_surface = replicate(!values.f_nan,6)
        endif 
    endif ;for FU-4



    ;Determine the energy channel used to create the context data 
    ;Johnson20 Table 3
    if fb eq '3' then begin 
        if ((campaign ge 1) and (campaign le 7)) then begin
            channel_type = 'surface'
            channel_for_context = 5
        endif 
        if ((campaign ge 8) and (campaign le 9)) then begin
            channel_type = 'collimated'
            channel_for_context = 1
        endif 
        if ((campaign ge 10) and (campaign le 20)) then begin
            channel_type = 'collimated'
            channel_for_context = 2
        endif 
        if ((campaign ge 21) and (campaign le 24)) then begin
            channel_type = 'collimated'
            channel_for_context = 3
        endif 
    endif else begin 
        if campaign eq 1 then begin
            channel_type = 'surface'
            channel_for_context = !values.f_nan  ;NOTE: SURFACE CHANNEL ON FU4 NEVER WORKED
        endif 
        if ((campaign ge 2) and (campaign le 7)) then begin
            channel_type = 'collimated'
            channel_for_context = 2
        endif 
        if ((campaign ge 8) and (campaign le 9)) then begin
            channel_type = 'collimated'
            channel_for_context = 1
        endif 
        if ((campaign ge 10) and (campaign le 20)) then begin
            channel_type = 'collimated'
            channel_for_context = 2
        endif 
;****NOTE: THE BELOW IS CORRECT. HOWEVER, THERE ARE MICROBURSTS FROM CAMPAIGNS > 24. 
;I'M NOT SURE WHAT THE VALUES ARE FOR THESE LATER CAMPAIGNS. 
;        if ((campaign ge 21) and (campaign le 24)) then begin
;            channel_type = 'collimated'
;            channel_for_context = 3
;        endif 
        if (campaign ge 21) then begin
          channel_type = 'collimated'
          channel_for_context = 3
        endif

    endelse 


    return, {campaign:campaign,$
             cadence:cadence[campaign_index],$
             channel_type_for_survey_data:channel_type,$
             channel_used_for_survey_calibration:channel_for_context,$
             g_factor_collimated:g_factor_collimated,$
             g_factor_surface:g_factor_surface,$
             energy_range_collimated:energy_range_collimated,$
             energy_range_surface:energy_range_surface} 



end