;Returns FIREBIRD II (FU3 and FU4) geometric factors, energy channel ranges, and cadence 
;for the collimated detectors and surface detector (FU3 only) based on the date.
;Note that the date has to fall within the start/stop times of a campaign. See list below or in Table IV
;in Johnson2020

;This comes from Arlo Johnson's 2020 paper (Tables 1,2,3)
;The FIREBIRD-II CubeSat mission: Focused investigations of relativistic electron burst intensity, range, and dynamics
;doi: 10.1063/1.5137905


;To calibrate from counts to flux (only for differential channels):
;   flux = counts/cadence/energy_width/geometric_factor


;Example usage:
;   x = firebird_get_calibration_counts2flux,'2017-12-05','3'


function firebird_get_calibration_counts2flux,date,fb


    if time_double(date) le time_double('2015-02-01') then begin
      print,'CHOSEN DATE OCCURS BEFORE START OF FIREBIRD MISSION ON 2015-02-01'
      return,-1
    endif
    


    ;----------------------------------------------------
    ;Campaign start and end dates for FU3 and FU4 
    ;----------------------------------------------------

    case fb of 
        
        '3': begin
            campaign_start = time_double(['2015-02-01','2015-03-21','2015-05-16',$
            '2015-07-03','2015-08-08','2015-11-15','2016-01-15',$
            '2016-05-20','2016-08-12','2016-12-21','2017-05-01',$
            '2017-07-01','2017-11-19','2018-02-27','2018-04-20',$
            '2018-06-25','2018-07-31','2018-09-17','2018-12-16',$
            '2019-01-24','2019-03-16','2019-05-05','2019-07-05',$
            '2019-09-10','2020-01-02'])
            campaign_end = time_double(['2015-02-21','2015-04-19','2015-06-15',$
            '2015-08-04','2015-09-04','2015-12-15','2016-02-03',$
            '2016-06-14','2016-09-07','2017-01-04','2017-05-21',$
            '2017-07-21','2017-12-14','2018-03-28','2018-05-13',$
            '2018-07-18','2018-08-20','2018-10-13','2019-01-10',$
            '2019-02-20','2019-04-10','2019-05-17','2019-07-29',$
            '2019-10-08','2020-01-27'])
        end

        '4': begin
            campaign_start = time_double(['2015-02-01','2015-03-21','2015-05-16',$
            '2015-07-03','2015-08-08','2015-11-15','2016-01-15',$
            '2016-06-09','2016-08-12','2016-12-21','2017-05-01',$
            '2017-07-01','2017-11-19','2018-02-27','2018-04-20',$
            '2018-06-25','2018-07-31','2018-09-17','2018-12-16',$
            '2019-01-24','2019-03-16','2019-05-05','2019-07-05',$
            '2019-09-10','2020-01-02'])
            campaign_end = time_double(['2015-02-21','2015-04-19','2015-06-15',$
            '2015-08-04','2015-09-04','2015-12-15','2016-02-03',$
            '2016-06-20','2016-09-07','2017-01-04','2017-05-21',$
            '2017-07-21','2017-12-14','2018-03-28','2018-05-13',$
            '2018-07-18','2018-08-20','2018-10-13','2019-01-10',$
            '2019-02-20','2019-04-10','2019-05-17','2019-07-29',$
            '2019-10-08','2020-01-27'])
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





    return, {campaign:campaign,$
             cadence:cadence[campaign_index],$
             g_factor_collimated:g_factor_collimated,$
             g_factor_surface:g_factor_surface,$
             energy_range_collimated:energy_range_collimated,$
             energy_range_surface:energy_range_surface} 



end