;------------------------------------------------
;ADD AUTO BIAS TO FLAG VALUES
;------------------------------------------------

;; AutoBias starts actively controlling the bias currents at V12 = -1.0 V,
;; ramping down the magnitude of the bias current so that when V12 = 0.0 V,
;; the bias current is very near to zero after starting out around -20
;; nA/sensor.

;; For V12 > 0.0 V, the bias current continues to increase (become more
;; positive), although at a slower rate, 0.2 nA/V or something like that.


;Auto Bias flag values. From 'rbsp?_efw_hsk_idpu_fast_TBD'
;Bit	Value	Meaning
;3	8	Toggles off and on every other cycle when AutoBias is;
;		active.
;2	4	One when AutoBias is controlling the bias, Zero when
;		AutoBias is not controlling the bias.
;1	2	One when BIAS3 and BIAS4 can be controlled by AUtoBias,
;		zero otherwise.
;0	1	One when BIAS1 and BIAS2 can be controlled by AUtoBias,
;		zero otherwise.





pro create_rbsp_autobias_list,sc,date

  sc = 'b'
  date = '2015-01-02'


  rbx = 'rbsp'+sc
  timespan,date


  rbsp_load_efw_hsk,probe=sc,/get_support_data



  ;Find times when auto biasing is active
  get_data,rbx+'_efw_hsk_idpu_fast_TBD',data=tbd
  tbd.y = floor(tbd.y)
  ab_flag = intarr(n_elements(tbd.x))


  ;Possible flag values for on and off
  ab_off = [1,2,3,8,10,11]
  ab_on = [4,5,6,7,12,13,14,15]

  goo = where((tbd.y eq 4) or (tbd.y eq 5) or (tbd.y eq 6) or (tbd.y eq 7) or (tbd.y eq 12) or (tbd.y eq 13) or (tbd.y eq 14) or (tbd.y eq 15))
  if goo[0] ne -1 then ab_flag[goo] = 1

  store_data,'ab_flag',data={x:tbd.x,y:ab_flag}


  if ~tdexists('ab_flag_interp',tr[0],tr[1]) then tinterpol_mxn,'ab_flag',times,/spline
  ;; tplot,['ab_flag','ab_flag_interp']

  get_data,'ab_flag_interp',data=ab_flag
  ab_flag = ab_flag.y
