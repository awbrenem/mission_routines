;Rotate the dvac PSP data products Efield sensor coordinates (V12, V34, Vz) 
;to SC coordinates. 
;Returns tplot variable with '_SC' appended to name

;See 6.4 in Coordinate transforms between FIELDS instruments by T. Dudok de Wit and D. Malaspina


pro psp_rotate_dvac_efieldsensor2sc,var


  ;Rotating V12, V34 into SC coord (x,y)
  rot_mat = [[0.64524,-0.82228,0.],$
             [0.76897,0.57577,0.],$
             [0,0,1]]



  get_data,var,data=dd
  Esc = reform(rot_mat ## dd.y)
  store_data,var + '_SC',dd.x,Esc


end




