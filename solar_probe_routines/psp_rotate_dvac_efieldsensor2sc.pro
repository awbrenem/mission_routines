;Rotate the dvac PSP data products Efield sensor coordinates (12, 34, z) 
;to SC coordinates. 
;Returns tplot variable with '_SC' appended to name

pro psp_rotate_dvac_efieldsensor2sc,var


  ;Rotation matrix from E-field sensor (uvw) to SC coord (Cindy emails in Sept 2020, and Marc Pulupa PPT slides)
  rot_mat = invert([[cos(55*!dtor), sin(55*!dtor), 0],$
              [-1*cos(40.*!dtor),sin(40.*!dtor), 0],$
              [0,0,1]])


  get_data,var,data=dd
  Esc = reform(rot_mat ## dd.y)
  store_data,var + '_SC',dd.x,Esc




end




