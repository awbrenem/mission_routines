;  mag_wave_angle calculates the angle of propogation 
;  relative to the background field vector given by b0
;
;  waveform is an [n,3] magnetic field waveform vector
;  
;  b0 is the background magnetic field. If the input is an [m,3] matrix, the average b0 will be 
;  taken as the background field vector, else if its a [3] vector, 
;  b0 is taken to be the field vector in a normal circumstance. 
;  
;
;
;

function mag_wave_angle, waveform, b0
  
  if n_elements(b0) gt 3 then begin
    
    btmp = fltarr(3)
    btmp[0] = mean(b0[*,0])
    btmp[1] = mean(b0[*,1])
    btmp[2] = mean(b0[*,2])
    
    b0 = btmp
    
  endif
  
  n_points = n_elements(waveform[*,0])
  
  og_avg_x = fltarr(40) 
  og_avg_y = fltarr(40)
  og_avg_z = fltarr(40)
  
  ;help, waveform
  
  for i=0,39 do begin
    
    og_avg_x[i] = mean(waveform[0+100*i:99+100*i,0])
    og_avg_y[i] = mean(waveform[0+100*i:99+100*i,1])
    og_avg_z[i] = mean(waveform[0+100*i:99+100*i,2])
    
  endfor
  
  ;ct = colortable(74, /REVERSE)
  
  ;p3d = PLOT3D(og_avg_x,og_avg_y,og_avg_z, 'o',/SYM_FILLED, RGB_TABLE=ct, SHADOW_COLOR='goldenrod',XTITLE='x',YTITLE='y')
  
  originalpoint = [og_avg_x[0],og_avg_y[0],og_avg_z[0]]
  
  og_arr = [[og_avg_x],[og_avg_y],[og_avg_z]]
  
  normal_arr = []

  for j=2,39 do begin
    
    vec1 = reform(og_arr[j,*])-originalpoint
    vec2 = reform(og_arr[j-1,*])-originalpoint
    
    normalvec = [vec2[1]*vec1[2]-vec2[2]*vec1[1],vec2[2]*vec1[0]-vec2[0]*vec1[2],vec2[0]*vec1[1]-vec2[1]*vec1[0]]
    
    norm_vec_norm = normalvec/sqrt(normalvec[0]^2+normalvec[1]^2+normalvec[2]^2)
    
    normal_arr = [[normal_arr],[norm_vec_norm]]
    
  endfor

  ;print, mean(normal_arr[0,*]), ' ', mean(normal_arr[1,*]), ' ', mean(normal_arr[2,*])
  ;print, b0/sqrt(b0[0]^2+b0[1]^2+b0[2]^2)
  
  normx = mean(normal_arr[0,*])
  normy = mean(normal_arr[1,*])
  normz = mean(normal_arr[2,*])
  
  
  norm_mag = sqrt(normx^2+normy^2+normz^2)
  b0_mag = sqrt(b0[0]^2+b0[1]^2+b0[2]^2)
  
  cosine = (normx*b0[0]+normy*b0[1]+normz*b0[2])/(b0_mag*norm_mag)
  
  degrees = 180/!PI*ACOS(cosine)
  
  degrees1 = 180 - degrees
  degrees2 = degrees
  
  if degrees1 lt degrees2 then degrees=degrees1
  
  
  return, degrees
end
 