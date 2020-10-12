;+
;FUNCTION: swrot_Vant2Esc.pro
;
;PURPOSE: Transforms measured SWAVES antenna voltages to electric field in
;         STEREO spacecraft coordinates using effective boom lengths and angles
;         as determined by the Graz wire-grid simulations and rheometry
;         measurements reported in "The electric antennas for the STEREO/WAVES
;         experiment", Bale, et. al.
;
;ARGUMENTS:
;       VEC   <-  SWAVES antenna voltage vector, Vx, Vy, Vz
;
;KEYWORDS: N/A
;
;RETURNS:  electric field vector in STEREO spacecraft coordinates
;
;CALLING SEQUENCE:
;       IDL> swrot_Vant2Esc,vec
;
;NOTES: (none)
;-
;CREATED BY:    Kris Kersten, July 2007
;
;MODIFICATION HISTORY:
;       7/2007  created, KK
;       4/2011  changed to double precision with built-in !DTOR degree to
;               radian conversion, KK
;
;INCLUDED MODULES:
;		swrot_Vant2Esc
;
;LIBS USED: (none)
;
;DEPENDENCIES: (none)
;-

function swrot_Vant2Esc,vec

; the following parameters are from the Graz wire-grid simulation
; and rheometry measurements with base caps
; NOTE, !DTOR=DOUBLE(pi/180), builtin double precision degree -> radian conv.
  thetaX=120.2D0*!DTOR
  thetaY=114.5D0*!DTOR
  thetaZ=124.5D0*!DTOR
  phiX=-135.0D0*!DTOR
  phiY=127.1D0*!DTOR
  phiZ=15.5D0*!DTOR
  hX=1.17D0
  hY=1.44D0
  hZ=0.97D0

; rotation using effective boom lengths and angles
  m2_inv=[[hX*cos(thetaX),hX*sin(thetaX)*sin(phiX),-hX*sin(thetaX)*cos(phiX)],$
          [hY*cos(thetaY),hY*sin(thetaY)*sin(phiY),-hY*sin(thetaY)*cos(phiY)],$
          [hZ*cos(thetaZ),hZ*sin(thetaZ)*sin(phiZ),-hZ*sin(thetaZ)*cos(phiZ)]]      
  m2=invert(m2_inv,/double)
  
  return,transpose(m2##vec)
end