pro EXformSC,sc,freq,dens,invMatrix
; print,'new'
;
;this routine calculates the "inverse" matrix, which transforms antenna 
;  voltages to electric fields in the spacecraft system, as in
;		ESCSY = invmatrix##VZYX, where VZYX are
;  the antenna signals, in order Z,Y,X and ESCSY is the electric field in
;  the spacecraft system, order X,Y,Z
;	NOW ORDER CHANGED TO VXYZ in Antenna_Vector and _Length
;
  Version = '30 Sep 2008'
  Version = '22 May 2017'			; low freq antenna length approx  
  Version = '26 June 2017'          ;  new antenna_length used
  Version = '30 June 2017'			; above saved as _OLD and comvention
;									; changed to order XYZ, i.e. UMN order
;									; both here and in Antenna_Vector and Length  
;      
; checked with antenna_check.pro  and antenna_test.pro
;
antL = fltarr(3)
Fmatrix =  fltarr(3,3) 
;print,’calling antenna_length from EXformSC ‘,sc,freq,dens
antenna_length, sc, freq, dens, AntL
;
antV = fltarr(3,3)
;print,’calling antenna_vector from XformSC
;	Antenna_Vector gives unit vectors on the right direction
antenna_vector,sc,antV 
;
;  the equation to be solved is E dot Avector = Vobs/L
;	to be solved for E,
;       so inverse of Avector * L is needed.
;  so E is invmatrix dot Vobs
;  note iA 0 = Z, 1=Y, 2=X for antennas        !!! NO MORE
;		now changed to 0=X, 1=Y, and 2=Z
;      approx 0 = R, 1 = T, and 2 = N for spacecraft directions, assuming
;      standard orientation
;
for iA = 0,2 do begin
  for iSC = 0,2 do begin
     Fmatrix[iSC,iA] = AntL[iA]*antV[iSC,iA]
  endfor
endfor
;
invmatrix = invert(Fmatrix)
;
;print,'Fmatrix ',Fmatrix
;print,'invMatx ',invmatrix
print,'f x i ',fmatrix#invmatrix
print,'i x f ',invmatrix#fmatrix
;print,'i##f  ',invmatrix##fmatrix
;antV = fltarr(3,3)
;invmatrix = -invert(Fmatrix)
;  print,' in XformSC, AntV ='
;  print,'X',antv[*,0]
;  print,'Y',antv[*,1]
;  print,'Z',antv[*,2]
 return
end
