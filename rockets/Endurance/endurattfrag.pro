;From Henry Freudenreich. Code to rotate Endurance data to geophysical system
;Ignore the "True" stuff. 

;Rotates to w, n, u coord

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/rotation_matrices/'
if keyword_set(fil) then restore,fil else restore,path + '47001_orig_att.sav'

restore,fil else restore,path + '47001_orig_att.sav'

;if keyword_set(rots) then begin
;  print,'rotating matrices'
;  rotate_amat,a11,a12,a13,a21,a22,a23,a31,a32,a33,rots,$
;    r11,r12,r13,r21,r22,r23,r31,r32,r33
;  a11=r11 & a12=r12 & a13=r13
;  a21=r21 & a22=r22 & a23=r23
;  a31=r31 & a32=r32 & a33=r33
;endif

;if keyword_set(addt) then begin
;  time=time+addt
;  print,'adding to attitude time'
;endif
c11=interpol(a11,atime,t)
c12=interpol(a12,atime,t)
c13=interpol(a13,atime,t)
c21=interpol(a21,atime,t)
c22=interpol(a22,atime,t)
c23=interpol(a23,atime,t)
c31=interpol(a31,atime,t)
c32=interpol(a32,atime,t)
c33=interpol(a33,atime,t)
endelse
;help,c11,c12,c13,c21,t
; normalize:
r=magger(c11,c21,c31)
c11=c11/r & c21=c21/r & c31=c31/r
r=magger(c12,c22,c32)
c12=c12/r & c22=c22/r & c32=c32/r
r=magger(c13,c23,c33)
c13=c13/r & c23=c23/r & c33=c33/r

; now rotate to NWU:
if keyword_set(reverse) then begin
  b1=-by    ; E
  b2= bx    ; N
  b3= bz    ; U
  bn=c11*b1+c21*b2+c31*b3   ; x
  bw=c12*b1+c22*b2+c32*b3   ; y
  bu=c13*b1+c23*b2+c33*b3   ; z
  if keyword_set(true) then begin
    Endure_SC_to_trueSC,t,bn,bw,bu,b1,b2,b3  ; from ATTITUDE to TRUE SC coords
    bn=b1 & bw=b2 & bu=b3
  endif
endif else begin
  if keyword_set(true) then begin
    Endure_SC_to_trueSC,t,bx,by,bz,b1,b2,b3,/rev ; from TRUE to ATTITUDE SC coords
    bw=-(c11*b1+c12*b2+c13*b3)
    bn=  c21*b1+c22*b2+c23*b3
    bu=  c31*b1+c32*b2+c33*b3
  endif else  begin
    bw=-(c11*bx+c12*by+c13*bz)
    bn=  c21*bx+c22*by+c23*bz
    bu=  c31*bx+c32*by+c33*bz
  endelse
endelse
