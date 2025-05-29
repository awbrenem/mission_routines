
function ellipticity_from_complex_field,bfl
;bfl = complexarr(3)
n = 201
phi = -lindgen(n)/float(n-1)*2*!pi
bf = [0.,0.,1.] ;background filed orientation in WHAMP
bc = fltarr(3,n)

for ic=0,2 do bc[ic,*] = real_part(bfl[ic]*complex(cos(phi),sin(phi)) )

bn = sqrt(total(bc^2,1))

;plot,phi*!radeg,bn
bmin = min(bn,i0)
bmax = max(bn,i1)

db = bc[*,1]-bc[*,0]
c = crossp(bc[*,0],db) & c = c/norm(c)
dot =(total(c*bf))

;return, sgn(dot)*bmin/bmax
return, signum(dot)*bmin/bmax
end
