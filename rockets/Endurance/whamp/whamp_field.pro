
function whamp_field, outs
; s= whamp_field(sol)
;compute EXB, E/B/c, Vgroup magnitude and angle for structure outs read by whamp_read_output.pro
; [x,y,z]  B background = [0.,0.,B0]
bm =[0.,0.,1.0]
n = n_elements(outs)
s = replicate( {wna:0.,wna1:0.,vg:0.,angvg:0., EXB:fltarr(3), ebc:0., ellipef:0., ellipbf:0., power:0., powerpara:0., $
    energy:0., energyes:0., energyem:0., energyblmag:0., energyelmag:0.},n)
for i=0l,n-1 do begin
    out = outs[i]
    ephase = atan(out.efl*1d0,/phase)
    exb = crossp( out.efl*1d0,out.bfl*1d0)
    s[i].wna = out.wna
    s[i].exb = sqrt( real_part(exb*conj(exb))*1d0)
    s[i].angvg = !radeg*atan(abs(out.vg[0]*1d0/out.vg[1]))
    s[i].vg = norm(out.vg)
    s[i].ebc = norm(out.efl)/norm(out.bfl)*1e6/out.c
    efl = out.efl*1e6/out.c
    psi = out.wna/!radeg
    eflelectrostatic = [efl[0]*sin(psi),efl[1]*0.0,efl[2]*cos(psi)]
    eflelectromagnetic = [efl[0]*cos(psi),efl[1],-efl[2]*sin(psi)]
;     a = norm(efl)^2
;     b = norm(eflelectrostatic)^2
;     c = norm(eflelectromagnetic)^2
;     print,'a=',a,' b+c=',b+c
;     stop
    bfl = out.bfl
    s[i].energy = norm(efl)^2 + norm(bfl)^2
    s[i].energyem = norm(eflelectromagnetic)^2 + norm(bfl)^2
    s[i].energyes = norm(eflelectrostatic)^2
    s[i].energyblmag = norm(bfl)^2
    s[i].energyelmag = norm(eflelectromagnetic)^2
    s[i].ellipbf = ellipticity_from_complex_field(out.bfl)
    s[i].ellipef = ellipticity_from_complex_field(out.efl)
    cs = fltarr(2,3)
    cs[0,*] = real_part(out.bfl)
    cs[1,*] = imaginary(out.bfl)
    ;q = kvec_ellip_phase(cs,bm,notes=notes)
    q = kvec_ellip_phase(cs,bm)
    s[i].wna1 = q.wna
    ;s[i].ellipbf = q.ellipticity
    s[i].power = q.power
    s[i].powerpara = q.powerpara
;cs[2,3] = [A,B]
endfor ;i
return,s
end
