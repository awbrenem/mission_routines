; .compile whamp_plot_new.pro
; apt = whamp_plot_new(sol, plasma)
function whamp_plot_new, sol, plasma
s= whamp_field(sol)
DIM=[1,1]*800
title = whamp_title(plasma, sol[0].wna )
x = real_part(sol.x)*sol.wc1/(2*!dpi)
xtitle='FREQUENCY (HZ)'
x = real_part(sol.x); *sol.wc1/(2*!dpi)
xtitle='FREQUENCY ($\Omega_P$)'
xi = IMAGINARY(sol.x)*sol.wc1/(2*!dpi)
xi = IMAGINARY(sol.x); *sol.wc1/(2*!dpi)
xr = cmb_minmax(x) ; xr[0] = 0.

wavelength = 3e5/real_part(sol.frq)/sol.CINDEX
; n = c*k/w = c/(l*f) ; w/k = c/n,  l = c/(n*f)
vphase =  3e5/sol.cindex ;wavelength*x

wcp = sol[0].wc1
wnadeg = sol[0].wna
;xc = pts_crange(101,xr=xr)
;solc = cold_plasma_solutions(xc*wcp,wnadeg, plasma)
;wavelengthc = 3e5/(xc*wcp/(2*!pi)*solc.r.cn)
;vphasec = 3e5/solc.r.cn
;;vphasec = wavelengthc*


ch = 0.02
nx = 1
ny = 5
pos = [.15,.1,.85,.85]
ip =0
;apt = plot(x, sol.ddxx,/ylog, ytitle='DDXX', xr=xr)
help, x, xi
symbol='circle'
sym_size =0.5
sym_filled=1
posa = sab_pos(nx,ny,ch=ch,pos=pos,ip)
ytitle='GR (Hz)'
ytitle='$GR (\Omega_P)$'
; , symbol=symbol, sym_filled=sym_filled, sym_size=sym_size
apt = scatterplot(x, xi, ytitle=ytitle, title=title, xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, $
	 current=0, pos = posa, dim=dim)
apt0 =plot(/over,xr,xr*0, linestyle=2)
ylog=1
posa = sab_pos(nx,ny,ch=ch,pos=pos,ip)
apt = scatterplot(ylog=ylog, x, wavelength, ytitle='wavelength (km)', xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, $
	 current=1, pos =posa)
;apt0 = plot(/over,xc,wavelengthc,'blue', thick=2)

posa = sab_pos(nx,ny,ch=ch,pos=pos,ip)
apt0 = scatterplot(/ylog, x, s.ELLIPEF, ytitle='ELLIPTICTY E-FIELD', xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, current=1, pos =posa)

posa = sab_pos(nx,ny,ch=ch,pos=pos,ip)
efl = abs(sol.efl)
apt01 = scatterplot(/ylog, x, reform(efl[0,*]), ytitle='|E-FIELD|', xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, sym_color='red', name='|Ex|', current=1, pos =posa)
apt02 = scatterplot(over=apt01, x, reform(efl[1,*]), xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, sym_color='blue', name='|Ey|')
apt03 = scatterplot(over=apt01, x, reform(efl[2,*]), xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, sym_color='green', name='|Ez|')
leg = LEGEND(TARGET=[apt01,apt02,apt03], POSITION=[.8,.9], /REL, /AUTO_TEXT_COLOR, HORIZONTAL_ALIGNMENT=0, vertical_ALIGNMENT=1.)

;apt = scatterplot(ylog=ylog, x, vphase, ytitle='Vphase (km/s)', xtitle='', xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, $
;	 current=1, pos =posa)
;apt0 = plot(/over,xc,vphasec,'blue', thick=2)


is=where( plasma.species.name eq 'e-' and plasma.species.frac gt 0)
dns =abs(sol.dns) ;& dns = total(reform(dns[is,*]),1)/total(sol.dn0[is])
dns = reform(dns[is,*])/total(sol.dn0[is])
posa = sab_pos(nx,ny,ch=ch,pos=pos,ip)
apt = scatterplot(ylog=ylog, x, dns, ytitle='$\deltaN_e/N_e$', xtitle=xtitle, xr=xr, symbol=symbol, sym_filled=sym_filled, sym_size=sym_size, $
	 current=1, pos =posa)

fout = 'whamp_'+ string(wnadeg, format='(i3.3)')
sxr =strtrim(xr,2)
fout = fout + '_' + sxr[0] + '_' + sxr[1]
stop
apt.Save, fout + ".png", BORDER=10, RESOLUTION=300 ;, /TRANSPARENT


return, apt
end

pro test_whamp_plot_new

if n_elements(frq) eq 0 then begin
frqr=[.5,10.]
;nfrq = (frqr[1]-frqr[0])/0.1 + 1
nfrq = (frqr[1]-frqr[0])/0.02 + 1
frq = lindgen(nfrq)/float(nfrq-1)*(frqr[1]-frqr[0]) + frqr[0]
;wnadeg = 89.90d0
;wnadeg = 85.0d0
wnadeg = 80.0d0
whamp_w_v0,frqr=frqr,wnadeg = wnadeg,/fixk, sol=sol, nfrq=nfrq
endif

files = file_search('whamp_w*.sv')
string_list, files
read,'input file #:', ifile
file = files[ifile]
restore,file,/v
sol = outs2

apt = whamp_plot_new(sol, plasma)

stop
end
