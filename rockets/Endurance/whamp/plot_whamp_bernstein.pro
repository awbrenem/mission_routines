;Plot WHAMP results for Endurance sounding rocket Bernstein waves
;

rbsp_efw_init
!p.charsize = 2

path = '/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/whamp/'
restore, path+'whamp_w_nT4.77e+04_FRQr-3.00-10.00fcp__WNA-089__beta_7.87e-08_3.85e-06_8.74e-06_.sv',/v



wcp = sol[0].wc1
wnadeg = sol[0].wna
xc = pts_crange(101,xr=xr)
solc = cold_plasma_solutions(xc*wcp,wnadeg, plasma)
wavelengthc = 3e5/(xc*wcp/(2*!pi)*solc.r.cn)
vphasec = 3e5/solc.r.cn
;vphasec = wavelengthc*







freqs = outs2.frq
f_fch = freqs / outs2.fcp

indexrf = outs2.cindex 
;c = 2.9979250e+08


wavelength = 3e5/real_part(freqs)/indexrf   ;km
vphase =  3e5/indexrf ;km/s
efield = outs2.efl
efield2 = sol.r.efl


ellipticity = sol.r.ellipk


!p.multi = [0,0,5]
plot,f_fch,indexrf,psym=4
plot,f_fch,wavelength
plot,f_fch,wavelength2




help,outs2,/str
** Structure <4100ac08>, 33 tags, length=2064, data length=2064, refs=1:
DABS            DOUBLE           420.62536
DDXX            DOUBLE       9.0743299e-13
ICON            LONG                 1   ; 1 means solution converged.
IERR            LONG                 0
IRK             LONG                 0
ITER            LONG                 4
FRQ             DCOMPLEX  (       2176.2980,   1.5850638e-08)  ; frequency in Hz
FCP             DOUBLE           726.73683    ; proton cyclotron frequency
WC1             DOUBLE           4566.2222    ; angular proton cyclotron frequency
C               DOUBLE       2.9979250e+08    ; velocity of clight in m/s
VREF            DOUBLE           4152.3768    ; thermal velocity of of first species
CINDEX          DOUBLE           512.51573   ; index of refraction  (c*k/w)
WNA             DOUBLE           89.950000   ; wave normal angle in degrees
X               DCOMPLEX  (       2.9946163,   2.1810698e-11)   ;  w/wcp or f/fcp
Z               DOUBLE       1.8551187e-05   ; parallel kvector  k*vref/wc1
P               DOUBLE         0.021258089   ; perpendicular kvector  k*vref/wc1
D               DCOMPLEX  (      -420.62500,     -0.55204393)
DX              DCOMPLEX  (   1.5478886e+14,       371747.85)
DZ              DCOMPLEX  (  -2.3842903e+16,  -2.4284574e+08)
DP              DCOMPLEX  (  -1.9843652e+16,      -47933230.)
EFL             DCOMPLEX  Array[3]  ; set to 1 mV/M
BFL             DCOMPLEX  Array[3]
VG              DOUBLE    Array[2]
XX1             DCOMPLEX  (       2.9946163,   2.1807132e-11)
PP1             DOUBLE         0.021258089
ZZ1             DOUBLE       1.8551187e-05
AMU             DOUBLE    Array[6]
DN0             DOUBLE    Array[6]
DNS             DCOMPLEX  Array[6]  ; wave perturbation density referenced to 1mV/M E-field
VS              DCOMPLEX  Array[3, 6] ; wave perturbation velocity referenced to 1mV/M E-field
JS              DCOMPLEX  Array[3, 6] ; current density contribution from each species
XES             DCOMPLEX  Array[6, 6] ; susceptiblity tensor
EDIE            DCOMPLEX  Array[6, 4] ; diectrict tensor








help, plasma,/st
** Structure <41906148>, 6 tags, length=616, data length=612, refs=1:
DENE            DOUBLE           246920.00
DENI            DOUBLE           246920.00
BMAG            FLOAT           47669.0
SPECIES         STRUCT    -> <Anonymous> Array[6]
LABEL           STRING    'DENe=2.47e+05cm!u-3!n B!d0!n=4.77e+04 nT  0.0200 H+, 0.980 O+, 1.00 e-,'
LABELS          STRING    Array[6]

help, plasma.species[0]
** Structure <41905f98>, 9 tags, length=80, data length=80, refs=2:
NAME            STRING    'H+'
AMU             DOUBLE           1.0000000
ZCHARGE         DOUBLE           1.0000000
FRAC            DOUBLE         0.019999999
TEMPPAR         DOUBLE       9.0000001e-05
AA              DOUBLE           1.0000000
BB              DOUBLE          0.50000000
DD              DOUBLE           1.0000000
VD              DOUBLE           0.0000000

