


path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/WHAMP/'

fn = 'boardsen_run_mar26_2024.dat'

;Cold plasma test run:

;98% O+ and 2% H+.
;ne = 246920
;ni = 246920
;bmag = 47669

;dns in M^-3;
;vs in km/s
;EFL in mV/M
;|BFL| is 1 nT
;vs[ic,ipar]  velocity of species S in m/s
;js[*,ipar]  current density of species S in A/m^2
;J = total(js,2) in A/m^2
;group velocity Vg in in units of the speed of light c


IDL> print,sol.plasma.labels[0]
0.0200 H+ T=9.00e-05 keV A=1.00 B=0.500 D=1.00 V=0.00 m/s
IDL> print,sol.plasma.labels[1]
0.980 O+ T=9.00e-05 keV A=1.00 B=0.500 D=1.00 V=0.00 m/s
IDL> print,sol.plasma.labels[2]
1.00 e- T=0.000200 keV A=1.00 B=0.500 D=1.00 V=0.00 m/s



;with ellipticity ~ 0.77
;Only RH solutions available - no LH solutions.  

restore,path+fn

help,sol.plasma.species[1],/st

help,sol,/st
* Structure <5016aa08>, 6 tags, length=111168, data length=111152, refs=1:
R               STRUCT    -> <Anonymous> Array[101]
L               STRUCT    -> <Anonymous> Array[101]
COEF            STRUCT    -> <Anonymous> Array[1]
PLASMA          STRUCT    -> <Anonymous> Array[1]
PAR             STRUCT    -> <Anonymous> Array[3]



;*****************
;Range of wave frequencies (101 of them). The first is at 5000 Hz

findex = 1

IDL> help,sol.r[findex],/st
** Structure <507bb408>, 14 tags, length=504, data length=504, refs=2:
  W               DOUBLE           4999.5307
   PSI             DOUBLE          0.78539816
   CN              DOUBLE           157.01654
   VG              DOUBLE        0.0075783181
   VG_PARA         DOUBLE        0.0071369968
   VG_PERP         DOUBLE        0.0025483686
   ELLIPK          DOUBLE          0.94813895
   EX_EY           DOUBLE           1.3413967
   BPARA_B         DOUBLE          0.51312860
   EFL             DCOMPLEX  Array[3]
   BFL             DCOMPLEX  Array[3]
   DNS             DCOMPLEX  Array[3]
   VS              DCOMPLEX  Array[3, 3]
   JS              DCOMPLEX  Array[3, 3]

sol.coef.R[findex]
;13422.603710862149
sol.coef.L[findex]
;-15317.993830257337
sol.coef.P[findex]
;-17767306.194819409
sol.coef.S[findex]
;-947.69505969759393
sol.coef.D[findex]
;14370.298770559744
sol.coef.RW[findex]
;-1.9242118482649648
sol.coef.LW[findex]
;2.6985148664981189
sol.coef.PW[findex]
;5342.9471242091740
sol.coef.SW[findex]
;0.38715150911657720
sol.coef.DW[findex]
;-2.3113633573815422

;index ~ 60 







IDL> help,x,/st
** Structure <4f767208>, 36 tags, length=1512, data length=1500, refs=1:
   FREQ            FLOAT           5000.00
   KMAG            DOUBLE           4.9951942
   THETA_KB        FLOAT              -NaN
   N               FLOAT           47.7006
   DENS            FLOAT           246920.
   BO              FLOAT           47669.0
   WAVELENGTH      DOUBLE           1.2578461
   RESANGLE        FLOAT           90.0000
   PHASEVEL        DOUBLE           6289.2306
   BPOL            FLOAT              -NaN
   CYCLO_COUNTERSTREAM_RES
                   STRUCT    -> <Anonymous> Array[1]
   CYCLO_COSTREAM_RES
                   STRUCT    -> <Anonymous> Array[1]
   LANDAU_RES      STRUCT    -> <Anonymous> Array[1]
   ENERGY_NOTES    STRING    Array[4]
   FPE             FLOAT       4.46226e+06
   FPH             FLOAT           14727.6
   FPHE            FLOAT           0.00000
   FPO             FLOAT           25773.4
   FCH             FLOAT           726.978
   FCHE            FLOAT           181.745
   FCO             FLOAT           45.4361
   FLHR            FLOAT           8506.71
   E_INERTIAL      FLOAT         0.0106860
   ION_INERTIAL    FLOAT          0.458835
   VA_ION          FLOAT           2091.29
   VA_ELECTRON     FLOAT           89608.8
   MEPP            FLOAT           22885.3
   EX2EY           FLOAT          0.770000
   EZ2EX           FLOAT              -NaN
   CB2E_OBLIQUEPROPAGATION_V1
                   FLOAT              -NaN
   CB2E_OBLIQUEPROPAGATION_V2
                   FLOAT              -NaN
   CB2E_PARALLELPROPAGATION
                   FLOAT              -NaN
   BETAANGLE_E_K   FLOAT              -NaN
   CP_PARAMS       STRUCT    -> <Anonymous> Array[1]
   NOTES           STRING    Array[11]


AMU             DOUBLE       0.00054462998
ZCHARGE         DOUBLE          -1.0000000
FRAC            DOUBLE           1.0000000
DEN             DOUBLE           246920.00
NAME            STRING    'e-'
FP2_X           DOUBLE    Array[3]
FC_X            DOUBLE    Array[3]
NV              LONG              2001
VRPARA          DOUBLE    Array[2]
VRPERP          DOUBLE    Array[2]
TEMPPAR         DOUBLE       0.00019999999
VTH             DOUBLE           264449.84
AA              DOUBLE           1.0000000
BB              DOUBLE          0.50000000
DD              DOUBLE           1.0000000
VD              DOUBLE           0.0000000





