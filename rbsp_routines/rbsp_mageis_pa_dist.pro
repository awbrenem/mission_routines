;Code comes from flux_calculation_jan6.pro


;SOURCES OF ERROR
;-MagEIS extrapolation to 3 deg - could be big
;-BARREL background subtraction
;-Isotropic assumption in Bremsstrahlung model. John Sample says this could make BARREL fluxes 3x too low. 
;-Uncertainties on BARREL fit parameters which, according to Michael McCarthy
;aren't small
;-Hiss freq spec modeled with Gaussian. 


;Calculate the flux of e- scattered into loss cone on RBSP

;Jan6th, 2014 at 21:00 and just before  (MagEIS df = 2-4 x 10^2


;Here's what I'm doing. By knowing the pitch angle diffusion rate I can find how far in
;pitch angle an electron will scatter in a single bounce period. (Half of) Electrons within this 
;angular distance from the loss cone will be lost in this time. If I take into account the
;entire 2*Pi solid angle, all electrons within this "thin slice" are lost. 


;Load MagEIS L3 file and find out the PA distribution
path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events_campaign2/Jan10/'
fn = 'rbspa_rel04_ect-mageis-L3_20140110_v8.1.0.cdf' & cdf2tplot,path+fn
;fn = 'rbspa_rel04_ect-mageis-L3_20140111_v8.1.0.cdf'
;fn = 'rbspb_rel04_ect-mageis-L3_20140110_v8.1.0.cdf'
;fn = 'rbspb_rel04_ect-mageis-L3_20140111_v8.1.0.cdf'


get_data,'FEDU',data=dd,dlim=dlim,lim=lim


;It takes RBSP-A about 15 min to cross the BARREL field of view (~1 RE
;mapped to magnetic equator). 
t0 = time_double('2014-01-10/21:50:30')
t = time_double('2014-01-10/21:58')
t1 = time_double('2014-01-10/22:05:30')

goo = where(dd.x ge time_double(t))
print,time_string(dd.x[goo[0]])
;2014-01-10/21:58:10



;--------------------------------------------------
;Determine range of MagEIS df at the loss cone
;--------------------------------------------------

;(1)--------------------------------------------------
;PLOT THE PA DISTRIBUTION FROM MAGEIS L3
;--------------------------------------------------

;MagEIS for time t; channel2 = 50 keV (from Blake13)
for i=0,10 do print,'PA= ',dd.v1[i],' : counts= ',dd.y[goo[0],2,i] ;counts for all pitch angles at 20:58:02 for 50 keV
;PA=        8.1818182 : counts=        4665.5674
;PA=        24.545455 : counts=        6794.9805
;PA=        40.909091 : counts=        8153.9766
;PA=        57.272727 : counts=        9130.0942
;PA=        73.636364 : counts=        9934.0703
;PA=        90.000000 : counts=        9187.5801
;PA=        106.36364 : counts=        9935.3261
;PA=        122.72727 : counts=        8785.6335
;PA=        139.09091 : counts=        7680.2422
;PA=        155.45455 : counts=        7704.1685
;PA=        171.81818 : counts=        5742.2964


;Plot range of MagEIS values vs PA for times t0-t1
 
plot,dd.v1,dd.y[goo[0],2,*],/ylog,yrange=[4000,12000],ystyle=1,ytitle=dlim.ysubtitle
for i=0,42 do oplot,dd.v1,dd.y[goo[0]-i,2,*],color=i
for i=0,42 do oplot,dd.v1,dd.y[goo[0]+i,2,*],color=i

;Values from this plot range from about 75 (same as Michael's
;estimate) to few hundred.


;(2)--------------------------------------------------
;JOE FENNELL'S 15 MIN AVERAGED FIT CENTERED ON 20:58
;--------------------------------------------------

;Joe Fennell has averaged the MagEIS flux over 15 minutes centered on
;this time. He then extracted the values at 4 deg PA using a fit
;function. For four deg loss cone I get a MagEIS value of 
;; m1 = 759.2   ;pm 100.46
;; m2 = 0.39326   ;pm 0.168
;; m3 = -1770.6  ;pm 105.79
;; m4 = 11.241   ;pm 1.542

m1 = [658.7,759.2,860]
m2 = [0.561,0.393,0.225]

m3 = [-1664.8,-1770.6,-1876.4]
m4 = [9.7,11.241,12.78]   


;x = indgen(180)*!dtor
x = 4.*!dtor
y = m1*sin(x)^(m2) - m3*sin(x)^(m4)
;; RBSP_EFW> print,y
;;       147.890      266.615      472.395


;(3)--------------------------------------------------
;Michael's fit of flux=k sin(alpha) at 20:58. Note that at
;this time the PA distribution is more perp than surrounding
;times. From his email, df = 75
;--------------------------------------------------

;(4)--------------------------------------------------
;Uncertainty in the functional form of PA distribution near loss cone
;This could be a very large error as there's hardly any
;experimental info on fine details of pitch angle distributions near
;the loss cone
;--------------------------------------------------








;Start with MagEIS differential flux at 3 deg PA, 50 keV
;See above comments (1)-(4) regarding this range of reasonable values

;df = [75,266,400]                  ;e-/cm2/s/sr/keV  (for 4 deg PA)
df = [25,710]




;Now find the "thin slice" solid angle term.
;To get this I need to know how far an electron will random walk in a single
;bounce period.

                                ;bounce-averaged diffusion coeff for 50 keV
;        Daa = 3d-4  ;Jan6 at 20:58
DaaA = 1d-4                     ;Jan6 at 20:58 but averaged over 15 min centered at this time. 


;Range of L-shells from t0 to t1
L0 = 4.753
L1 = 4.432


;Range of Bmag from t0 to t1  (mlat within 1 deg of mag eq)  
Bo_mageq0 = 257.
Bo_mageq1 = 309.


;But, we have to account for the enhanced flux caused by 
;conservation of flux in narrowing flux tube
;Find ratio of B2/B1 = A1/A2. (Note that 70 km used here may be a bit
;low for 50 keV electrons. Probably 100 km is a better value.)


.compile ~/Desktop/community/dipole.pro
dip0 = dipole(L0)
dip1 = dipole(L1)

radius0 = dip0.r - 6370.
radius1 = dip1.r - 6370.

boo = where(radius0 le 70.)
boo = boo[0]
Bo_70km0 = dip0.b[boo]          ;Bo_70km0 = 53282 nT
boo = where(radius1 le 70.)
boo = boo[0]
Bo_70km1 = dip1.b[boo]          ;Bo_70km1 = 52964 nT

                                ;Bo_mageq = 160.   ;nT  at 20:05 UT 
Bratio0 = Bo_70km0/Bo_mageq0    ;207
Bratio1 = Bo_70km1/Bo_mageq1    ;171



;Find loss cone angle
ang0 = asin(sqrt(1/Bratio0))/!dtor ;3.99 deg
ang1 = asin(sqrt(1/Bratio1))/!dtor ;4.38 deg




                                ;bounce period for 30 keV e- at loss cone
Ekev = 50.
;; Tb0 = 5.62d-2 * L0 * (1-0.43*sin(ang0*!dtor))/sqrt(EkeV/1000.)         ;=1.16 s
;; Tb1 = 5.62d-2 * L1 * (1-0.43*sin(ang1*!dtor))/sqrt(EkeV/1000.)         ;=1.08 s

;Manually set bounce period to 1 sec for calculation convenience. This
;is a very small source of error compared to other terms. 
Tb0 = 1.
Tb1 = 1.


;Random walk distance in single bounce <Daa> = dangle^2/2tb  (Kennel69 eqn3)
;	dangle = sqrt(Daa*2*Tb)/!dtor  ;=0.8 deg
dangleA0 = sqrt(DaaA*2*Tb0)/!dtor         ;=0.87 deg
dangleA1 = sqrt(DaaA*2*Tb1)/!dtor         ;=0.84 deg

                                ;"thin slice" solid angle
;	SA = 2*!pi*sin(ang*!dtor) * dangle*!dtor 
SAA0 = 2*!pi*sin(ang0*!dtor) * dangleA0*!dtor 
SAA1 = 2*!pi*sin(ang1*!dtor) * dangleA1*!dtor 


;Now I can solve for the #/s/cm2/keV of electrons lost (from MagEIS)
;	nm = df * SA   
nmA0 = df * SAA0 
nmA1 = df * SAA1   


;Of these electrons in the "thin slice" half will random walk to the DLC and the other
;half to higher pitch angles. So, 
;	nm = nm/2. ;=242 e-/s/cm2/keV
nmA0 = nmA0/2.                  ;=251 e-/s/cm2/keV
nmA1 = nmA1/2.                  ;=242 e-/s/cm2/keV


;BARREL e-/s/cm2 50 keV (from Michael McCarthy's email)
nb = 13                         ;at 20:58 (MMcCarthy's email)
;nb = [25,1]
;Lots of uncertainty associated with this number



;--------------------------------------------------
;Compare BARREL and RBSP observations
;--------------------------------------------------

                                ;What should be observed by BARREL is the MagEIS flux
                                ;enhanced by Bratio due to field line
                                ;focusing effect

                                ;nb2 = nm * Bratio   
nb2A0 = nmA0 * Bratio0
nb2A1 = nmA1 * Bratio1


                                ;print,nb2/nb 
print,nb2A0/nb
print,nb2A1/nb
print,'*****'
                                ;adjust for non-isotropic flux (John
                                ;Sample indicates that this is likely
                                ;a factor of 3)
                                ;print,nb2/nb/3.
print,nb2A0/nb/3.
print,nb2A1/nb/3.

;; RBSP_EFW> print,nb2A1/nb/3.
;;        1.1185967       3.9672897       5.9658491

;;      0.37286557       10.589382

stop
end
