 PRO antenna_length, sc, freq, dens,AntL;
;	returns antenna effective length in meters;	Antenna effectivee length is NEGATIVE for transform to E
;	high freq part from Oswald, "final report", avr HGA position
;	order Z = 0, Y = 1, X = 2
;	high freq changed to 63/130 times unloaded length on 18 Sep 2007
;	some problem with numbering, changed 14 Mar 2010
; version = '22 May 2017'			; low freq length set 3. m version = '26 June 2017'			; full frequency range version = '30 JUne 2027'			; removed Graz convention, A same as B, neg Leff;									;		made negative at return  B = 'B'AntL	= fltarr(3);stop, ' length, got to here 1';print,' in A_L. freq ',freq;  	AntL(2) = 1.11    AntL(1) = 1.66  	
AntL(0) = 1.35
;;	coupling as function of frequency;CB = 67.e-12CA = 62.e-12RA = (7.4177/(dens^.7) + .45)*1.e6		; approx from plot_antenna_resistance;	correct for capacitive coupling calprint,'RA ',raAntL[0:2] = ((62. + 67.)/67.)*AntL[0:2]print,'corrected ',AntL[*];complex Zant,ZbaseZant = RA/complex(1., (2.*!pi*freq*CA*RA));  this takes Rbase = infinity, only preamp input resistorZbase = 1./complex(1./33.e6, (2.*!pi*freq*CB))print,'Zbase,Zant ',Zbase,ZantZratio = abs(Zbase/(Zant+Zbase))print,'Zratio ',ZratioAntL[0:2] = -Zratio*AntL[0:2];   return
;end