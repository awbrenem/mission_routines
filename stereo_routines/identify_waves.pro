;+
;*****************************************************************************************
;
;   FUNCTION :   identify_waves
;   PURPOSE  :   Identifies all the wave from the input structure (from get_TDS.pro). Adds
;				 wave identification to the input structure.
;
;
;	;------------
;	;DUST IMPACTS
;	;------------
;	; freq - typically 50 to >1000 Hz
;	; bandwidth - usually large. Can get as low as ~10, but often >1000
;	; bandwidth stddev - usually from 1-10. Can be <1
;	; df/f - almost always 1-5, can be < 1
;	; quality = -1
;	; amp -> few to a few hundred mV/m
;
;	;---------------------------
;	;non-identified dust impacts
;	;---------------------------
;	;Sometimes dust impacts are not identified by wave_quality.pro
;	;An easy eay to identify them is with the following tests.
;	;NOTE THAT WAVE_QUALITY.PRO CANNOT TEST FOR THIS KIND OF DUST B/C 
;	;IT ONLY LOOKS AT A SINGLE CHANNEL AT A TIME.
;
;	;------------------------
;	;broadband whistler waves
;	;------------------------

	;df_f < 1
	;qual < 100
	;freq < 500
	;bw < 100

;	;--------------------
;	;narrowband whistlers
;	;--------------------

	;df_f < 1
	;qual > 100
	;freq < 500
	;bw < 100


;	;------------------------------------------------------------------------------------------
;	;noise/choppy waves - essentially the same as broadband whistlers but with a bandwidth stddev > 1. 
;	;		These are less like whistler mode time-domain signals and more like noise.
;	;		Some waves identified as noise may actually be ion acoustic or Langmuir waves.
;	;-------------------------------------------------------------------------------------------

	;f < 10000
	;qual < 50
	;df_f > 1


;	;--------------------------------
;	;Langmuir waves
;	;--------------------------------

;If "CheckPlasma" keyword is set:
;	*Langmuir f>fpe  and  Ion Acoustic f<fpe
;	*If Plasma data is unavailable it classifies as 'noise' 

	;df_f < 1
	;f > 10000
	;qual > 100

;	;--------------------------------
;	;Ion Acoustic waves
;	;--------------------------------

;If "CheckPlasma" keyword is set:
;	*Langmuir f>fpe  and  Ion Acoustic f<fpe
;	*If Plasma data is unavailable it classifies as 'noise' 

	;200 < freq < 10000
	;bw > 100
	;quality < 1000



            
;
;   CALLS:	get_TDS.pro
;               
;   REQUIRES:   the output file from wave_quality.pro 
;               
;   INPUT:		
;   EXAMPLES:    wave_struct = identify_waves(t0,t1,file=wave_file,/checkplasma)
;               
;
;   KEYWORDS:   t0 -> initial time to consider. Defaults to '2007-01-01/00:00:00'
;				t1 -> final time to consider. Defaults to '2100-01-01/00:00:00'
;				file -> file to read. Output from wave_quality.pro
;				checkplasma -> compare wave freq to plasma freq to help identify difference
;						b/t Langmuir and ion acoustic waves. If not set then 
;						Langmuirs have freq >10 kHz and Ion Acoustics <10 kHz
;
;
;   CHANGED:  
;
;   NOTES:    
;
;STEPS
;	1. IDENTIFY QUAL=0. THESE ARE THRUSTER EVENTS AND CALIBRATION EVENTS. 
;	2. IDENTIFY DUST IMPACTS, Q=-1 (WAVE_QUALITY.PRO SEEMS TO BE DOING A VERY GOOD JOB OF IDENTIFYING THESE)
;	3. IDENTIFY MISSED DUST IMPACTS. DONE BY COMPARING B/T CHANNELS, SOMETHING I CANNOT DO IN WAVE_QUALITY.PRO
;	4. IDENTIFY ALL WHISTLERS
;	5. IDENTIFY LARGE AMPLITUDE, NARROWBAND WHISTLERS
;	6. IDENTIFY LANGMUIR WAVES
;	7. IDENTIFY ACOUSTIC WAVES
;
;
;   CREATED:  01/31/2011   v1.0.0
;   CREATED BY:  Aaron W. Breneman
;
;    LAST MODIFIED:  07/28/2011   v1.1.0
;    MODIFIED BY: Sam Schreiner
;    LAST MODIFIED:  11/28/2011   v1.1.0
;    MODIFIED BY: Aaron Breneman
;
;	-Changed method for determining Langmuirs and Ion Acoustics by
;	comparing plasma frequency (obtained from density data) to wave frequency.
;
;*****************************************************************************************
;-
;*****************************************************************************************


function identify_waves,t0,t1,file=file,checkPLASMA=checkPLASMA

init_path	;set up general path to data

if not keyword_set(t0) then t0 = '2007-01-01/00:00:00'
if not keyword_set(t1) then t1 = '2100-01-01/00:00:00'
if not keyword_set(file) then file = !data.stereo.wavequality + 'STEREO_TDS_STA_q_waveclass_20070101_to_20070515_OLD.txt'

tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1

tmp0 = t0
tmp1 = t1

;.compile get_TDS.pro	;sometimes necessary?
x = get_TDS(tmp0,tmp1,file=file)

type = strarr(n_elements(x[0].times),4)  ;wave type for each burst capture
type2 = type

;--------------------------------------
;CALIBRATION STATE OR THRUSTER CAPTURES
;--------------------------------------
tx = where(x[0].quality eq 0.)
ty = where(x[1].quality eq 0.)
tz = where(x[2].quality eq 0.)
tp = where(x[3].quality eq 0.)

thruster = [tx,ty,tz,tp]


;----------------------------------
;previously identified dust impacts
;----------------------------------
;Captures with quality=-1 are reliably dust. However, wave_quality.pro does miss many dust impacts.

dx = where(x[0].quality eq -1)
dy = where(x[1].quality eq -1)
dz = where(x[2].quality eq -1)
dp = where(x[3].quality eq -1)


;---------------------------
;non-identified dust impacts
;---------------------------

;Sometimes dust impacts are not identified by wave_quality.pro
;An easy way to identify them is with the following tests.
;NOTE THAT WAVE_QUALITY.PRO CANNOT TEST FOR THIS KIND OF DUST B/C 
;IT ONLY LOOKS AT A SINGLE CHANNEL AT A TIME.


;Compare the amp of signal in one direction to the amps of the other two
;directions. If this signal is >5x larger than it is almost certainly dust.

;dust in x-antenna
amprat_xy = x[0].amp/x[1].amp
amprat_xz = x[0].amp/x[2].amp

dx2 = -1 & dy2 = -1 & dz2 = -1

x0 = where((amprat_xy ge 5.) or (amprat_xz ge 5.))  
x1 = where(x[1].quality ne -1.) ;make sure it's not already classified as dust
if x0[0] ne -1 and x1[0] ne -1 then dx2 = setintersection(x0,x1)

;dust in y-antenna
amprat_yx = x[1].amp/x[0].amp
amprat_yz = x[1].amp/x[2].amp

y0 = where((amprat_yx ge 5.) or (amprat_yz ge 5.))  
y1 = where(x[1].quality ne -1.) ;make sure it's not already classified as dust
if y0[0] ne -1 and y1[0] ne -1 then dy2 = setintersection(y0,y1)

;dust in z-antenna
amprat_zx = x[2].amp/x[0].amp
amprat_zy = x[2].amp/x[1].amp

z0 = where((amprat_zx ge 5.) or (amprat_zy ge 5.))  
z1 = where(x[0].quality ne -1.) ;make sure it's not already classified as dust
if z0[0] ne -1 and z1[0] ne -1 then dz2 = setintersection(z0,z1)

;dust in pseudo if we have dust in Ex or Ey
dp2 = [dx2,dy2]


;	;------------------------
;	;broadband whistler waves
;	;------------------------
	;df_f < 1
	;qual < 100
	;freq < 500
	;bw < 100


tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1

wh = make_array(n_elements(x[0].freq),4,/LONG,value=-1)
FOR ii=0,3 DO BEGIN
  x1 = where(x[ii].freq le 500.)        
  x2 = where(x[ii].bandwidth le 100.)
  x3 = where((x[ii].quality gt 0.) and (x[ii].quality lt 100.))   
  x4 = where(x[ii].deltaf_f le 1)

  if x1[0] ne -1 and x2[0] ne -1 then tmp = setintersection(x1,x2)
  if tmp[0] ne -1 and x3[0] ne -1 then tmp = setintersection(tmp,x3)
  if tmp[0] ne -1 and x4[0] ne -1 then tmp = setintersection(tmp,x4)
  
  wh[0:n_elements(tmp)-1,ii] = tmp
  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1
ENDFOR


;	;--------------------
;	;narrowband whistlers
;	;--------------------
	;df_f < 1
	;qual > 100
	;freq < 500
	;bw < 100


law = make_array(n_elements(x[0].freq),4,/LONG,value=-1)
FOR ii=0,3 DO BEGIN
  x1 = where(x[ii].freq le 500.)        
  x2 = where((x[ii].bandwidth ge 0.) and (x[ii].bandwidth le 100.))
  x3 = where((x[ii].quality ge 100.) and (x[ii].quality le 1000.))   
  x4 = where(x[ii].deltaf_f le 1)
  x5 = where(x[ii].amp ge 1)

  if x1[0] ne -1 and x2[0] ne -1 then tmp = setintersection(x1,x2)
  if tmp[0] ne -1 and x3[0] ne -1 then tmp = setintersection(tmp,x3)
  if tmp[0] ne -1 and x4[0] ne -1 then tmp = setintersection(tmp,x4)
  if tmp[0] ne -1 and x5[0] ne -1 then tmp = setintersection(tmp,x5)
  law[0:n_elements(tmp)-1,ii] = tmp
  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1
ENDFOR


;------------------------------------------------------------------------------------------
;noise/choppywaves - essentially the same as broadband whistlers but with a bandwidth stddev > 1. 
;		These are less like whistler mode time-domain signals and more like noise.
;		Some waves identified as noise may also be identified as Langmuir waves. Therefore
;-------------------------------------------------------------------------------------------
	;f < 10000
	;qual < 50
	;df_f > 1


n = make_array(n_elements(x[0].freq),4,/LONG,value=-1)
FOR ii=0,3 DO BEGIN
  x2 = where(x[0].deltaf_f ge 1.)
  x3 = where(x[0].quality lt 50.)
  x4 = where(x[0].freq le 10000.)


  if x2[0] ne -1 and x3[0] ne -1 then tmp = setintersection(x2,x3)
  if tmp[0] ne -1 and x4[0] ne -1 then tmp = setintersection(tmp,x4)

  n[0:n_elements(tmp)-1,ii] = tmp
  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1
ENDFOR




;	;--------------------------------
;	;Langmuir waves
;	;--------------------------------

;If "CheckPlasma" keyword is set:
;	*Langmuir f>fpe  and  Ion Acoustic f<fpe
;	*If Plasma data is unavailable it classifies as 'noise' 

	;df_f < 1
	;f > 10000
	;qual > 100

lang = make_array(n_elements(x[0].freq),4,/LONG,value=-1)
FOR ii=0,3 DO BEGIN
  x1 = where((x[ii].freq gt 10000.))        
  x2 = where((x[ii].bandwidth ge 100.))
  x3 = where((x[ii].quality ge 100.))   
  x4 = where(x[ii].deltaf_f lt 1)
  x5 = where(x[ii].amp ge 1)


  if x1[0] ne -1 and x2[0] ne -1 then tmp = setintersection(x1,x2)
  if tmp[0] ne -1 and x3[0] ne -1 then tmp = setintersection(tmp,x3)
  if tmp[0] ne -1 and x4[0] ne -1 then tmp = setintersection(tmp,x4)
  if tmp[0] ne -1 and x5[0] ne -1 then tmp = setintersection(tmp,x5)
  lang[0:n_elements(tmp)-1,ii] = tmp
  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1
ENDFOR


;	;--------------------------------
;	;Ion Acoustic waves
;	;--------------------------------

;If "CheckPlasma" keyword is set:
;	*Langmuir f>fpe  and  Ion Acoustic f<fpe
;	*If Plasma data is unavailable it classifies as 'noise' 

	;200 < freq < 10000
	;bw > 100
	;quality < 1000


iaw = make_array(n_elements(x[0].freq),4,/LONG,value=-1)
FOR ii=0,3 DO BEGIN
  x1 = where((x[ii].freq ge 200.) and (x[ii].freq le 10000.))        
  x2 = where((x[ii].bandwidth ge 100.))
  x3 = where((x[ii].quality le 1000.))   
  x5 = where(x[ii].amp ge 1)

  if x1[0] ne -1 and x2[0] ne -1 then tmp = setintersection(x1,x2)
  if tmp[0] ne -1 and x3[0] ne -1 then tmp = setintersection(tmp,x3)
  if tmp[0] ne -1 and x5[0] ne -1 then tmp = setintersection(tmp,x5)
  iaw[0:n_elements(tmp)-1,ii] = tmp
  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1 & x4 = -1 & x5 = -1 & x6 = -1
ENDFOR


;;	;-----------------------------
;;	;Dust captures
;;	;-----------------------------
;
;dust = make_array(n_elements(x[1].freq),4,/LONG,value=-1)
;FOR ii=0,3 DO BEGIN
;
;  x1 = where(x[ii].quality eq -1)
;;  x1 = where((x[ii].freq ge 200.) and (x[ii].freq le 10000.))        
;;  x2 = where((x[ii].bandwidth ge 100.))
;;  x3 = where((x[ii].quality le 1000.))   
;;  x5 = where(x[ii].amp ge 1)
;
;  tmp = x1
;  dust[0:n_elements(tmp)-1,ii] = tmp
;  tmp = -1 & x1 = -1
;ENDFOR
	



;;########################################################################
;;---------------------------------------
;; Langmuir and Ion Acoustic waves (new in Version 1.1.0)
;;---------------------------------------
;; bandwidth > 1000 Hz
;; bandwidth stddev - ranges from 0 to >100 (PROBS THE BEST MEASURE FOR LANGMUIR WAVES)
;; df/f ~ 1. Interesting, this is remarkably consistent - often ~1.23 even though the range of freq and bw vary
;;			Can sometimes be < 1
;; quality > 0
;; amp -> few to a hundred mV/m 0-p
;
;
;;find the indicies of the combined langmuir/ion acoustic waves
;langORion = make_array(n_elements(x[1].bandwidth),4,/LONG,value=-1)
;lang = langorion & ion = langorion
;
;for ii=0,3 do begin
;  x1 = where(x[ii].bandwidth ge 1000.)
;  x2 = where(x[ii].quality gt 0.)
;  x3 = where(x[ii].bandwidth_stddev ge 0.)
;
;  if x1[0] ne -1 and x2[0] ne -1 then tmp = setintersection(x1,x2)
;  if tmp[0] ne -1 and x3[0] ne -1 then tmp = setintersection(tmp,x3)
;  langORion[0:n_elements(tmp)-1,ii] = tmp
;
;  tmp = -1 & x0 = -1 & x1 = -1 & x2 = -1 & x3 = -1
;endfor
;
;
;if keyword_set(checkPlasma) then begin	
;;---------- Langmuir f>fpe  and  Ion Acoustic f<fpe ----------
;
;  pla_times = [''] &   pla_n = [0.0]
;
;  year0 = strmid(min(x.times),0,4)
;  year1 = strmid(max(x.times),0,4)
;  num_years = fix(year1)-fix(year0)+1
;
;  for kk=0,num_years-1 do begin		;read in plasma density data
;    pla_file = !data.stereo.wavequality + '/SWEA_mom_STA_at_burst_times_'+ strtrim(string(fix(year0)+kk),1) +'_complete.txt'
;    pla_struct = get_sw_data('mom',file=pla_file,/removebad)
;    pla_times = [pla_times,pla_struct.burstTS]
;    pla_n = [pla_n,pla_struct.n]
;  endfor
;
;
;  ind0 = where(pla_times eq x[1].times[0])	;trim density data to the same length as burst data
;  ind1 = where(pla_times eq x[1].times[n_elements(x[1].times)-1])
;  pla_times = pla_times[ind0:ind1]
;  pla_n = pla_n[ind0:ind1]
;  fpe = 8980 * sqrt(pla_n)	;plasma frequency
;
;
;  for ii=0,3 do begin
;    good_plasma = where(finite(fpe) eq 1 and fpe gt 0)	;no events with bad plasma data or zero density
;;     lm = where(abs((x[ii].freq - fpe))/fpe le .1)
;;     ia = where(abs((x[ii].freq - fpe))/fpe gt .1)
;    lm = where(x[ii].freq ge fpe)
;    ia = where(x[ii].freq lt fpe)
;    lm2 = setintersection(lm,good_plasma)
;    ia2 = setintersection(ia,good_plasma)
;    lm3 = setintersection(lm2,langORion[*,ii])
;    ia3 = setintersection(ia2,langORion[*,ii])
;
;    lang[0:n_elements(lm3)-1,ii] = lm3
;    ion[0:n_elements(ia3)-1,ii] = ia3
;
;    lm=-1 & & ia=-1 & lm2=-1 & ia2=-1 & lm3=-1 & ia3=-1
;  endfor
;
;  lm=-1 & & ia=-1 & lm2=-1 & ia2=-1 & lm3=-1 & ia3=-1
;
;endif else begin
;;---------- Langmuir f>10 kHz  and  Ion Acoustic f<10 kHz ----------
;
;  for ii=0,3 do begin
;    lm = where(x[ii].freq ge 10000.)	;langmuir
;    ia = where(x[ii].freq lt 10000.)	;ion acoustic
;
;    if lm[0] ne -1 and langORion[0,ii] ne -1 then tmp = setintersection(lm,langORion[*,ii])
;    lang[0:n_elements(tmp)-1,ii] = tmp
;    if ia[0] ne -1 and langORion[0,ii] ne -1 then tmp = setintersection(ia,langORion[*,ii])
;    ion[0:n_elements(tmp)-1,ii] = tmp
;    
;    tmp=-1 & lm=-1 & ia=-1
;  endfor
;
;endelse
;########################################################################




;--------------------------------
;FIND WHERE THERE IS MISSING DATA
;--------------------------------

mx1 = where(finite(x[0].freq) eq 0.)
my1 = where(finite(x[1].freq) eq 0.)
mz1 = where(finite(x[2].freq) eq 0.)
mp1 = where(finite(x[3].freq) eq 0.)


;--------------------------
;NOW CLASSIFY ALL THE WAVES
;--------------------------
;NOTE: THE ORDER YOU PUT THESE IN IS IMPORTANT SINCE SOME WAVES CAN SATISFY MULTIPLE CATEGORIES.



FOR ii=0,3 DO BEGIN
  if n[0,ii] ne -1 then type[n[where(n[*,ii] ne -1),ii],ii] = 'noise'
  if lang[0,ii] ne -1 then type[lang[where(lang[*,ii] ne -1),ii],ii] = 'langmuir'
  if iaw[0,ii] ne -1 then type[iaw[where(iaw[*,ii] ne -1),ii],ii] = 'ion acoustic'
  if wh[0,ii] ne -1 then type[wh[where(wh[*,ii] ne -1),ii],ii] = 'whistler broadband'
  if law[0,ii] ne -1 then type[law[where(law[*,ii] ne -1),ii],ii] = 'whistler narrowband'
ENDFOR


if dx[0] ne -1 then type[dx,0] = 'dust'
if dy[0] ne -1 then type[dy,1] = 'dust'
if dz[0] ne -1 then type[dz,2] = 'dust'
if dp[0] ne -1 then type[dp,3] = 'dust'

if dx2[0] ne -1 then type[dx2,0] = 'dust2'
if dy2[0] ne -1 then type[dy2,1] = 'dust2'
if dz2[0] ne -1 then type[dz2,2] = 'dust2'
if dp2[0] ne -1 then type[dp2,3] = 'dust2'


if mx1[0] ne -1 then type[mx1,0] = 'missing data'
if my1[0] ne -1 then type[my1,1] = 'missing data'
if mz1[0] ne -1 then type[mz1,2] = 'missing data'
if mp1[0] ne -1 then type[mp1,3] = 'missing data'

if thruster[0] ne -1 then type[thruster,0] = 'thruster/calibration'
if thruster[0] ne -1 then type[thruster,1] = 'thruster/calibration'
if thruster[0] ne -1 then type[thruster,2] = 'thruster/calibration'
if thruster[0] ne -1 then type[thruster,3] = 'thruster/calibration'



;-----------------------
;FIX BAD IDENTIFICATIONS
;-----------------------

;occasionally burst captures get mislabeled dust or noise when they are actually
;Langmuir or ion acoustic waves. These will always have a power peak at low freqs (<<1 kHz)
;Note that b/c of low frequency noise it's not always possible to find the frequency
;of the interesting wave

;This doesn't work well
;boo = where((type[*,0] eq 'dust') and (x[1].freq ge 1000.) and (x[1].freq lt 10000.))
;boo = where((type[*,0] eq 'dust') and (x[1].freq ge 10000.))


;bad1 = where((type[*,0] eq 'dust') and (x[1].bandwidth_stddev gt 10)); and (x[1].freq lt 10000.))
;bad2 = where((type[*,1] eq 'dust') and (x[2].bandwidth_stddev gt 10)); and (x[2].freq lt 10000.))
;bad3 = where((type[*,2] eq 'dust') and (x[3].bandwidth_stddev gt 10)); and (x[3].freq lt 10000.))
;bad4 = where((type[*,3] eq 'dust') and (x[4].bandwidth_stddev gt 10)); and (x[4].freq lt 10000.))

;if bad1[0] ne -1 then type[bad1,0] = 'ion acoustic or langmuir'
;if bad2[0] ne -1 then type[bad2,1] = 'ion acoustic or langmuir'
;if bad3[0] ne -1 then type[bad3,2] = 'ion acoustic or langmuir'
;if bad4[0] ne -1 then type[bad4,3] = 'ion acoustic or langmuir'

;bad1 = where((type[*,0] eq 'dust2') and (x[1].bandwidth_stddev gt 10)); and (x[1].freq lt 10000.))
;bad2 = where((type[*,1] eq 'dust2') and (x[2].bandwidth_stddev gt 10)); and (x[2].freq lt 10000.))
;bad3 = where((type[*,2] eq 'dust2') and (x[3].bandwidth_stddev gt 10)); and (x[3].freq lt 10000.))
;bad4 = where((type[*,3] eq 'dust2') and (x[4].bandwidth_stddev gt 10)); and (x[4].freq lt 10000.))

;if bad1[0] ne -1 then type[bad1,0] = 'ion acoustic or langmuir'
;if bad2[0] ne -1 then type[bad2,1] = 'ion acoustic or langmuir'
;if bad3[0] ne -1 then type[bad3,2] = 'ion acoustic or langmuir'
;if bad4[0] ne -1 then type[bad4,3] = 'ion acoustic or langmuir'


;This doesn't work very well.
;bad1 = where((type[*,0] eq 'noise') and (x[1].bandwidth_stddev gt 10)); and (x[1].freq lt 10000.))
;bad2 = where((type[*,1] eq 'noise') and (x[2].bandwidth_stddev gt 10)); and (x[2].freq lt 10000.))
;bad3 = where((type[*,2] eq 'noise') and (x[3].bandwidth_stddev gt 10)); and (x[3].freq lt 10000.))
;bad4 = where((type[*,3] eq 'noise') and (x[4].bandwidth_stddev gt 10)); and (x[4].freq lt 10000.))


;for i=0,200 do print,x[1].times[bad1[i]] + ' : ' + strtrim(x[1].freq[bad1[i]],2) + ' : ' + type[bad1[i]]

;for i=0,299 do print,x[1].times[i],type[i,0],strtrim(x[1].quality[i],2),strtrim(x[1].bandwidth[i],2),$
;	strtrim(x[1].bandwidth_stddev[i],2),strtrim(x[1].freq[i],2),format='(a23,4x,a23,4x,f10.3,4x,f10.3,4x,f10.3,4x,f10.3)'

;for i=0,500 do print,x[1].times[i],type[i,0],type[i,1],type[i,2],type[i,3],format='(a23,4x,a23,4x,a23,4x,a23,4x,a23)'


;openw,lun,'~/Desktop/cindy.txt',/get_lun
;for i=0,200 do printf,lun,x[1].times[i],type[i,0],type[i,1],type[i,2],type[i,3],format='(a23,4x,a23,4x,a23,4x,a23,4x,a23)'
;close,lun
;free_lun,lun



tt1 = x[0]
str_element,tt1,'type',type[*,0],/ADD_REPLACE
tt2 = x[1]
str_element,tt2,'type',type[*,1],/ADD_REPLACE
tt3 = x[2]
str_element,tt3,'type',type[*,2],/ADD_REPLACE
tt4 = x[3]
str_element,tt4,'type',type[*,3],/ADD_REPLACE


x2 = [tt1,tt2,tt3,tt4]

return,x2

end

