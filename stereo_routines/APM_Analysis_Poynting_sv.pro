 pro APM_Analysis_Poynting_sv
;
;	modified from Apm_analysis
;	copied from EandDensity2.pro on 30 June 2016
;	Note that this program uses 8 Hz data,  typical name =
;	magfilename='STB_L1_MAG_RTN_20130805_170130.TXT'  rather than
;		1 min MAGPLASMA data used in APM_Analysis

;	in contrast to APM_analysis_poynting_readin, this uses CDAWeb data
;   UCLA data is years first, and CDAWeb data is days first
;		As density data are not available at 8 Hz, density must be
;		entered by hand, but it doesnt matter much, except the program
;		will bomb if denst is 0.
;	APM data read in lines 346-            mag data lines 469 -
;	UR8 time of begin,end of event is Tseries[0], Tseries[2047]
;	UR8 time of first,last selected APM sample is
;
;
;	Temp_stuff.txt is used to store stuff for plotting by plotboth4
;	stuff_plot = fltarr(??) is plotted by plotboth4
;
;	There are 3 sets of momenta, used for averages:
;		VarA-VarF moments of full *phys sample set
;		VarAx-VarFx  moments from N1-N2
;		VarAxF-VarCxF  moments of selected bandpass filtered data
;
;	plotboth1		All_spect    comparison of spectra, APM, LRS
;	plotboth2 		ninestrip	 data from 3 APM's. 5 LRS. dens
;	plotboth2A		hodogram
;	plotboth3		correlations, XY, ZY, for APM and LRS    'LRS_Corr.ps'
;	plotboth4 		miscellany to be plotted
;	plotboth5		plot ninestrip_pub with x-y
;	plotboth6		plot 3 strips of ninestrip plus magnetic field
;	plotboth7
;	plotboth8		Waveforms,  E, B and poynting vector
;	plotboth9		APM_E,B_figure.ps
;	plotboth10
;	plotboth11		a different hodogram, for proposal
;
version = '13 Dec 2016'		; copied from APM_Analysis
version = '20 Feb 2017'		; changed sign  from potential to E in hodogrmams,
;							 Did this in first freq correction
;							  and changing back again in ninestrip and ninestrip_pub
version = '19 May 2017'     ; error in effective length corrected. A big factor
version = '22 May 2017'		; low f length set at 3. m, an approx
version = '26 May 2017'     ;  remove A,B,Cphys where calibraton is needed
;
;  now (10 sep 2015) it calls  PRO vmatrix, vector, N1, N2, eigenval, eigenvect
;  it loads LRSBurst 8 channels, lines 276-286.  In line 419 it call HK_density which
;	line 1033 or device Efield, calls antenna, which just provides the electrical
;		directions of the antennas, and antennaL, just length
;  2 Mar 2008 new normalization applied to all Fourier spectra
;	some fiddling of ranges, etc. to 1600
; 17 Feb new normalization of power, see NORMALIZATION below
;
;	Coordinate systems
;		ABCPhys[size]  is telemetered data corrected for frequency response
;				in line 515
;		ALL_APM[3,size] is ABCPhys combined
;		All_Phys_SC[3,size] is above transformed to Spacecraft system
;		All_Phys_RTN[3,size] is transformed to RTN system
;		All_Phys_B[3,size] is system with 0 along B,
;		APMV[3,size] is ABCPhys minus averages			;line 880
;		ABCWave[nfilt] are fft's of ABCPhys				;line 920
;			used for determining frequency of chosen sample
;			bandpass filter line 1070
;		ALL_Phys_BF[3,size] is bandpass filtered ALL_Phys_B
;			but I think this did not meet expectations for hodogram
;
;
;**************
; ihc = 1
; goto, plotboth4
;**************
;
@my_colors.inc
;
StreamId	= 1L
Err		= long(0)
Errors		= long(0)
Itype		= long(0)
fp		= fltarr(4)
junk		= 'junk'
nxt			= intarr(7)
;
size=2048
;
UR8		= dblarr(1)
UR8ATT		= -1.D03
UR8ATRD		= -2.D03
SCET		= 1.d00
scett		= 1.d00
select		= 'c'
pquestion	= 'q'
sc			= 'A'
yyyymmdd	= 20061101
hhmmss		= 120000
ant_char	= ['Ex','Ey','Ez','Ex','Ey','Ez','Ex','Ex']
ant_temp	= 'ex_apm'
Aseries 	= lonarr(size)
Bseries 	= lonarr(size)
Cseries 	= lonarr(size)
Dseries 	= lonarr(size)
Eseries 	= lonarr(size)
Fseries 	= lonarr(size)
Gseries 	= lonarr(size)
Hseries 	= lonarr(size)
Pseries		= lonarr(size)
Tseries		= fltarr(2050)
Btime		= fltarr(1030)
Bnow		= fltarr(4)
Bnowt		= fltarr(4)
Bnext		= fltarr(3)
nmin		= 0
sec1 		= 0
DOY			= 1
denst		= 8.22			; for 20130805,1701
denst		= 4.22			; for 20130805,2058
denst		= 2.57			; for 20130806,1008
vswt  		= 4.56
tideg		= 3.21
timeax 		= fltarr(2048)
rmscts		= fltarr(8)
rmsmv		= fltarr(8)
xsave		= fltarr(2)
ysave		= fltarr(2)
tot_avr		= fltarr(8)
infinity 	= -1
magfilename = 'pjk'
count 		= 0
eps0 = 8.8419e-12
 mu0 = 4.*!pi*1.E-7
;
print,'alive...'
;
print,'choose spacecraft, type A or a or B or b'
read,sc
if sc eq 'a' then sc = 'A'
if sc eq 'b' then sc = 'B'
scname = "Stereo_A"
if sc eq "B" then scname = "Stereo_B"
;
select = 'c'
;print,'type: r for realtime data or c for choose from archive'
;read,select
;
  Itype = 2
;  print,'type 0 to run through a file without stopping'
;  print,'type 1 to run through after a time to be entered'
;  print,'type 2 to do all events after a time to be entered with pause to print'
;  print,'at the pause, to come later, print "p" to print, else return'
;  read, Itype
;  print,Itype
;
;help,/source_files,output = temp
Iset = 0
if Itype gt 0 then begin
  print,'type start time as yyyymmdd,hhmmss'
  read, yyyymmdd, hhmmss
  print, 'start time ',yyyymmdd,hhmmss
  hhmmssst = hhmmss
  make_SCET,yyyymmdd,hhmmss,UR8Start
  print,'UR8 start time',UR8start
  Iset = 1
endif
;
Err = TM_Select_Domain(StreamId,"STEREO", scname, "SWaves", "LRSburst")
if err ne 0 then print, 'Domain - Error = ',Err
;
	  if err ne 0 then begin
		bsize = 1000L
		Err = TM_Error_Stack(StreamId,"name",errbuff,bsize,sizet)
		print,'Error after Find_event = ',"name" + errbuff
		Err = 	 TM_Error_Stack(StreamId,"description",errbuff,bsize,sizet)
		print,'Error description ' , errbuff
		Err = TM_Error_Stack(StreamId,"stacktrace",errbuff,bsize,sizet)
		print,'Error stacktrace ', errbuff
		stop
	  endif
;
    Err = TM_Select_Stream_Timerange(StreamId,UR8start,-1.d00)
    print,'timerange,streamid err = ',err,streamid
    err = TM_Get_Item_char(StreamId,"magFilename", timefile, 200L,sizet)
    print,'size,err,filetitle',sizet,err,timefile
    Err = TM_Get_Item_char(StreamId, "filepath", fpath, 200L, sizet)
    print,'size,err,filepath',sizet,err,fpath
    fixday = fix(UR8start)
    print,'fixday',fixday
    hrday0 = double(fixday)

print,'Connected to LRSburst Stream, STREAMID = ',Streamid
;stop
;
if Itype gt 0 then begin
  print,'going to set position at ',UR8start
  errSP=0
  if Itype gt 1 and Iset gt 0 then errSP = TM_Set_Position  (StreamId,UR8start)
  Iset = 0
  print,'set position err = ',errSP
endif
;
if select eq 'r' then filetitle = 'realtime'
;printf,3,filetitle
;printf,4,filetitle
;
;print,'calling systime'
J_time = systime(1,/Julian)
print,'J_time',J_time
CALDAT,J_time,month,day,year,hour,minute,second
timeS=string(hour,minute,second,year,month,day,format="(i2,':',i2.2,':',i2.2,' ',i4,'.',i2.2,'.',i2.2)")
print,'Getting the first event at: ',timeS

delt = 1./64.
deltms = 1000.*delt
fundfr = 1./(delt*size)
freq = fundfr*indgen(size/2)
print,'freq',freq[1],freq[size/2-1]
;
;fcomp = complexarr(size)
print,'fourcheck,size=',size
Afourx = complexarr(size)
Bfourx = complexarr(size)
Cfourx = complexarr(size)
Dfourx = complexarr(size)
Efourx = complexarr(size)
Ffourx = complexarr(size)
Gfourx = complexarr(size)
Hfourx = complexarr(size)
ERfour = complexarr(size)
ETfour = complexarr(size)
ENfour = complexarr(size)
DENfour =complexarr(size)
Aphys  = fltarr(size)
Bphys  = fltarr(size)
Cphys  = fltarr(size)
;AFphys = fltarr(size)
;BFphys = fltarr(size)
;CFphys = fltarr(size)
Dphys  = fltarr(size)
Ephys  = fltarr(size)
Fphys  = fltarr(size)
Gphys  = fltarr(size)
Hphys  = fltarr(size)
Pphys  = fltarr(size)
Denphys = fltarr(size)
Pfourx = complexarr(size)
ATT_Matrix  = dblarr(3,3)
ATT_MX_read = dblarr(3,3)
Ayyyy		= 1927
Adoy		= 11
ATTsec		= 86400
Aflag		= 0
apmV	= fltarr(3)
corner = fltarr(2)
xtemp	= fltarr(2)
ytemp	= fltarr(2)
pdens	= 0.
Binout	= 1
;
apmcorr = complexarr(size)
apmgain,0.,gapmc, gapmmag, gapmph
apmcorr(0) = conj(gapmc)
LRScorr = complexarr(size)
; write files for plotting gains
; for the fft, n = 1 and size-1 are complex conjugates
;openw, 6, 'lrsgain.dat'
;openw, 7, 'apmgain.dat'
openw,8, 'temp_stuff.txt'
; openw,66,'EandDesntiy2_Stuff.txt'
	for n = 1,size/2 do begin
	  f = n*fundfr
	  apmgain,f,gapmc, gapmmag, gapmph
	  apmcorr(n) = conj(gapmc)
	  apmcorr(size-n) = gapmc
	  PALRSgain,f,glrsc, glrsmag, glrsph
;	  LRScorr(n) = conj(glrsc)
;	  LRScorr(size-n) = glrsc
	  LRScorr(size-n) = conj(glrsc)
	  LRScorr(n) = glrsc
;	  print,f,gapmc,gapmmag,gapmph
;	  printf,7,f,gapmc,gapmmag,gapmph
;	  printf,6,f,gLRSc,glrsmag,glrsph
	endfor
LRScorr(0) = complex(1., 0.)
; **************

nremove=63		; remove the lowest nremove frequencies
LRScorr(1:nremove) = complex(1., 0.)
lowrem = size-nremove
LRScorr(lowrem:size-1) = complex(1., 0.)
;
timesc = -1.
;
;
;	open UCLA OR CDAWEB magnetic data file and read header
;		now it should be 8 or 22 samples per sec, not magplasma
;
;magfilename='./mag_data/ST'+sc+'_'+string(yyyymmdd,format='(I8)')+'.dat'
;magfilename='./mag_data/ST'+sc+'_'+string(yyyymmdd)+'.dat'
;	for the files from CDAWeb
;magfilename='STA_L2_MAGPLASMA_1M_CARR2054.TXT'
;	starts 20070304
;magfilename='STA_L2_MAGPLASMA_1M_CARR2099TRY.TXT'
;magfilename='STA_L2_MAGPLASMA_1M_CARR2099.TXT'
;ATTfilenam ='Att_pointing_Ahead_2010_RTN.txt'
;
;magfilename='STB_L2_MAGPLASMA_1M_CARR2140.TXT'
;magfilename='STB_L1_MAG_RTN_20130805_133100.TXT'
;magfilename='STB_L1_MAG_RTN_20130805_170130.TXT'
;magfilename='STB_L1_MAG_RTN_20130805_1728.TXT'
;magfilename='STB_L1_MAG_RTN_20130805_205800.TXT'
;magfilename='STB_L1_MAG_RTN_20130807_094800.TXT'
;magfilename='STB_L1_MAG_RTN_20130820.TXT'
;magfilename='STB_L1_MAG_RTN_20130806_100800_CDAW.TXT
;magfilename='STB_L1_MAG_RTN_20130805.TXT'		; year first
;magfilename='STB_L1_MAG_RTN_20130805_100000_CDAW.TXT'
;magfilename='STB_L1_MAG_RTN_20130827.TXT'
magfilename='STB_L1_MAG_RTN_20130827_160000.TXT'
; starts 20130804
ATTfilename='ATT_Pointing_Behind_2013_RTN.txt'
;resultfilename='APM_Analysis_results_B_Carr2140.txt'
;result2filename='APM_Analysis_results2_B_Carr2140.txt'
resultfilename='APM_Analysis_results_trial.txt'
result2filename='APM_Analysis_results2_trial.txt'
;resultfilename='APM_Analysis_results_B_Carr2140_2Hz.txt'
;result2filename='APM_Analysis_results2_B_Carr2140_2Hz.txt
;
;magfilename='STB_L2_MAGPLASMA_1M_CARR2140X.TXT'
print,'mag file name ',magfilename
openr,9,magfilename,error=magerror
print,' at open magerror and count ',magerror
if (magerror ne 0)  then goto,getmag
;for i = 0,51 do begin  then 45, then 55
; this uses CDAWeb data,  use 47
;for UCLA use 71 or 72  last line sghould read "data:"
if magerror eq 0 then begin
  for i = 0,47 do begin
    readf,9,junk
    print,junk
  endfor
endif
print,'last magfile line  ',junk
;
;	open ATTitude file and advance
 openr,22,ATTfilename,error=ATTerrpr
 print,'UR8start ',UR8start
 while UR8ATRD lt UR8Start do begin
    ATT_Matrix[*,*] = ATT_MX_Read[*,*]
    UR8ATT = UR8ATRD
    readf,22,Ayyyy,ADoy,Attsec,Aflag,ATT_MX_Read
	OK = TM_UR8_from_Ydoy(UR8ATRD,Ayyyy,Adoy,1000.*ATTsec)
;    print,'ATT ',OK,UR8Att,Att_MX_Read[0,*]
 endwhile
     print,'ATT ',OK,UR8Att,Att_Matrix[0,*]
     print,'ATT ',OK,UR8AtRD,Att_MX_Read[0,*]
;
;while count gt infinity do begin
newevent: count = count + 1

;	err = TM_Find_Event(Streamid)
;	if err ne 0 then print, 'EVENT - Error = ',Err

	print,'LRSburst event number: ',count,format="(a,i6)"

	J_time = systime(1,/Julian)
	CALDAT,J_time,month,day,year,hour,minute,second
	timeS=string(year,month,day,hour,minute,second,$
	format="(i4,'/',i2.2,'/',i2.2,' ',i2,':',i2.2,':',i2.2)")
	print,' acquired at: ',timeS
	timeP = timeS

;GET LRS AND APM DATA

	Errors = Err
	last = 'no'
	in = 0
	for in=0,7 do begin
;
	  err = TM_Find_Event(Streamid)
   	  if err ne 0 then print, 'EVENT - Error = ',Err
	  Err = TM_Get_Item_char(StreamId, "spacecraft_long", scname, 200L, sizet)
	  Err = TM_Get_Item_char(StreamId, "spacecraft_short", sc, 2L, sizet)
	  print,'sc name,err,size= ' ,sc,err,sizet
;
	  if err ne 0 then begin
		Err = TM_Error_Stack(StreamId,"name",errbuff,bsize,sizet)
		print,'Error after Find_event = ',"name" + errbuff
		Err = 	 TM_Error_Stack(StreamId,"description",errbuff,bsize,sizet)
		print,'Error description ' , errbuff
		Err = TM_Error_Stack(StreamId,"stacktrace",errbuff,bsize,sizet)
		print,'Error stacktrace ', errbuff
	  endif
;
          Err = TM_Get_Item_I4(StreamId, "channel", Ich, 1L, size1)
    	  Err = TM_Get_Item_R4(StreamId,'Temperature_Converter', TCvt, 1L, sizet)
;*********  a special section to search for events when B is parallel to X or Z antennas
;	while not(EOF() do begin
;	    	err = TM_Find_Event(Streamid)
;   	  		if err ne 0 then print, 'EVENT - Error = ',Err
;	  while SCETMAG LT scet do begin
;	    readf,9, daym,mon,nyear,hr,mmin,sec1,usec,denst,vswt,tideg,$
;	    Bnowt,format='(I2,x,I2,x,I4,I3,x,I2,x,I2,x,I3,6F14.5)'
;	    Bnow[0:2] = Bnowt[0:2]
;	    Bnow[3] = sqrt(Bnow[0]^2 + Bnow[1]^2 + Bnow[2]^2)
;	    err = TM_UR8_from_YMD(scetmag,nyear,mon,daym,hr,mmin,sec1,usec)
;
;	  ENDwhile
;	endwhile
;*********  end of special section
;
      if ich eq 0 then begin
	    Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Aseries, 2048L, sized)
	    for nt = 0,2047 do begin
	    	if abs(Aseries[nt]) gt 32000 then $
	    		Aseries[nt] = (Aseries[nt-1]+Aseries[nt+1])/2
	    endfor
	  endif
;
; 0,A is APM X, 1,B is APM Y, 2,C is APM Z, 3,D is APM X-Y
; 4,E is LRS X, 5,F is LRS Y, 6,G is LRS Z, 7,H is LRS X_Y
;
if ich eq 1 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Bseries, 2048L, sized)
if ich eq 2 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Cseries, 2048L, sized)
if ich eq 3 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Dseries, 2048L, sized)
if ich eq 4 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Eseries, 2048L, sized)
if ich eq 5 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Fseries, 2048L, sized)
if ich eq 6 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Gseries, 2048L, sized)
if ich eq 7 then Err = TM_Get_Item_I4(StreamId, "Time_Series_Raw", Hseries, 2048L, sized)

;	    print,'after raw, size=',sized
;	    print,'After time_series_raw, error is:',Err,format="(a,i)"
;	    if err ne 0 then stop
	    Err = TM_Get_Item_char(StreamId, 'antenna',ant_temp,10L,sizech)
;	    print,'in,channel',in,ich,sizech,ant_temp
	    ant_char[ich] = ant_temp
	    Err = TM_Get_Item_char(StreamId, "Event_is_Last", last, 10L, sizech)
;	    print,'last test ',last
;	    print,' the returned size,err was:',sizech,format="(a,i)"
;	    if last eq 'yes' then goto, process
	    if ich eq 7 then goto, process
	endfor
;
	process: in = 0
;	Err = TM_Get_Item_R8(StreamId, "CCSDS_SCET_UR8", SCET, 1, tsize)
	Err = TM_Get_Item_R8(StreamId, "Time_series_SCET_UR8",Tseries, 2048, tsize)
	print,'series time check',Tseries[0],Tseries[tsize-1],Tseries[tsize-1]-tseries[0]
	scet = Tseries[0]
	Err = TM_UR8_to_string_comparison( SCET, timesc)
	print,'at process, sc event time ',scet,' ',timesc
	start = scet - double(2047.)/64./864.d02
;REMOVE 1553 GLITCHES IF IFIX EQ 1

	ifix = 0
	if ifix eq 1 then begin
;	  find an unglitched sample
	  n = 0
	  while Aseries(n) do n = n+1
	  nun = n
	  print, 'nun,A',nun,Aseries(nun)
	  Aseries(0) = Aseries(n)
;
	  for n = 1,size-1 do begin
	    if Aseries(n) then Aseries(n) = Aseries(n-1)
	  endfor
	endif
;
; For full APM sample
;		convert to SC coordinates
;		convert to RTN coordinates
;		convert to B-V coordinates
; For selected sample only,  in Plotboth2A
;		bandpass filter
;		Transform to B system
;		plot
;		save results
;
	timeax = findgen(2048)*32./2048.			; timeax is in sec
;

SCETMAG = -1.d00
;stop
;
;stop
;	*************
;	temporary
;
;	GET MAGNETIC FIELD FROM UCLA-IGPP FILES IN /MAG_DATA
;
	scetmag = -1.d00
;
 NBtime = 32
;
 if magerror eq 0 then begin
;
;	for the files from CDAWeb
	print,' mag,APM time ',scetmag,scet
	Scetmag = 0.D00
	scet = Tseries[0]
;   hunt for sample time matching APM event start
	  readf,9, daym,mon,nyear,hr,mmin,sec1,usec,Bnowt,$
	  format='(I2,x,I2,x,I4,I3,x,I2,x,I2,x,I3,4F24.5)'
	while SCETMAG LT scet do begin
	  print,' '
	  print,'going to read 9 to hunt, scetmag,scet ',scetmag,scet
	  Bnow[0:2] = Bnowt[0:2]
	  print,'Bnowt ',Bnowt
	  Bnow[3] = sqrt(Bnow[0]^2 + Bnow[1]^2 + Bnow[2]^2)
;		probably for magdata
;	  readf,9, daym,mon,nyear,hr,mmin,sec1,usec,denst,vswt,tideg,$
;	  Bnowt,format='(I2,x,I2,x,I4,I3,x,I2,x,I2,x,I3,6F14.5)'
	  readf,9, daym,mon,nyear,hr,mmin,sec1,usec,Bnowt,$
	  format='(I2,x,I2,x,I4,I3,x,I2,x,I2,x,I3,4F24.5)'
	  err = TM_UR8_from_YMD(scetmag,nyear,mon,daym,hr,mmin,sec1,usec)
	  print,'tchk ',hr,mmin,sec1,usec
	  print,'mag chk ',scetmag-scet,Bnowt
;	  if denst gt 0. then stop
	ENDwhile
	print,'scet,mag,diff ',scet,scetmag,scetmag-scet,format='(A16,3F16.8)'
;	load B data, assumed to be 8 samples per sec so 256 samples
;	?? the file of B data is usually larger than an event, 32 sec
;		But the attempt here seems to be to load only the event data
;		32 sec should be 256 B samples
	sizeB = 257
	Btime = dblarr(sizeB)
	Bevent = fltarr(4,SizeB)
	Btime[0] = scetmag					; first scetmag later than Tseries[0]
	print, 'scetmag ',scetmag,format='(A20,F14.8)'
	print,'sizeB ',sizeB
	Btimepl = fltarr(sizeB)
	Bevent[0:2,0] = Bnowt[0:2]
	Btimepl[0] = 86400.*(Btime[0]-Tseries[0])
;	calculate number of magnetic field samples
    No_Bsamp = sizeB-1
	for nload = 1,No_Bsamp-1 do begin
	  readf,9, daym,mon,nyear,hr,mmin,sec1,usec,Bnowt,$
	  format='(I2,x,I2,x,I4,I3,x,I2,x,I2,x,I3,4F24.5)'
	  err = TM_UR8_from_YMD(scetmag,nyear,mon,daym,hr,mmin,sec1,usec)
	   print,'tchk2 ',hr,mmin,sec1,usec
	  Btime[nload] = scetmag
	  Bevent[0:2,nload] = Bnowt[0:2]
	  Btimepl[Nload] = 86400.d00*(SCETmag-Tseries[0])
	  print,'load ',nload,BtimePL[Nload]
	endfor
	NBload = Nload-1
	print,'nbload ',nbload
;	now B samples corresponding to the whole selection have been loaded
;		there are Nload of them, which should be about (N2-N1)/8
;
    print,'Btime 1st last',Btime[0],Btime[nload-1],No_Bsamp-1
	Bnorm = fltarr(3)
	print,'chk denst  ',denst
	if denst gt 0. then begin
 ; 	  print,'bnow check',scetmag,tseries[0],bnow
	  Bnorm[0:2] = Bnow[0:2]/Bnow[3]
;	 	in MagPlasma, B is in RTN, Binout = 0 is B away from sun
;		0 is out, 1 is in
	  Binout = 0
;     R direction is out, so Binout = 0 for out
      if (Bnow[0] - Bnow[1]) lt 0. then Binout = 1
      Binoutsv = binout
;
	  NB = NBtime
	  Bmax = max(Bevent[0,0:NB-1],min=Bmin0)>$
	  max(Bevent[1,0:NB-1],min=Bmin1) >$
	  max(Bevent[2,0:NB-1],min=Bmin2)
	  Bmin = Bmin0<Bmin1<Bmin2
	endif
 endif
;
 getmag:print,'magerror ',magerror		; entry for no magnetic data
;
;CORRECT FOR FREQUENCY AND PLOT THE FREQUENCY CORRECTED DATA ON THE SCREEN, and B

;*******
;Aarons additions
white = 0
black = 1
;*******

	!p.background=white
	!p.color=black
	!p.multi=[0,1,9]
	!x.range=[timeax[0],timeax[2047]]
	!y.range=[0,0]
	nulltick= replicate('  ',30)
;
	print,'freq check',size,freq[1],freq[size/2-1]
;
	y0 = .86
	y1 = .95
	dy = .1
;
;	Afourx = (2.5/6.5536e4)*fft(float(Aseries[0:size-1]))/APMcorr
	Afourx = -(2.5/6.5536e4)*fft(float(Aseries))/APMcorr
	Aphys = float(fft(Afourx,/INVERSE))
	Bfourx = -(2.5/6.5536e4)*fft(Bseries[0:size-1])/APMcorr
	Bphys = float(fft(Bfourx,/INVERSE))
	Cfourx = -(2.5/6.5536e4)*fft(Cseries[0:size-1])/APMcorr
	Cphys = float(fft(Cfourx,/INVERSE))
	Dfourx = -(2.5/6.5536e4)*fft(Dseries[0:size-1])/LRScorr
	Dphys = float(fft(Dfourx,/INVERSE))
	Efourx = -(2.5/6.5536e4)*fft(Eseries[0:size-1])/LRScorr
	Ephys = float(fft(Efourx,/INVERSE))
	Ffourx = -(2.5/6.5536e4)*fft(Fseries[0:size-1])/LRScorr
	Fphys = float(fft(Ffourx,/INVERSE))
	Gfourx = -(2.5/6.5536e4)*fft(Gseries[0:size-1])/LRScorr
	Gphys = float(fft(Gfourx,/INVERSE))
	Hfourx = -(2.5/6.5536e4)*fft(Hseries[0:size-1])/LRScorr
	Hphys = float(fft(Hfourx,/INVERSE))
;	print,'var check',aphys[0],aphys[1],aphys[size-2],aphys[size-1]
;	print,'four check',afourx[1],afourx[size-1]
	varA = moment(Aphys(0:size-1), sdev=Arms)
	varB = moment(Bphys(0:size-1), sdev=Brms)
	varC = moment(Cphys(0:size-1), sdev=Crms)
	varD = moment(Dphys(0:size-1), sdev=Drms)
	varE = moment(Ephys(0:size-1), sdev=Erms)
	varF = moment(Fphys(0:size-1), sdev=Frms)
;	print,'A var',varA,Arms,format='(A8,5f10.3)'
;	print,'A var',varA,Arms
;	print,'B var',varB,Brms
;	print,'C var',varC,Crms
;	print,'D var',varD,Drms
;	print,'E var',varE,Erms
;	print,'F var',varF,Frms
	tot_avr[0] = VarA[0]
	tot_avr[1] = VarB[0]
	tot_avr[2] = VarC[0]
	apmV1 = [varA[0], varB[0], varC[0]]
	All_APM			= fltarr(3,size)
;	print,'all_apm size = ',size
	All_APM[0,*] = Aphys[0:size-1]
	All_APM[1,*] = Bphys[0:size-1]
	All_APM[2,*] = Cphys[0:size-1]
;
;********  temporary
  	Temp = 23.7
	if sc eq 'B' then Temp = 34.5
;
nepoch = 3
if scet gt 9075.d00 then nepoch = 4
if scet gt 9086.d00 then nepoch = 5
if scet gt 9097.d00 then nepoch = 6
if scet gt 9112.d00 then nepoch = 7
if scet gt 9115.d00 then nepoch = 8
;
B = [-6.8, 8.7, -.2]
;getmag,'B',scet,32.,Blist
;
;	  high pass filter or band pass filter
;
 Ifilt = 0
 if Ifilt ne 0 then begin
   order = 4
   hpfilter = Digital_Filter(1./64., 1., 50., order)
;	example   AFphys = Convol(Aphys,hpfilter)
;	print,'hpfilter ',hpfilter
; the IDL routine does not normalize????.  So do it
   Sumbp = 0.
   for nhp = 0,2*order do begin
     sumhp = sumhp + hpfilter[nhp]
   endfor
; print,'sumbp ',sumbp
   hpfilter[*] = hpfilter[*]/sumhp
;
   AFphys[0:size-1] = Convol(Aphys[0:size-1],hpfilter)
   BFphys[0:size-1] = Convol(Bphys[0:size-1],hpfilter)
   CFphys[0:size-1] = Convol(Cphys[0:size-1],hpfilter)
 endif
;
;	determine ratio of X to Y for APM
;	minimize sum(Aphys - alpha*Bphys)^2 with averages subtracted
;
	Sumsqxy = 0.
	Sumsqyy = 0.
	for n = 0,size-1 do begin
	  Sumsqxy = Sumsqxy + (Aphys[n]-varA[0])*(Bphys[n]-varB[0])
	  Sumsqyy = sumsqyy + (Bphys[n]-varB[0])^2
	endfor
	alpha = Sumsqxy/Sumsqyy
;	print,'averages,x,y',varA[0],varB[0]
;	print,'ratio X to Y',alpha
	xyratioAPM = alpha
;
;	determine ratio of Z to Y for APM
;	minimize sum(Cphys - alpha*Bphys)^2 with averages subtracted
	Sumsqzy = 0.
	Sumsqyy = 0.
	for n = 0,size-1 do begin
	  Sumsqzy = Sumsqzy + (Cphys[n]-varC[0])*(Bphys[n]-varB[0])
	  Sumsqyy = sumsqyy + (Bphys[n]-varB[0])^2
	endfor
	alpha = Sumsqzy/Sumsqyy
;	print,'averages,z,y',varC[0],varB[0]
;	print,'ratio APM Z to Y',alpha
	zyratioAPM = alpha
;
;  XformSC converts volt and calculates matrix to transform from Ant to SC system
;		but does not make the transform
;	call XformSC to calculate the inverse matrix for the
;	coordinate transformation from antenna voltages to E in the
;	spacecraft system.  XformSC calls antenna_length
;	and antenna_vector to calculate the matrix
 ctr_freq	= 1.				; FAKE
 print,'calling XformSC ',sc,ctr_freq,denst
 XformSC, sc, ctr_freq, denst, invmatrix
; print,'inverse '
; print,invmatrix
 All_Phys_SC = fltarr(3,size)
;
 for nxf = 0,size-1 do begin
 	All_Phys_SC[*,nxf] = InvMatrix#All_APM[*,nxf]
 endfor
;
 ant_check = 0
;	a section to check orientation of antennas in RTN and Bsys
 if ant_check ne 0 then begin
   antenna_vector,sc,antV
;	in spacecraft coordinates
;  print,' in XformSC, AntV ='
;  print,'Z',antv[*,0]
;  print,'Y',antv[*,1]
;  print,'X',antv[*,2]

; TRANSFORM TO RTN SYSTEM
;
  Ant_RTN = fltarr(3,3)
  for nxf = 0,2 do begin
     Ant_RTN[*,nxf] = Att_Matrix#AntV[*,nxf]
     print,'Ant_Rtn ',nxf,Ant_RTN[*,nxf]
   endfor
;  stop
 endif
;
;  CHECK if new RTN matrix is needed
;		not done yet
;
; TRANSFORM TO RTN SYSTEM
;
  All_Phys_RTN = fltarr(3,size)
  for nxf = 0,size-1 do begin
     All_Phys_RTN[*,nxf] = Att_Matrix#All_Phys_SC[*,nxf]
;    print,All_APM[*,nxf]
;    print,All_Phys_RTN[*,nxf]
  endfor
;
;	transform to B-V system, first axis along B and outward from sun
;	second axis perp to B in plane of B and Vsw, taken as radial
;		so parallel to
;	third axis perp to 1 and 2
;		 VSW taken as radial
;
 CBmatrix = fltarr(3,3)
 CBmatrix[*,2] = Bnorm[0:2]
 if binout eq 1 then CBmatrix[*,2] = -Bnorm[0:2]
;
 ALL_Phys_B = fltarr(3,size)
 Vtemp1 = fltarr(3)
 Vtemp2 = fltarr(3)
 cbnorm = sqrt(Bnow[1]^2 + Bnow[2]^2)
 CBmatrix[*,0] = [0.,-Bnow[2], Bnow[1]]/cbnorm
 print,'Bnow,cbnorm ',Bnow,cbnorm
; dot product of 2 with 0 is zero, as is dot product with R in RTN system
 Vtemp2 = CBmatrix[*,0]
 Vtemp1 = CrossP(Bnorm,Vtemp2)
 print,'vtemp2 ',vtemp2
 print,'Bnorm  ',bnorm
 print,'vtemp1 ',vtemp1
;
 CBmatrix[*,1] = Vtemp1[*]/norm(Vtemp1)
print,'B',Bnow[*]
print,'cBnorm ',cBnorm
print,'Binout ',binout
print,'VM[0] ',CBmatrix[*,0]
print,'VM[1] ',CBmatrix[*,1]
print,'VM[2] ',CBmatrix[*,2]
;
;	check for right handed coordinate system
;
Vtemp3 = fltarr(3)
vtemp3 = crossp(CBmatrix[*,0],CBMatrix[*,1])
;print,'Vtemp3 ',vtemp3
;print,'CBmatrix ',CBmatrix[*,2]
RHchk = Vtemp3[0]*CBmatrix[0,2]+Vtemp3[1]*CBmatrix[1,2]+Vtemp3[2]*CBmatrix[2,2]
if RHchk lt 0. then begin
  print,'left handed system, RHchk',RHchk
  print,'B',Bnow[*]
  print,'cBnorm ',cBnorm
  print,'Binout ',Binout
;	reverse axis 0
  CBmatrix[*,0] = -CBmatrix[*,0]
  print,'new axes '
  print,'chmatrix[*,0] ',CBmatrix[*,0]
  print,'chmatrix[*,1] ',CBmatrix[*,1]
  print,'chmatrix[*,2] ',CBmatrix[*,2]
endif
;if RHchk lt 0. then stop
All_Phys_temp = fltarr(3)
;print,'binout',binout
;
;	check transform from RTN to B-V system,  That is, apply the transformation
;	to B.  Then B should be all 3rd component, Bz
;
;    cxform,Bnow[0:2],ALL_Phys_Temp,CBMatrix
;    print,'chek B RTN ,All_Phys_temp, it gave: ',All_Phys_temp
;   chek B RTN       1.19807 -7.91624e-09 -7.45058e-09
;
;  PRO cxform, oldvec, newvec, xmatrix
;      cxform,ALL_Phys_RTN[0:2,0],ALL_Phys_Temp[0:2],CBMatrix
;
 for nxf = 0,size-1 do begin
   cxform,All_Phys_RTN[0:2,nxf],ALL_Phys_Temp,CBMatrix
   ALL_Phys_B[0:2,nxf] = All_Phys_Temp[0:2]
;   print,'temp ',all_phys_temp
 endfor
 print,'CBmx ',CBmatrix
 print,'RTN ',All_Phys_RTN[*,1024]
 print,' B  ',ALL_PHYS_B[*,1024]
  print,' B  ',ALL_PHYS_B[*,1025]
 print,' B  ',ALL_PHYS_B[*,1124]
;
;
ihc = 1			; here 1 is screen
;
plotboth1: print,'ihc = ',ihc
	!p.background=white
	!p.color=black
	!p.multi=[0,1,3]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_Check.ps',xsize= 18.,ysize=20.,yoffset=4.,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
n1B = 50
n2B = 500
	plot,timeax[N1b:N2B],Bphys[N1B:N2B]-VarB[0],charsize=2.,$
	  Ytitle='ABCphys ',xstyle=17,xrange=[0.,0.],$
	  xticklen=.1,/nodata,/ynozero,/noerase
	oplot,timeax[ N1B:N2B],Aphys[ N1B:N2B]-VarA[0],color=red
	oplot,timeax[ N1B:N2B],Bphys[ N1B:N2B]-VarB[0],color=blue
	oplot,timeax[ N1B:N2B],Cphys[ N1B:N2B]-VarC[0],color=green
;
	plot,timeax[ N1B:N2B],Bphys[ N1B:N2B]-VarB[0],charsize=2.,$
	  Ytitle='ABCphys ',xstyle=17,$
	  xticklen=.1,/nodata,/ynozero,/noerase
	plot,timeax[ N1B:N2B],Bphys[ N1B:N2B]-VarB[0],charsize=2.,$
	  Ytitle='ABCphys ',xstyle=17,$
	  xticklen=.1,/nodata,/ynozero,/noerase
;
	timeS = systime(0)
	xyouts, 0.,.0,scname+' APM_Check.ps, version '+version+ ',  output at '+ timeS,$
    charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth1
;
    N1 = 1
	N2 = size-1
;
ihc = 0			; here 1 is screen but do hardcopy first
;
plotboth2: print,'at plotboth2, ihc = ',ihc
;
;	now 9 strip screen plot if B is present, 8 strips otherwise
;
	ninestrip:  print,'ninestrip'
	!p.background=white
	!p.color=black
	!p.multi=[0,1,9]
	y0 = .86
	y1 = .95
	dy = .1
;
	if ihc eq 1 then erase
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_ninestrip.ps',xsize= 25.,ysize=13.5,/landscape,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;	plot differences, first get data
	alrsxy = 1.05
	alrszy = 1.
	ratio_data_name = 'APM_antenna_ratios_B_Carr2140.txt'
	openr,37,Ratio_data_name
	readf,37,n_entry_no
	print,'entry no ',n_entry_no
	APM_density = fltarr(n_entry_no)
	xyAPMratio = fltarr(n_entry_no)
	ZyAPMratio = fltarr(n_entry_no)
	  readf,37,APM_Density
;	  print,APM_density
	  readf,37,xyAPMratio
	  readf,37,zyAPMratio
	close,37
	xyratio = interpol(xyAPMratio,APM_density,denst)
	zyratio = interpol(zyAPMratio,APM_density,denst)
;	print,'denst ',denst
	print,'denst,xyratio ',denst,xyratio
;	print,'density',APM_density
;	print,'zy ',zyAPMratio
;
	for np = 0,7 do begin
	  plottitle = ''
	  if np eq 0 then begin
		Pphys = Aphys[0:2047]
		plottitle = scname+' '+timeSC
		plcolor = red
;		print,'at 570 ihc,np,Aphys[1] ',ihc,np,Aphys[1]
  endif
;	  if ihc eq 1 then stop
	  if np eq 0 or np eq 3 then plcolor = red
	  if np eq 1 then Pphys = Bphys[0:2047]
	  if np eq 1 then Pseries = Bseries[0:2047]
	  if np eq 1 or np eq 4 then plcolor = blue
	  if np eq 2 then Pphys = Cphys[0:2047]
	  if np eq 2 then Pseries = Cseries[0:2047]
	  if np eq 2 or np eq 5 then plcolor = green
	  if np eq 3 then Pphys = Dphys[0:2047]
	  if np eq 3 then Pseries = Dseries[0:2047]
	  if np eq 4 then Pphys = Ephys[0:2047]
	  if np eq 4 then Pseries = Eseries[0:2047]
	  if np eq 5 then Pphys = Fphys[0:2047]
	  if np eq 5 then Pseries = Fseries[0:2047]
;	  plot differencesa
	  if np eq 6 then Pphys = Bphys[0:2047]-Aphys[0:2047]/xyratio
	  if np eq 7 then Pphys = Bphys[0:2047]-Cphys[0:2047]/zyratio
;	  print,'chk timeax, np ',np,timeax[0],timeax[2047]
;	  print,'chk oplot,timeax,pphys ',timeax[0],timeax[2047],$
;	  	Pphys[0],Pphys[2047]
;	  if np eq 0 then plottitle =  scname+' '+timeSC
	  if np eq 0 then begin
	    print,'starting np = 0 check '
	    print,'ihc ',ihc
		print,'title ',plottitle
		print,'pphys ',Pphys[0],Pphys[2047]
	    print,'time ',timesc
	  endif
;
;	print,'calling plot ihc,np = ',ihc,np
;	if (ihc eq 1) and (np eq 1) then stop      ; stops with Ex strip OK
;
;*******   Here change from E to potential
	ymax = max(-Pphys[0:size-1],min = ymin)<1.5
	Ymin = ymin>(-1.)
	if np ne 7 then plot,timeax,-Pphys,charsize=2.,Ytitle=ant_char[np],$
	  title=plottitle,xstyle=17,position=[.1, y0, .99, y1],yrange=[ymin,ymax],$
	  xticklen=.1,xtickname=nulltick,thick=2,/nodata,/ynozero,/noerase
;
	if np eq 7 then plot,timeax[0:2047],-Pphys[0:2047],charsize=2.,$
	  Ytitle=ant_char[np],xstyle=17,position=[.1, y0, .99, y1],$
	  xticklen=.1,thick=2,/nodata,/ynozero,/noerase
;
;	if (ihc eq 1) and (np eq 1) then stop 		; after 16 sec, only Ey strip shows
;	print,'return from plot,calling oplot, ihc,np = ',ihc,np
	  oplot,timeax,-Pphys[0:size-1],color=plcolor
;	  	  if (ihc eq 1) and (np eq 1) then stop
;	print,'return from oplot ihc,np = ',ihc,np
	  var = moment(float(Pseries[0:size-1]), sdev=ctrms)
;	  rmscts[0] = ctrms
;	print,'rms in counts ',ctrms
	  var = moment(Pphys[0:size-1], sdev=vrms)
;	  rmsmV[0] = 1000.*vrms
;	  print,'screen,avr in V, rms in mV',var, 1000.*vrms
	  prms = string(format='(f8.1,a4,f7.2,a3)',ctrms,' cts',$
	1000.*vrms,' mV')
	  xyouts, .70, .2*y0+.8*y1, prms,/normal
	  y0 = y0-dy
	  y1 = y1-dy
	endfor
;

;	  xticklen=.05,xtickname=nulltick,/ynozero,color=magenta
;	  xticklen=.05,xtickname=nulltick,/ynozero,color=cyan
;
 ;
	timeS = systime(0)
	xyouts,0.,0.,scname+'APM_ninestrip.ps, from APM_Analysis, version '+version+$
 ',  output at '+ timeS,charsize=.7,/normal
; 	      stop	OK
;

;    N1 = 1
;	N2 = size-1
    if ihc eq 1 then begin
	    print,'   type return to plot next event'
        print,'or TYPE C OR c TO RECORD THE POSITION OF THE CURSOR and HISTOGRAM'
        print,'or Type p to make ninestrip_pub or type Q or q to stop'
 		dispose = 'r'
		xsave2 = 0.
		ysave2 = 0.
    	READ, DISPOSE
    	IF (DISPOSE EQ 'Q') OR (DISPOSE EQ 'q') then stop
    	if dispose eq 'c' then dispose = 'C'
    	if dispose eq 'p' then goto, plotboth5
    	if dispose ne 'C' then goto, newevent
;
    	IF (DISPOSE EQ 'C') OR (DISPOSE EQ 'c') THEN begin
          YMAX = 0.
;
          print,'set the cursor in the upper plot for the range to be'$
          +'fitted, lower end first, and click'
;
          delt = 32./2048.
          print,'delt ',delt
 select:  FOR nch = 0,1 DO BEGIN
            CURSOR,XSAVET,YSAVET,/data
            print,'ihc,nxh,cursor ',ihc,nch,xsaveT,ysaveT
;			printf,?,xsave2,ysave2
			xsave[nch] = XsaveT
			ysave[nch] = YsaveT
;
            NSAMPLE = fix(XSAVET/DELT + .5)
            if nch eq 0 then N1 = NSAMPLE
            if nch eq 0 then tsamp1 = XsaveT
            if nch eq 1 then N2 = NSAMPLE
            if nch eq 1 then tsamp2 = XsaveT
            N2 = N2<size-1
            wait, .5
          ENDFOR
        	print,N1,N2,tsamp1,tsamp2
        	print,'delt ',delt
          						; end of cursor
          if (N2-N1) lt 5 then begin
			print,'bad selection,, set cursor again'
			goto, select
          endif
;          print,'now fit ',xsave[0],xsave[1],N1,N2
;
        endif 				; end "If DISPOSE --"
;
	endif					; If IHC eq 1
;
      n1f = n1
      n2f = n2
;   print,'at VarAx, N1,N2 = ',N1,N2
;	this will be changed to B system
 	varAx = moment(Aphys[N1:N2])
 	varBx = moment(Bphys[N1:N2])
 	varCx = moment(Cphys[N1:N2])
  	varDx = moment(Dphys[N1:N2])
   	varEx = moment(Ephys[N1:N2])
   	varFx = moment(Fphys[N1:N2])
;
if ihc ne 1 then begin
	device,/close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 0 then print,'going to plotboth2 from line 814 '
if ihc ne 0 then goto,plotboth2
;
; stop  got to here OK
;
ihcA = 1			; here 1 is screen
;
plotboth2A: print,'at plotboth2A, ihcA = ',ihcA
	!p.background=white
	!p.color=black
	!p.multi=[0,3,2]
print,'at plotboth2A N1,N2 ',N1,N2
 print,' B 2nd ',ALL_PHYS_B[*,1024]
print,' B 2nd ',ALL_PHYS_B[*,1025]
 print,' B 2nd ',ALL_PHYS_B[*,1124]
;	if ihc eq 0 then stop
;
	ihc = 1				; necessary to avoid close 15 lines above
;
if ihcA ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_Hodogram.ps',xsize= 23., ysize=14.2,$
	yoffset=25.5,/landscape,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
for n = 0,size-1 do begin
  APMV = [Aphys[n]-varAx[0],Bphys[n]-varBx[0],Cphys[n]-varCx[0]]
endfor
;
;
     xYmax = max(Bphys[N1:N2],min=xYmin)			; for plot range
     print,'Bphys min, max ',xYmin,xYmax
;
	  nfilt = N2-N1+1
      fltfreq = fltarr(nfilt/2)
      fltfundfr = 1./(delt*nfilt)
	  fltfreq = fltfundfr*indgen(nfilt/2)
;
;	  determine observed frequency
;
      n1f = n1
      n2f = n2
      Awave = fft(Aphys[N1f:N2f])
      Bwave = fft(Bphys[N1f:N2f])
      Cwave = fft(Cphys[N1f:N2f])
      Mwave = fft(Aphys[N1f:N2f])
      Ewave = fft((Bphys[N1f:N2f] - Aphys[N1f:N2f]/xyratio))
      fltfreq = fltarr(nfilt/2)
;
; print,'Awave, N1,N2 ',N1f,N2f
; print,'nfilt ',nfilt
; print,'Asize ',n_elements(Awave)
; print,Awave[*]
;print,'Awave[1:20]',4096.*abs(Awave[1:20])^2
;print,'Bwave[1:20]',4096.*abs(Bwave[1:20])^2
;print,'Cwave[1:20]',4096.*abs(Cwave[1:20])^2
;print,'xyratio ',xyratio
;print,'Ewave[1:20]',4096.*abs(Ewave[1:20])^2
	  Xmax = max(4096.*abs(Awave[1:20])^2,indX)
  	  Ymax = max(4096.*abs(Bwave[1:20])^2,indY)
  	  Zmax = max(4096.*abs(Cwave[1:20])^2,indZ)
  	  Emax = max(4096.*abs(Ewave[1:20])^2,indE)
  	  findex = indX
 print,'Xmax,Y,Z',xmax,Ymax,Zmax,Emax
;
;Mwave[0:nfilt-1] = Awave[0:nfilt-1]
;  	  TpowerM = moment(abs(Aphys[N1f:N2f]))
;  	  Tpower  = (N2f-N1f+1)*TpowerM[1]
  	  if (Ymax gt Xmax) and (IndY ge Indx) then begin
  	  	findex = indY
 ; 	  	TpowerM = moment(abs(Bphys[N1f:N2f]))
 ;	    Tpower  = (N2f-N1f+1)*TpowerM[1]
 	    Mwave = fft(Bphys[N1f:N2f])
 	  endif
  	  if (Zmax gt Xmax) and (Zmax gt Ymax) and (IndZ ge findex) then begin
  	    findex = indZ
  	    Mwave = fft(Cphys[N1f:N2f])
  	  endif
;  	  print,'Mwave[*]',4096.*abs(Mwave[*])^2
  	  TpowerM = moment(abs(Mwave[1:(N2f-N1f-1)]))
  	  Tpower = (N2f-N1f+1)*TpowerM[1]
;  	  print,'Tpower from moment ',Tpower,TpowerM[1]
  	  findex = findex+1
  	  TpowerC = 0.					;
  	  for ncl = 1,nfilt-1 do begin
  	    TpowerC = TpowerC + abs(Mwave[ncl])^2
  	  endfor
;  	  print,'endchk ',abs(Mwave[nfilt-1])^2
;  	  print,'TpowerC ',TpowerC
; 	  print,'N1,N2 ',N1f,N2f
  	  print,'i ',indx,indy,indz,inde
  	  print, 'choice',findex
;		calculate signal power and Bckgpower
	  Sigpower = (Abs(Mwave[findex-1]))^2 + (Abs(Mwave[findex]))^2 + $
	     (Abs(Mwave[findex+1]))^2
;	  Sigpower = 2.*Sigpower
	  Bckgpower = TpowerC - Sigpower
	  print,'Total',TpowerC
	  print,'Signl',Sigpower
	  print,'Bckg ',Bckgpower
	  Bckgpower = Bckgpower>0.
;	  stop
  	  Ypara = fltarr(3)
  	  delx = fltfundfr
  	  ypara[1] = 4096.*abs(Mwave[findex])^2
  	  ypara[0] = 4096.*abs(Mwave[findex-1])^2
  	  ypara[2] = 4096.*abs(Mwave[findex+1])^2
  	  parabola_fit,Ypara,delx,ctr_freq,max
  	  fnyquist = .5/delt
  	  print,'nyquist,findex, fltfundfr ',fnyquist,findex, fltfundfr
;  	  print,'ypara ',ypara[0:2]
  	  ctr_freq = ctr_freq + findex*fltfundfr
  	  print,'ctr_freq ',ctr_freq
;  	  stop
;
;	calculate powers in peak and rest
;	  stop
;
;			Calculate rotation for X,Y, unfiltered data
;			VarAx is moment for unfiltered data, N1 to N2,
 	varAxB = moment(All_Phys_B[0,N1:N2])
 	varBxB = moment(All_Phys_B[1,N1:N2])
 	varCxB = moment(All_Phys_B[2,N1:N2])
;
      Del_Ang 	= 0.
      Angsum	= 0.
      print,'chk ',VarAxB[0],varBxB[0],varCxB[0]
;			In the B system, axis 0 is along B, 2 is in V-B plane
;			so 1,2 is plane perp to B
;		it seems than IDL atan(y,x) is arc tan (y/x)
;		Atan with 2 arguments ranges -pi to pi
;		Have to eliminate jumps at -pi and pi where y changes sign and x is neg
;		Angsum is accumulated rotation, Angsv is last angle, Angxy is this one
;		Ysv is last value of y
      Angsv = atan((All_phys_B[1,N1]-VarBxB[0]),(All_phys_B[0,N1]-VarAxB[0]))
;  	  print,'Dx,Angxy,Del,Sv,tot ',n,Aphys[1]-VarAx[0],Bphys[1]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
      Ysv = All_phys_B[1,N1]-VarBxB[0]
      for n = N1+1,N2 do begin
      print,'N1,N2 ',N1,N2
      print,All_Phys_B[0,n],All_Phys_B[1,n],All_Phys_B[2,n]
      print,All_Phys_B[0,n+1],All_Phys_B[1,n+1],All_Phys_B[2,n+1]
        Angxy = atan((All_PHYS_B[1,n]-VarBxB[0]),(All_phys_B[0,n]-VarAxB[0]))
        Del_ang = Angxy - Angsv
        ynew = All_phys_B[1,N1]-VarBxB[0]
		if (ysv*ynew gt 0.) and ((All_phys_B[0,n]-VarAxB[0]) gt 0.) then $
		    Angsum = Angsum + Del_ang
;	    print,'Dx,Angxy,Del,Sv,tot ',n,Aphys[n]-VarAx[0],Bphys[n]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
	    Angsv = Angxy
	    ysv = ynew
      endfor
      angsum = angsv
     print,'XY rotation ',Angsum,findex,Angsum/findex, ' radians'
;
;		for plot limits
	  xBmax0 = max(All_Phys_B[0,N1:N2],min=XBmin0)
	  xBmax1 = max(All_Phys_B[1,N1:N2],min=XBmin1)
	  xBmax2 = max(All_Phys_B[2,N1:N2],min=XBmin2)
	  print,'Xbmin,XBmax ',XBmin0,XBmax0,XBmin1,XBmax1,XBmin2,XBmax2
;
;	make square plots, all with same range, Hrange
;
  Hrange = (XBmax0-XBmin0)>(XBmax1-XBmin1)>(XBmax2-XBmin2)
  Dhr    = .5*Hrange
  Plim0 = [.5*(XBmax0+xBmin0)-Dhr,.5*(XBmax0+xBmin0)+Dhr]
  Plim1 = [.5*(XBmax1+xBmin1)-Dhr,.5*(XBmax1+xBmin1)+Dhr]
  Plim2 = [.5*(XBmax2+xBmin2)-Dhr,.5*(XBmax2+xBmin2)+Dhr]
;
	  plot,All_Phys_B[0,N1:N2],All_Phys_B[1,N1:N2],xtitle='X Bsys',ytitle='Y Bsys',$
	    charsize=1.8,xrange=Plim0,yrange=Plim1,title=scname,/ynozero,/nodata
	  oplot,All_Phys_B[0,N1:N2],All_Phys_B[1,N1:N2]
	  Ycenter = [varBxB[0],varBxB[0]]
	  Xcenter = [varCxB[0],varCxB[0]]
;	  print,'Ycenter ',Ycenter
;	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,color=red
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(Angsum/findex,format='(F4.1)')
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(a)'
;
;			Calculate rotation for Z,Y, unfiltered data
      Del_Ang 	= 0.
      Angsum	= 0.
;      Angsv = atan((All_phys_B[0,N1:N2]-VarAxB[0]),(All_phys_B[1,N1:N2]-VarBxB[0]))
;		it seems than IDL atan(y,x) is arc tan (y/x) , here y is Z axis, x is ~par Vsw
      Angsv = atan((All_phys_B[2,N1]-VarCxB[0]),(All_phys_B[1,N1]-VarBxB[0]))
      for n = N1+1,N2 do begin
        Angzy = atan((All_phys_B[2,n]-VarCxB[0]),(All_phys_B[1,n]-VarBxB[0]))
        Del_ang =  Angzy - Angsv
        ynew = All_phys_B[1,N1]-VarBxB[0]
		if (ysv*ynew gt 0.) and ((All_phys_B[1,n]-VarBxB[0]) gt 0.) then $
		    Angsum = Angsum + Del_ang
;	    print,'Dx,Angxy,Del,Sv,tot ',n,Aphys[n]-VarAx[0],Bphys[n]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
	    Angsv = Angzy
;	    print,'Angzy,Del,Sv ',ysv*ynew,(All_phys_B[1,n]-VarBxB[0]),Angzy,Del_ang,Angsv
	    ysv = ynew
      endfor
      print,'ZY rotation ',angsum,findex,Angsum/findex, ' radians'
;
	  plot,ALL_Phys_B[1,N1:N2],ALL_Phys_B[2,N1:N2],xtitle='~par Vsw',ytitle='E par B',$
		charsize=1.6,xrange=Plim1,yrange=Plim2,title=timeSC,/ynozero,/nodata
	  oplot,All_phys_B[1,N1:N2],All_Phys_B[2,N1:N2]
	      print,'all_Phys_B el',n_elements(ALl_Phys_B)
    print,'vals ',All_PHys_B[*,N1]
	  Ycenter = [varCxB[0],varCxB[0]]
	  Xcenter = [varBxB[0],varBxB[0]]
;	  print,'Ycenter ',Ycenter
;	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,color=red
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  N2 = N2<size-1
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(Angsum/findex,format='(F4.1)')
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(b)'
;
;	  xrange = [xmin,xmax]
;	  oplot, freq[1:size/2-1], 4096.*abs(Afourx[1:size/2-1])^2,color=red
;	  oplot, freq[1:size/2-1], 4096.*abs(Bfourx[1:size/2-1])^2,color=green
;	  oplot, freq[1:size/2-1], 4096.*abs(Cfourx[1:size/2-1])^2,color=blue
;	  oplot, freq[1:size/2-1], 4096.*abs(Dfourx[1:size/2-1])^2,color=magenta
;   	  oplot, freq[nremove:size/2-1], 4096.*abs(Efourx[1:size/2-1])^2,color=cyan
;	  oplot, freq[nremove:size/2-1], 4096.*abs(Ffourx[1:size/2-1])^2,color=black
;	  oplot, freq[nremove:size/2-1], 4096.*abs(Hfourx[1:size/2-1])^2,color=cyan
;
;	  determine observed frequency
;
      Awave = fltarr(nfilt)
      Bwave = fltarr(nfilt)
      Cwave = fltarr(nfilt)
      Awave = fft(Aphys[N1f:N2f])
      Bwave = fft(Bphys[N1f:N2f])
      Cwave = fft(Cphys[N1f:N2f])
      fltfreq = fltarr(nfilt/2)
;
      fltfundfr = 1./(delt*nfilt)
	  fltfreq = fltfundfr*indgen(nfilt/2)
;
;check plot
;
;	   band pass filter
;
;	The Nyquist frequency is .5/delt
;   peak freq in terms of nyquist
;   ctr/nyquist = Ctr_freq/fnyquist
;    filter from 1/2 center freq to 5*center freq
;    or .5 * findex/nfilt to 5 * findex/nfilt
     flow = .5*ctr_freq/fnyquist
	 fhigh = 5.*ctr_freq/fnyquist
 Ifilt = 1
 if Ifilt ne 0 then begin
   order = 8
   bpfilter = Digital_Filter(flow, fhigh, 40., order)
; the IDL routine does not normalize????.  So do it
; yes it does with ??
   bpsize = n_elements(bpfilter)
;   print,'bpsize ',bpsize
   Sumbp = 0.
   for nbp = 0,2*order do begin
     sumbp = sumbp + bpfilter[nbp]
   endfor
; print,'sumbp ',sumbp
   bpfilter[*] = bpfilter[*]/sumbp
	print,'f ',ctr_freq,flow,fhigh
;   stop
;
;	example   AFphys = Convol(Aphys,hpfilter)
;	print,'bpfilter ',bpfilter
   All_phys_BF = fltarr(3,size)
   Atemp1 = fltarr(size)
   Atemp2 = fltarr(size)
;
;	make All_Phys_B filtered
;
   Atemp2[*] = All_Phys_B[0,0:size-1]
   print,'_B 4th',All_Phys_B[*,1024]
   Atemp1[0:size-1] =  Convol(Atemp2,bpfilter)
   All_Phys_BF[0,0:size-1] = Atemp1[0:size-1]
;   print,'At 1320 ',ATemp1[*]
;
   ;stop
   print,'_BF ',All_Phys_BF[0,0],All_Phys_BF[0,1],All_Phys_BF[0,size-1]
   Atemp2[*] = All_Phys_B[1,0:size-1]
   ATemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_BF[1,0:size-1] = Atemp1[0:size-1]
;
   Atemp2[*] = All_Phys_B[2,0:size-1]
   Atemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_BF[2,0:size-1] = Atemp1[0:size-1]
; endif
;
;	Make RTN filtered
;
   All_phys_RTNF = fltarr(3,size)
;
   Atemp2[*] = All_Phys_RTN[0,0:size-1]
;   print,'_B 4th',All_Phys_B[*,1024]
   Atemp1[0:size-1] =  Convol(Atemp2,bpfilter)
   All_Phys_RTNF[0,0:size-1] = Atemp1[0:size-1]
;   print,'At ',ATemp1[*]
;
   ;stop
;   print,'_BF ',All_Phys_BF[0,0],All_Phys_BF[0,1],All_Phys_BF[0,size-1]
   Atemp2[*] = All_Phys_RTN[1,0:size-1]
   ATemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_RTNF[1,0:size-1] = Atemp1[0:size-1]
;
   Atemp2[*] = All_Phys_RTN[2,0:size-1]
   Atemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_RTNF[2,0:size-1] = Atemp1[0:size-1]
;
 	varAxRTNF = moment(All_Phys_RTNF[0,N1:N2])
 	varBxRTNF = moment(All_Phys_RTNF[1,N1:N2])
 	varCxRTNF = moment(All_Phys_RTNF[2,N1:N2])
;
;	plot the three waveforms in one panel
	xspmin = float(floor(tsamp1))
	xspmax = float(floor(Tsamp2)) + 1.
	  plot,timeax[N1:N2],1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),$
	  title='samples '+string(N1,format='(I5)')+' to'+string(N2,format='(I5)')$
	  ,xrange=[xspmin,xspmax],charsize=2.0,ytitle='mV',$
	    xtitle='sec',/nodata
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),color=red
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[1,N1:N2]-varBxRTNF[0]),color=blue
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[2,N1:N2]-varCxRTNF[0]),color=green
;
;		put antenna designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .07*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'R',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'T',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'N',color=green
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(c)'
;
;
  varAxF = moment(ALL_PHys_BF[0,N1:N2])
  varBxF = moment(All_PHys_BF[1,N1:N2])
  varCxF = moment(All_Phys_bF[2,N1:N2])
;
;			Calculate rotation for X,Y, filtered data
;				for plot limits
	print,'for plot limit at 1395 N1,N2= ',N1,N2
	  xBFmax0 = max(All_Phys_BF[0,N1:N2],min=XBFmin0)
	  xBFmax1 = max(All_Phys_BF[1,N1:N2],min=XBFmin1)
	  xBFmax2 = max(All_Phys_BF[2,N1:N2],min=XBFmin2)
;
;
;	make square plots, all with same range, Hrange
;
  Frange = abs(XBFmax0-XBFmin0)>abs(XBFmax1-XBFmin1)>abs(XBFmax2-XBFmin2)
  Dhr    = .5*Frange
  PFlim0 = [.5*(XBFmax0+xBFmin0)-Dhr,.5*(XBFmax0+xBFmin0)+Dhr]
  PFlim1 = [.5*(XBFmax1+xBFmin1)-Dhr,.5*(XBFmax1+xBFmin1)+Dhr]
  PFlim2 = [.5*(XBFmax2+xBFmin2)-Dhr,.5*(XBFmax2+xBFmin2)+Dhr]
;
 print,'nfilt ',nfilt
; print,'_BF 2 ',All_phys_BF[0:nfilt+2*order]
 ;stop
  	print,'moments[0] ',VarAxB[0],VarBxB[0],VarCxB[0]
      Del_Ang 	= 0.
      Angsum	= 0.
;      print,'chk s',VarAxF[0],varBxF[0],varCxf[0]
;		it seems than IDL atan(y,x) is arc tan (y/x)
      Angsv = atan((ALL_phys_BF[1,N1]-varBxF[0]),(All_Phys_BF[0,N1]-varAxF[0]))
      Ysv = (ALL_phys_BF[1,N1]-varBxF[0])
;      for n = order,nfilt+order do begin
;	stuff_plot = fltarr(??) is plotted by plotboth4
	  stuff_plot = fltarr(N2-N1)		; first point not plotted
      for n = N1+1,N2 do begin
        Angxy = atan((ALL_PHYS_BF[1,n]-varBxF[0]),(ALL_PHYS_BF[0,n])-varAxF[0])
        Del_ang = Angxy - Angsv
        ynew = (ALL_phys_BF[1,n]-varBxF[0])
        if (ysv*ynew gt 0.) or ((ALL_PHYS_BF[0,n]-varAxF[0]) gt 0.) then $
        Angsum = Angsum+Del_Ang
;        if (ALL_PHYS_BF[2,n]-VarCxF[0])*(ALL_phys_BF[1,n]-VarBxF[0]) gt 0. $
;        	then Angsum = Angsum+Del_Ang
;	    print,'Dx,Angxy,Del,Sv,tot ',n,Aphys[n]-VarAx[0],Bphys[n]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
;		print,'Dx,Angxy,Del,Sv,tot ',n,ynew,All_phys_BF[0,n]-varAxF[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
;		stuff_plot[n-N1-1] = Angsum
	    Angsv = Angxy
	    ysv = ynew
      endfor
;     print,'XY rotation ',Angsum/findex, ' radians'
	  cyclecnt = (N2-N1+1)*delt*ctr_freq
;      rotxy = angsum/findex
      rotxy = angsum/cyclecnt
;
 	  nstsp = nfilt/20					; NO OF POINTS TO COLOR for start and finish

	  plot,All_Phys_BF[0,N1:N2],All_Phys_BF[1,N1:N2],xrange=PFLIM0,$
	  	yrange=PFLIM1,xtitle='X Bsys',ytitle='Y Bsys',charsize=1.8
	  oplot,All_Phys_BF[0,N1:N1+nstsp],All_Phys_BF[1,n1:N1+nstsp],color=green,thick=2
	  oplot,All_Phys_BF[0,N2-nstsp:N2],$
	     ALL_Phys_BF[1,N2-nstsp:N2],color=red,thick=2
	  Ycenter = [varBxF[0],varBxF[0]]
	  Xcenter = [varAxF[0],varAxF[0]]
	  print,'Ycenter ',Ycenter
	  print,'Xcenter ',Xcenter
	  oplot,Xcenter,Ycenter,psym=4,color=red
  	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(rotxy,format='(F6.1)')
  	  xplace = .2*!x.crange[0] + .8*!x.crange[1]
	  yplace = .9*!y.crange[0] + .1*!y.crange[1]
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(d)'
	  Binout = Binoutsv
	  if Binout eq 1 then xyouts,xplace,yplace,'B in'
	  if Binout eq 0 then xyouts,xplace,yplace,'B out'
;	  add 2 if angle B to R is greater than 63 deg
	  Btrans = sqrt(Bnow[0]^2 + Bnow[1]^2)
	  if abs(Bnow[2]) gt 2.*Btrans then Binout = Binout + 2
;

;			Calculate rotation for Z,X, filtered data
      Del_Ang 	= 0.
      Angsum	= 0.
;      print,'chk ',VarAx[0],varBx[0],varCx[0]
;		it seems than IDL atan(y,x) is arc tan (y/x)
      Angsv = atan((All_Phys_BF[2,N1]-varCxF[0]),(All_phys_BF[0,N1]-varAxF[0]))
      Ysv = (All_Phys_BF[2,N1]-VarBxF[0])
      for n = order+1,nfilt+order do begin
        Angxy = atan((All_PHys_BF[2,n]-varCxF[0]),(ALL_phys_BF[0,n]-varAxF[0]))
        Del_ang = Angxy - Angsv
        ynew = (All_Phys_BF[2,n]-varCxF[0])
;        if (All_phys_BF[1,n]-VarBxF[0])*(All_phys_BF[0,n-1]-VarAxF[0]) gt 0. then $
		 if (ysv*ynew gt 0.) and ((All_Phys_BF[0,n]-varAxF[0]) gt 0.) then $
        	Angsum = Angsum+Del_Ang
;	    print,'Dx,Angxy2Del,Sv,tot ',n,Aphys[n]-VarAx[0],Bphys[n]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
	    Angsv = Angxy
	    ysv = ynew
      endfor
;
;      print,'ZY rotation ',Angsum/findex, ' radians'
;      rotzy = angsum/findex
      rotzy = angsum/cyclecnt
	  plot,All_Phys_BF[0,N1:N2],All_phys_BF[2,N1:N2],xtitle='X Bsys',$
	  ytitle='E par B',charsize=1.8,xrange=PFLIM0,yrange=PFLIM2,$
	  title='input hhmmss '+string(hhmmssst,format='(I6)'),/ynozero,/nodata
	   oplot,All_Phys_BF[0,N1:N2],All_phys_BF[2,N1:N2]
	   oplot,All_phys_BF[0,N1:N1+nstsp],All_phys_BF[2,N1:N1+nstsp],$
	   color=green,thick=2
	   oplot,All_phys_BF[0,N2-nstsp:N2],$
	     All_phys_BF[2,N2-nstsp:N2],color=red,thick=2
	  Xcenter = [varAxF[0],varAxF[0]]
	  Ycenter = [varCxF[0],varCxF[0]]
;	  print,'Ycenter ',Ycenter
;	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,color=red
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(Angsum/findex,format='(F6.1)')
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(e)'
;
;	plot B instead of spectrum
;
;	for the files from CDAWeb
;		already loaded,  now select corresponding samples
	print,' mag,APM time ',scetmag,scet,format='(A18,2F14.6)'
	print,'N1,N2 ',N1,N2
	print,'tsamp1,2 ',tsamp1,tsamp2
;
;		find mag samples corresponding to selection
;
;	Tsamp1 is in sec = 1./86400 of aday
;		Tseries is in days
;		Btime is in days
;	1 msec = 1./8.64D07 of a day
    print,'Tseries[0],Btime[0] ',Tseries[0],Btime[0]
	scetmag = 0.d00
    scett = Tseries[0] + Tsamp1/8.6400D04			; day
    nsel = 0
	while SCETMAG LT scett do begin
	  print,' '
	  scetmag = Btime[nsel]
	  print,'going to get Btime, scetmag,scett ',scetmag,scett
;	  print,'tchk ',hr,mmin,sec1,usec
;	  print,'mag chk ',scetmag-scet,mmin,denst,tideg,Bnowt
	  nsel = nsel+1
	ENDwhile
	print,'at 1542,search Btime select, nsel= ',nsel
	print,'scet,mag,diff ',scett,scetmag,scetmag-scett,format='(A16,3F16.8)'
;	load B data, assumed to be 8 samples per sec so 256 samples
	Nsel2 = Nsel + (N2-N1)/8
;
 print,'nsel,Btime[0],Btime[nsel] ',nsel,Btime[0],Btime[nsel]
;  N1 and N2 are the first and last APM sample numbers for selection
;  nsel and nsel2 are same for B samples in the selection
 print,'going to plot times ',Btimepl[nsel],Btimepl[nsel2]
 for npl = nsel-1,nsel2 do begin
   Btimepl[npl] = (Btime[npl]-Tseries[0])*8.64D04 			; days to sec
;   print,'pl ',npl,Btimepl[npl],format='(A10,I5,F14.8)'
 endfor
 print,'BtinePL ',Btimepl[nsel-1],Btimepl[nsel2]
 print,'got to 1483,nsel,nsel2 ', nsel,nsel2
;	NBload is the number of B samples in the full event, 32 sec
;   BtimePL
;
 	varAxBev = moment(Bevent[0, 0:NBload-1])
 	varBxBev = moment(Bevent[1, 0:NBload-1])
 	varCxBev = moment(Bevent[2, 0:NBload-1])
 	varAxBsm = moment(Bevent[0,nsel:nsel2])
 	varBxBsm = moment(Bevent[1,nsel:nsel2])
 	varCxBsm = moment(Bevent[2,nsel:nsel2])
	Thismax = max(Bevent[0:2,Nsel-1:nsel2],min=thismin)
	plot,BtimePL[ 0:NBload-1],BEVENT[1, 0:NBload-1]-VarBxBev[0],$
	yrange=[0,0.],xrange=[xspmin,xspmax],xtitle='sec',$
	ytitle='B (nT)',charsize=2,/nodata
	oplot,BtimePL[ nsel:nsel2],BEVENT[0, nsel:nsel2]-VarAxBev[0],color=red
	oplot,BtimePL[ nsel:nsel2],BEVENT[1, nsel:nsel2]-VarBxBev[0],color=blue
	oplot,BtimePL[ nsel:nsel2],BEVENT[2, nsel:nsel2]-VarCxBev[0],color=green
	xplace = .15*!x.crange[0] + .85*!x.crange[1]
	xyouts,xplace,yplace,'(f)'
;
	timeS = systime(0)
	xyouts, 0., 0.,scname+' APM_Hodogram.ps, from APM_Analysis_Poynting, version '$
	+version+ ',  output at '+ timeS,charsize=.7,/normal
;
if ihcA ne 1 then begin
	device, /close
	set_plot, mydevice
endif
;
;**************
; ihc = 1
; goto, plotboth4
;**************
ihcA = 1-ihcA
;
print,'got to read dispose '
print,'IhcA ',ihcA
;stop
;if ihcA eq 1 then read,dispose
if ihcA ne 1 then goto,plotboth2A
;

	err = TM_UR8_to_Ydoyh(Scet,yyyy,doy,hh)
	err = TM_UR8_to_YMD(Scet,yyyy,mmon,dd,hh,min,sec,msec)
	yyyymmddp = 10000*yyyy + 100*mmon + dd
;	print,'yyyy chk ',yyyymmddp
	hhmmssp = 10000*HH + 100*MIN + SEC
	SCETI4 = [yyyymmddp,hhmmssp]
;
 print,'min,denst ',nmin,denst
;
;now I need eigenvectors, ratio of eigenevalues, direction of principal
;		eigenvector wrt B

    eigenvectAnt 	= fltarr(3,3)
    eigenvalAnt	= fltarr(3)
    print,'at 1594 Vmatrix, size =',size
;    	print,all_apm[2,0:4],all_apm[2,size-1]
	vmatrix, All_APM, N1, N2, eigenvalAnt, eigenvectAnt
	print,'Ant eigvals ',eigenvalAnt
	print,'principal in Ant coords ',eigenvectAnt[*,0]
	print,'   next   in Ant coords ',eigenvectAnt[*,1]
;
    eigenvectSC 	= fltarr(3,3)
    eigenvalSC	= fltarr(3)
;    print,'at Vmatrix, size =',size
;    	print,all_apm[2,0:4],all_apm[2,size-1]
	vmatrix, All_Phys_SC, N1, N2, eigenvalSC, eigenvectSC
	print,'SC eigvals ',eigenvalSC
	print,'principal in SC coords ',eigenvectSC[*,0]
	print,'   next   in SC coords ',eigenvectSC[*,1]
;
    eigenvectB 	= fltarr(3,3)
    eigenvalB	= fltarr(3)
    print,'all_Phys_B el',n_elements(ALl_Phys_B)
    print,'vals ',All_PHys_B[*,N1]
	vmatrix, All_Phys_B, N1, N2, eigenvalB, eigenvectB
;	vmatrix, All_Phys_B, N1, N2, eigenvalB, eigenvectB
	print,'B eigvals ',eigenvalB
	print,'principal in B coords ',eigenvectB[*,0]
	print,'   next   in B coords ',eigenvectB[*,1]
	Eperp = sqrt(EigenvalB[1] + EigenvalB[2])
	eccent = Eperp/sqrt(eigenvalB[0])
;
    eigenvectRTN	= fltarr(3,3)
    eigenvalRTN	= fltarr(3)
	vmatrix, All_Phys_RTN, N1, N2, eigenvalRTN, eigenvectRTN
	print,'RTN eigvals ',eigenvalRTN
	print,'principal in RTN coords, chosen ',eigenvectRTN[*,0]
	print,'   next   in RTN coords, chosen ',eigenvectRTN[*,1]
;
;  transform principal eigenvector from antenna potential
;		to electric field in antenna system
;   in eigenvect, first index is component
; calculate dot product OF E with B
  EdotBRTN = (eigenvectRTN[0,0]*bnorm[0]+EigenvectRTN[1,0]*bnorm[1]+$
  	EigenvectRTN[2,0]*bnorm[2])
print,'EdotB RTN',EdotBRTN
;
 ESC = fltarr(3)
 Vtemp = fltarr(3)
 ERTN = fltarr(3)
 ESC = eigenvectAnt[*,0]
 XformSC, sc, ctr_freq, denst, invmatrix
 Vtemp[*] = invmatrix#ESC	     ;Vtemp is now principal Eigv in SC system
 Ertn[*] = Att_Matrix#VTemp
 Vtempnorm = sqrt(Vtemp[0]^2+Vtemp[1]^2+Vtemp[2]^2)
 ERTnorm = sqrt(ERTN[0]^2+ERTN[1]^2+ERTN[2]^2)
;print,'B rtn ',bnorm
;print,'ESCnorm ',sqrt(ESC[0]^2+ESC[1]^2+ESC[2]^2)
;print,'Vtempnorm ',Vtempnorm
;print,'ERTNorm ',ERTnorm
;
; calculate dot product OF E with B
  EdotB1 = (ERTN[0]*bnorm[0]+ERTN[1]*bnorm[1]+ERTN[2]*bnorm[2])/ERTnorm
print,'EdotB ',EdotB1
;
;antenna_test
;
; print,'stuff for checking '
;    print,'Bnorm in RTN ',Bnorm
;    print,'Bnorm in Ant ',
    eigenvectSC	= fltarr(3,3)
    eigenvalSC	= fltarr(3)
	vmatrix, All_APM, N1, N2, eigenvalSC, eigenvectSC
 	print,'principal in SC coords ',eigenvectSC[*,0]
	print,'   next   in SC coords ',eigenvectSC[*,1]
; eigenvalB[1]/eigenvalB[0],
 print,'check 39 ',resultfilename
 print,'check 41 ',result2filename
 print,'check N1,N2 ',N1,N2
 Eenergy = 0.
 Benergy = 0.
 Duration = (Tseries[N2]-Tseries[N1])
 for nen = N1,N2 do begin
   Eenergy = Eenergy + (All_phys_RTNF[0,nen]-varAxRTNF[0])^2
   Eenergy = Eenergy + (All_phys_RTNF[1,nen]-varBxRTNF[0])^2
   Eenergy = Eenergy + (All_phys_RTNF[2,nen]-varCxRTNF[0])^2
 endfor
 Eenergy = .5*eps0*Eenergy
 for nen = nsel,nsel2 do begin
   Benergy = Benergy + (BEVENT[0,nen]-VarAxBsm[0])^2
   Benergy = Benergy + (BEVENT[1,nen]-VarBxBsm[0])^2
   Benergy = Benergy + (BEVENT[2,nen]-VarCxBsm[0])^2
 endfor
; E is V/m but B is nT
 Benergy = .5*1.e-18*Benergy/mu0
 print,'E.B,energy ',Eenergy,Benergy
;
;	calculate average Poynting vector
;		find APM samples matching B samples
;
  PoYAvr = fltarr(3)
  TempE		= fltarr(3)
  TempB		= fltarr(3)
; now nsel-1 and nsel2 are the limits on the Binput data
;	I need the corresponding limits on the APM input, N1, N2
;  find APM data corresponding to B data
;
  NPplot = nsel2-nsel+2
  poynting = fltarr(3,NPplot)
  for nPoY = nsel-1,nsel2 do begin
    nselP = nPoY
    scetmag = Btime[nPoY]
    print,'SCETT,scetmag ',SCETT,scetmag
	while SCETT LT scetmag do begin
;	  print,' '
	print,'nselP ',nselP,nsel,nsel2
;	  print,'going to get Btime, scetmag,scett ',scetmag,scett,format='(A30,2F14.7)'
	  scett = Tseries[nselP] 						; day
	  nselP = nselP+1
	ENDwhile
;	BEVENT[0,NSEL-1:nsel2]-VarAxBev[0]

;    TempB[0] = Bevent[0,nPoY]-VarAxBev[0]
;    TempB[1] = Bevent[1,nPoY]-VarBxBev[0]
;    TempB[2] = Bevent[2,nPoY]-VarCxBev[0]
    TempB[0] = Bevent[0,nPoY]-VarAxBsm[0]
    TempB[1] = Bevent[1,nPoY]-VarBxBsm[0]
    TempB[2] = Bevent[2,nPoY]-VarCxBsm[0]
    TempE[0] = All_phys_RTNF[0,nselP]-varAxRTNF[0]
    TempE[1] = All_phys_RTNF[1,nselP]-varBxRTNF[0]
    TempE[2] = All_phys_RTNF[2,nselP]-varCxRTNF[0]
    print,'***',nselP,TempB[*],	(Btime[nPoY]-Tseries[0])*8.64D04
    PoYAvr = crossp(tempE,TempB)
;    print,'Poynting ',PoYavr
    Poynting[*,nPoY-nsel+1] = PoYAvr[*]
;    print,npoy-nsel+1,PoYAvr
  endfor
  print,'var APM ',varAxRTNF[0]
  print,'var B   ',varAxBev[0]
;
 crd_sys = 'RTN'
; I need:  SCETI4, Eenergy,Benergy, Eforce, Bforce, PoyntingV(3),Binout,Rotxy
;  also proton speed??, obliquity,  PoyntingV.dot.B
; SCETI$ is 16, 7 more at 12.3 is 84, 6 more, so need two files  (106 total)
 openw,39,resultfilename,/append
   printf,39,SCETI4,ctr_freq,Tot_avr[0:2],1000.*Arms,1000.*Brms,1000.*Crms,$
   findex,rotxy,rotzy,N1,N2,$
   binout,format='(I9,I7.6,F6.2,3F5.2,3F6.2,I3,2F6.1,2I5,I2)'
 openw,41,result2filename,/append
    printf,41,SCETI4,eccent,EdotB1,EdotBRTN,VSWt,xyratioAPM,zyratioAPM,denst,$
    eigenvalB[1]/eigenvalB[0],eigenvalB[2]/eigenvalB[0],sqrt(Bckgpower/Sigpower)$
    ,crd_sys,format='(I9,I7.6,F6.3,2F7.3,F6.0,2f7.3,f6.2,2F6.3,F5.2,A4)'
 close,39
 close,41
;
 erase
 ihc = 1
; goto,ninestrip
;
ihcA = 1			; here 1 is screen
;
;stop
;
ihc = 1			; here 1 is screen
;
plotboth7: print,'ihc = ',ihc
	!p.background=white
	!p.color=black
	!p.multi=[0,1,1]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='Alias_detail.ps',xsize= 18., ysize=18.
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
	lf = size/2-50
	  plot, freq[lf:size/2-1], 4096.*abs(Bfourx[lf:size/2-1])^2,yrange=[0.,0.],$
	  xrange=[0.,0.]
;	  oplot, freq[1:size/2-1], 4096.*abs(Bfourx[1:size/2-1])^2,color=blue
	timeS = systime(0)
	xyouts, 0., -.03,scname+' timeseq.ps, version '+version+ ',  output at '+ timeS,$
charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth7
;
;stop
;
 print,'N1,N2 at 1771 ',N1,N2
; stop
plotboth3: print,'at plotboth3, ihc = ',ihc
;  printf,66,'at plotboth3, scet= ',scet
	!p.background=white
	!p.color=black
	!p.multi=[0,2,2]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='LRS_Corr.ps',xsize= 18., ysize=18.
endif
;
;		PLOT correlations, Y-X, Z-X, FOR APM AND LRS SEPARATELY
;
;
	print,'starting lrs_corr.ps on screen'
;
	!x.range = [0.,0.]
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE
;
;	N1 = 0
;	N2 = 2047
	eigenval2 	= fltarr(2)
	eigenvect2 	= fltarr(2,2)
;
	time1 = n1*delt
	time2 = n2*delt
	toptitle = SCname+'  SCET '+timeSC+string(time1,format='(F5.1)')+$
	string(time2,format='(F5.1)')
	plot,Bphys[n1:n2],Aphys[n1:n2],xtitle='Y APM',ytitle = 'X APM',$
	Title=toptitle,	xrange = [0.,0.],$
	charsize=1., psym = 4, symsize=.4,/ynozero
;
	vmatrix2, Bphys, Aphys, N1, N2, eigenval2, eigenvect2
;	print,'B,A eigenvalues',eigenval2
;	print,'B,A, eigenvect',eigenvect2
;	print,'B/A ratio',eigenvect2[0,0]/eigenvect2[1,0]
;	print,'avr,B,A',varBx[0],varAx[0]
	xe1 = -sqrt(eigenval2[0])*eigenvect2[0,0] + varBx[0]
	xe2 =  sqrt(eigenval2[0])*eigenvect2[0,0] + varBx[0]
	ye1 = -sqrt(eigenval2[0])*eigenvect2[1,0] + varAx[0]
	ye2 =  sqrt(eigenval2[0])*eigenvect2[1,0] + varAx[0]
	xeigAB = [xe1,xe2]
	yeigAB = [ye1,ye2]
;	print,'xeigAB',xeigAB
;	print,'yeigAB',yeigAB
	oplot,xeigAB,yeigAB,color=black
	xplace = .9*!x.crange[0] + .1*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'X/Y '+string(eigenvect2[1,0]/eigenvect2[0,0])
;	if n2 ne 0 then stop
;
;	 xstyle=17, charsize=1., psym = 4, symsize=1.,/ynozero
;
	BCTitle = 'N1, N2'+string(N1,format='(I6)')+string(N2,format='(I6)')
	plot,Bphys[n1:n2],Cphys[n1:n2],xtitle='Y APM',ytitle = 'Z APM',$
	title = BCtitle,charsize=1.,psym = 4, symsize=.4,xrange = [0.,0.],/ynozero
;
	vmatrix2, Bphys, Cphys, N1, N2, eigenval, eigenvect
;	print,'B,C eigenvalues',eigenval
;	print,'B,C, eigenvect',eigenvect
;	print,'avr,B,C',varBx[0],varCx[0]
	xe1 = -sqrt(eigenval[0])*eigenvect[0,0] + varBx[0]
	xe2 =  sqrt(eigenval[0])*eigenvect[0,0] + varBx[0]
	ye1 = -sqrt(eigenval[0])*eigenvect[1,0] + varCx[0]
	ye2 =  sqrt(eigenval[0])*eigenvect[1,0] + varCx[0]
	xeigBC = [xe1,xe2]
	yeigBC = [ye1,ye2]
;	print,'xeigBC',xeigBC
;	print,'yeigBC',yeigBC
	oplot,xeigBC,yeigBC,color=black
	xplace = .9*!x.crange[0] + .1*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'Z/Y '+string(eigenvect[1,0]/eigenvect[0,0])
;
;	EDTitle = string(' sec'
	plot,Ephys[n1:n2],Dphys[n1:n2],xtitle='Y LRS',ytitle = 'X LRS',$
	charsize=1. ,psym = 4, symsize=.4,xrange = [0.,0.],/ynozero
;
	vmatrix2, Ephys, Dphys, N1, N2, eigenval, eigenvect
	print,'in 1852 ',N1,N2
;;	stop
;	print,'E,D eigenvalues',eigenval
;	print,'E,D, eigenvect',eigenvect
;	print,'avr,E,D',varEx[0],varDx[0]
	xe1 = -sqrt(eigenval[0])*eigenvect[0,0] + varEx[0]
	xe2 =  sqrt(eigenval[0])*eigenvect[0,0] + varEx[0]
	ye1 = -sqrt(eigenval[0])*eigenvect[1,0] + varDx[0]
	ye2 =  sqrt(eigenval[0])*eigenvect[1,0] + varDx[0]
	xeigED = [xe1,xe2]
	yeigED = [ye1,ye2]
;	print,'xeigED',xeigED
;	print,'yeigED',yeigED
	oplot,xeigED,yeigED,color=black
	xplace = .9*!x.crange[0] + .1*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'X/Y eig'+string(eigenvect[1,0]/eigenvect[0,0])
	yplace = yplace - .1*(!y.crange[1] - !y.crange[0])
	xyouts,xplace,yplace,'X/Y pwr'+string(alrsxy)
;
	plot,Ephys[n1:n2],Fphys[n1:n2],xtitle='Y LRS',ytitle = 'Z LRS',$
	charsize=1. ,psym = 4, symsize=.4,xrange = [0.,0.],/ynozero
;
	vmatrix2, Ephys, Fphys, N1, N2, eigenval, eigenvect
;	print,'E,F eigenvalues',eigenval
;	print,'E,F, eigenvect',eigenvect
;	print,'avr,E,F',varEx[0],varFx[0]
	xe1 = -sqrt(eigenval[0])*eigenvect[0,0] + varEx[0]
	xe2 =  sqrt(eigenval[0])*eigenvect[0,0] + varEx[0]
	ye1 = -sqrt(eigenval[0])*eigenvect[1,0] + varFx[0]
	ye2 =  sqrt(eigenval[0])*eigenvect[1,0] + varFx[0]
	xeigEF = [xe1,xe2]
	yeigEF = [ye1,ye2]
;	print,'xeigEF',xeigEF
;	print,'yeigEF',yeigEF
	oplot,xeigEF,yeigEF,color=black
	xplace = .9*!x.crange[0] + .1*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'Z/Y eig'+string(eigenvect[1,0]/eigenvect[0,0])
	yplace = yplace - .1*(!y.crange[1] - !y.crange[0])
	xyouts,xplace,yplace,'Z/Y pwr'+string(alrszy)
;
	timeS = systime(0)
	xyouts, 0., -.03,scname+' LRS_Corr.ps, from APM_Analysis,version '+version+$
	 ',  output at '+ timeS,charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
;	stuff_plot = fltarr(??) is plotted by plotboth4 ;
;  printf,66,'after plotboth3, scet= ',scet
;
 endif
;
ihc = 1			; here 1 is screen
;
plotboth4: print,'ihc = ',ihc
	!p.background=white
	!p.color=black
	!p.multi=[0,1,1]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='Miscellany.ps'
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
	nfilt = n_elements(Stuff_plot)
	print,'nfilt = ',nfilt
	plot, Stuff_plot,xrange=[0,nfilt-1]
;
	timeS = systime(0)
	xyouts, 0., .0,scname+' miscellany.ps, version '+version+ ',  output at '+ timeS,$
charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth4
;
;stop
 goto, ninestrip
;
;
ihc = 0			; here 1 is screen but do hardcopy first
;
plotboth5: print,'at plotboth5, ihc = ',ihc
;
;	now 3 of 9 strip screen plot plus x-y
;
	!p.background=white
	!p.color=black
	!p.multi=[0,1,4]
	y0 = .77
	y1 = .95
	dy = 1.2*(y1-y0)
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_ninestrip_pub.ps',xsize= 25.,ysize=18.5,$
	/landscape,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;		get density response data for calculating response to density fluct
;
	alrsxy = 1.05
	alrszy = 1.
	ratio_data_name = 'APM_antenna_ratios_B_Carr2140.txt'
	openr,37,Ratio_data_name
	readf,37,n_entry_no
;	print,'entry no ',n_entry_no
	APM_density = fltarr(n_entry_no)
	xyAPMratio = fltarr(n_entry_no)
	ZyAPMratio = fltarr(n_entry_no)
	  readf,37,APM_Density
;	  print,APM_density
	  readf,37,xyAPMratio
	  readf,37,zyAPMratio
	close,37
	xyratio = interpol(xyAPMratio,APM_density,denst)
	zyratio = interpol(zyAPMratio,APM_density,denst)
;	print,'denst ',denst
;	print,'denst,xyratio ',denst,xyratio
;	print,'density',APM_density
;	print,'zy ',zyAPMratio
;
  for np = 0,3 do begin
	  plottitle = ''
	  if np eq 0 then begin
		Pphys = Aphys[0:2047]
		plottitle = scname+' '+timeSC
		plcolor = red
;		print,'at 570 ihc,np,Aphys[1] ',ihc,np,Aphys[1]
      endif
;	  if ihc eq 1 then stop
	  if np eq 0 or np eq 3 then plcolor = red
	  if np eq 1 then Pphys = Bphys[0:2047]
	  if np eq 1 then Pseries = Bseries[0:2047]
	  if np eq 1 or np eq 4 then plcolor = blue
	  if np eq 2 then Pphys = Cphys[0:2047]
	  if np eq 2 then Pseries = Cseries[0:2047]
	  if np eq 2 or np eq 5 then plcolor = green
	  if np eq 3 then plcolor = magenta
;
;	  if np eq 0 then begin
;	    print,'starting np = 0 check '
;	    print,'ihc ',ihc
;		print,'title ',plottitle
;		print,'pphys ',Pphys[0],Pphys[2047]
;	    print,'time ',timesc
;	  endif
;
	ymax = max(-Pphys[0:size-1],min = ymin)<1.5
	Ymin = ymin>(-1.)
	if np ne 3 then plot,timeax,-Pphys,charsize=2.2,Ytitle=ant_char[np],$
	  title=plottitle,xstyle=17,position=[.1, y0, .99, y1],yrange=[ymin,ymax],$
	  xticklen=.1,xtickname=nulltick,thick=2,/nodata,/ynozero
;
;	if np eq 3 then Y?
	if np eq 3 then Pphys = Bphys[0:2047]-Aphys[0:2047]/xyratio
	if np eq 3 then plot,timeax[0:2047],Pphys[0:2047],charsize=2.2,$
	  Ytitle='x-y dipole',xstyle=17,position=[.1, y0, .99, y1],$
	  xticklen=.1,/nodata,/ynozero
;
	  oplot,timeax,Pphys[0:size-1],color=plcolor
	  var = moment(float(Pseries[0:size-1]), sdev=ctrms)
		print,'rms in counts ',ctrms
	  var = moment(Pphys[0:size-1], sdev=vrms)
	  prms = string(format='(f8.1,a4,f7.2,a3)',ctrms,' cts',$
		1000.*vrms,' mV')
	  xyouts, .70, .2*y0+.8*y1, prms,/normal
;
	  xline = [tsamp1,tsamp1]
	  yline = !Y.Crange
;	  print,'xline,yline ',xline,yline
	  oplot, xline,yline,color=red
	  xline = [tsamp2,tsamp2]
	  yline = !Y.Crange
;	  print,'xline,yline ',xline,yline
	  oplot, xline,yline,color=red
;
	  y0 = y0-dy
	  y1 = y1-dy
	  if np eq 2 then begin
	  	y0 = y0-.2*dy
	    y1 = y1-.2*dy
	  endif
  endfor
;
 ;
	timeS = systime(0)
	xyouts,0.,0.,scname+'APM_ninestrip_pub.ps, from APM_Analysis_Poynting'$
	+' version '+version+$
 ',  output at '+ timeS,charsize=.7,/normal
; 	      stop	OK
;
	print,'finished plotboth5, ihc= ',ihc
;
if ihc ne 1 then begin
	device,/close
	set_plot, mydevice
endif
ihc = 1-ihc
;
 if ihc ne 0 then goto,plotboth5
;

;
ihc = 1			; here 1 is screen, do first
;
plotboth6: print,'at plotboth6, ihc = ',ihc
;
;	now 3 of 9 strip screen plot plus 3 components of magnetic field
;
	!p.background=white
	!p.color=black
	!p.multi=[0,1,4]
	y0 = .77
	y1 = .95
	dy = 1.2*(y1-y0)
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_ninestrip_Mag_pub.ps',xsize= 25.,ysize=18.5,$
	/landscape,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;		get antenna ratio data for interpolation
	alrsxy = 1.05
	alrszy = 1.
	ratio_data_name = 'APM_antenna_ratios_B_Carr2140.txt'
	openr,37,Ratio_data_name
	readf,37,n_entry_no
;	print,'entry no ',n_entry_no
	APM_density = fltarr(n_entry_no)
	xyAPMratio = fltarr(n_entry_no)
	ZyAPMratio = fltarr(n_entry_no)
	  readf,37,APM_Density
;	  print,APM_density
	  readf,37,xyAPMratio
	  readf,37,zyAPMratio
	close,37
	xyratio = interpol(xyAPMratio,APM_density,denst)
	zyratio = interpol(zyAPMratio,APM_density,denst)
	print,'denst,xyratio ',denst,xyratio
;	print,'density',APM_density
;	print,'zy ',zyAPMratio
;
;	determine section of B data corresponding to APM event
;
    B_pub_cnt = 0
    scet = Tseries[0]
;
	for np = 0,2 do begin
	  plottitle = ''
	  if np eq 0 then begin
		Pphys = Aphys[0:2047]
		plottitle = scname+' '+timeSC
		plcolor = red
		print,'at 570 ihc,np,Aphys[1] ',ihc,np,Aphys[1]
  	  endif
;	  if ihc eq 1 then stop
	  if np eq 0 or np eq 3 then plcolor = red
	  if np eq 1 then Pphys = Bphys[0:2047]
	  if np eq 1 then Pseries = Bseries[0:2047]
	  if np eq 1 or np eq 4 then plcolor = blue
	  if np eq 2 then Pphys = Cphys[0:2047]
	  if np eq 2 then Pseries = Cseries[0:2047]
	  if np eq 2 or np eq 5 then plcolor = green
;
	print,'******* ',np
	  if np eq 0 then begin
	    print,'starting np = 0 check '
	    print,'ihc ',ihc
		print,'title ',plottitle
		print,'pphys ',Pphys[0],Pphys[2047]
	    print,'time ',timesc
	  endif
;
	ymax = max(Pphys[0:size-1],min=ymin)<1.5
	Ymin = ymin>(-1.)
	print,'y ',ant_char[np]
	plot,timeax,Pphys,charsize=2.2,Ytitle=ant_char[np],$
	  title=plottitle,xstyle=17,position=[.1, y0, .99, y1],yrange=[ymin,ymax],$
	  xticklen=.1,xtickname=nulltick,thick=2,/ynozero,/nodata
	oplot,timeax,Pphys,color=plcolor
;
      y0 = y0-dy
	  y1 = y1-dy

  endfor

	  y0 = y0-.2*dy
	  y1 = y1-.2*dy
;
	  print,'Tseries[0] ',Tseries[0],format='(A11,F14.7)'
	  print,'Btime[0]   ',Btime[0],format='(A11,F14.7)'

;	  print,'check Btime ',Btime[0],Btime[1],Btime[1]-Btime[0],$
;	    format='(A14,3F14.7)'
;	  print,'check Tseries ',Tseries[0],Tseries[1],Btime[0],$
;	    format='(A14,3F14.7)'
;	  print,'check diff ',Btime[0]-Tseries[0],(Btime[0]-Tseries[0])*8.64D04
;	  print,'check Btimepl ',Btimepl[0],Btimepl[1],Btimepl[sizeB-1]
;	  print,'check Btimepl[20,40,60] ',btimepl[20],btimepl[40],Btimepl[60]
;	  print,'check Btimepl[19,20,21] ',Btimepl[19],Btimepl[20],Btimepl[21]
;	  print,'check Bevent[40] ',Bevent[0:2,40]
	  thismax = max(Bevent[0:2,0:256-1],min=thismin)
;	  print,'whole event, min,max ',thismin,thismax
;
	  plot,Btimepl[0:sizeB-1],Bevent[1,0:sizeB-1],charsize=2.2,$
	  Ytitle='B (nT)',xstyle=17,yrange=[thismin,thismax],$
	  position=[.1, y0, .99, y1],xticklen=.1,/nodata
;
	  oplot,Btimepl,Bevent[0,0: NBload-1],color=red
	  oplot,Btimepl,Bevent[1,0: NBload-1],color=blue
	  oplot,Btimepl,Bevent[2,0: NBload-1],color=green
	  var = moment(Bevent[0,0: NBload-1], sdev=ctrms)
;	  print,'rms in counts ',ctrms
	  var = moment(Bevent[0,0: NBload-1], sdev=vrms)
	  prms = string(format='(f8.1,a4,f7.2,a3)',ctrms,' cts',$
		1000.*vrms,' nT')
	  xyouts, .70, .2*y0+.8*y1, prms,/normal
	  xline = [tsamp1,tsamp1]
	  yline = !Y.Crange
;	  print,'xline,yline ',xline,yline
	  oplot, xline,yline,color=red
	  xline = [tsamp2,tsamp2]
	  yline = !Y.Crange
;	  print,'xline,yline ',xline,yline
	  oplot, xline,yline,color=red
;
;
;	  xticklen=.05,xtickname=nulltick,/ynozero,color=magenta
;	  xticklen=.05,xtickname=nulltick,/ynozero,color=cyan
;
;
	timeS = systime(0)
	xyouts,0.,0.,scname+'APM_ninestrip_Mag_pub.ps, from APM_Analysis_Poynting'$
	+', version '+$
	version+',  output at '+ timeS,charsize=.7,/normal
; 	      stop	OK
;
;
if ihc ne 1 then begin
	device,/close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth6
;
;
ihc = 1			; here 1 is screen, do first
;
plotboth8: print,'at plotboth8, ihc = ',ihc
;
	!p.background=white
	!p.color=black
	!p.multi=[0,1,3]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_E,B,Poyn_figure.ps',xsize= 18.,ysize=20.,$
	yoffset=4.,/portrait,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;	plot the three E waveforms in one panel
	xspmin = float(floor(tsamp1))
	xspmax = float(floor(Tsamp2)) + 1.
	ymax0 = max(1000.*(All_Phys_RTNF[0,N1:N2]-varAxRTNF[0]),min=ymin0)
	ymax1 = max(1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),min=ymin1)
	ymax2 = max(1000.*(All_Phys_RTNF[2,N1:N2]-varCxRTNF[0]),min=ymin2)
	print,'maxs ',ymax0,ymax1,ymax2
	print,'mins ',ymin0,ymin1,ymin2
	YMAX = YMAX0>YMAX1>YMAX2
	ymin = ymin0<ymin1<ymin2
	print,'ymin,max ',N1,N2,ymin,ymax
	  plot,timeax[N1:N2],1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),$
	  xstyle=17,title=scname+' '+timesc,subtitle=$
	  'sampled '+string(tsamp1,format='(f4.1)')+' to'+string(tsamp2,format='(f4.1)')+' s'$
	  ,xrange=[xspmin,xspmax],yrange=[ymin,ymax],charsize=2.0,ytitle='mV',$
	    xtitle='sec',xticks=6,thick=5,/nodata
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),thick=3,color=red
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[1,N1:N2]-varBxRTNF[0]),thick=3,color=blue
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[2,N1:N2]-varCxRTNF[0]),thick=3,color=green
;
	var0 = moment((All_Phys_RTNF[0,N1:N2]-varAxRTNF[0]),sdev=vrms0)
	var1 = moment((All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),sdev=vrms1)
	var2 = moment((All_Phys_RTNF[2,N1:N2]-varCxRTNF[0]),sdev=vrms2)
	totrmsE = sqrt(vrms0^2 + vrms1^2 + vrms2^2)
	print,'totrmsE = ',totrmsE
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .04*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'X',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Y',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Z',color=green
;	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
;	  xyouts,xplace,yplace,'(c)'
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,' RTN'

;
	Thismax = max(Bevent[0,nsel:nsel2]-VarAxBsm[0],min=thismin)
	Thismax = max(Bevent[1,nsel:nsel2]-VarBxBsm[0],min=thismin1)>thismax
	thismin = thismin<thismin1
	Thismax = max(Bevent[2,0: NBload-1]-VarCxBsm[0],min=thismin2)>thismax
	thismin = thismin<thismin2
	totrms = 0.
	var0 = moment((Bevent[0,nsel:nsel2]-VarAxBsm[0]),sdev=vrms0)
	var1 = moment((Bevent[1,nsel:nsel2]-VarAxBsm[0]),sdev=vrms1)
	var2 = moment((Bevent[2,nsel:nsel2]-VarAxBsm[0]),sdev=vrms2)
	totrmsB = sqrt(vrms0^2 + vrms1^2 + vrms2^2)
	print,'totrmsB = ',totrmsB
	plot,BtimePL[ 0: NBload-1],BEVENT[1, 0: NBload-1]-VarBxBev[0],$
	yrange=[thismin,thismax],xrange=[xspmin,xspmax],xtitle='sec',$
	xstyle=17, ytitle='B (nT)',charsize=2.0,xticks=6,thick=5,/nodata
	oplot,BtimePL[nsel:nsel2],BEVENT[0,nsel:nsel2]-VarAxBsm[0],color=red,thick=3
	oplot,BtimePL[nsel:nsel2],BEVENT[1,nsel:nsel2]-VarBxBsm[0],color=blue,thick=3
	oplot,BtimePL[nsel:nsel2],BEVENT[2,nsel:nsel2]-VarCxBsm[0],color=green,thick=3
;
	xplace = .3*!x.crange[0] + .7*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'RMS '+string(totrmsB,format='(F7.4)')+' nT'
;
;    Poynting[*,nPoY-nsel+1] = PoYAvr[*]
	ymax = max(Poynting[*,*],min=ymin)
	print,'poynting 0 ',max(poynting[0,*],min=ymin)
;	ymax1 = ymax
;	print,'poynting 1 ',max(poynting[1,*],min=ymin1)>ymax1
;	ymin = ymin<ymin1
;	print,'poynting 2 ',max(poynting[2,*],min=ymin1)
;	ymin=ymin<ymin1
    plot,BtimePL[nsel:nsel2],Poynting[0,*],xstyle=17,$
    yrange=[ymin,ymax],xrange=[xspmin,xspmax],xtitle='sec',$
	ytitle='Poynting V',charsize=2.,xticks=6,thick=5,/nodata
    oplot,BtimePL[nsel:nsel2],Poynting[0,*],color=red,thick=3
    oplot,BtimePL[nsel:nsel2],Poynting[1,*],color=blue,thick=3
    oplot,BtimePL[nsel:nsel2],Poynting[2,*],color=green,thick=3
	var = moment(poynting[0,*],sdev=vrms)
	varpl = fltarr(nsel2-nsel+1)
	varpl = [var[0],var[0]]
	oplot,[BtimePL[nsel],BtimePL[nsel2]],varpl,color=red
	xplace = .3*!x.crange[0] + .7*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'average Pr '+string(var[0],format='(F8.4)')
;

print,'xrange ',xspmin,xspmax
print,'poy lims ',BtimePL[0],Btimepl[ NBload-1]
;print,'last POy ',poynting[*, NBload-1]
print,'size ',sizeB,n_elements(Poynting)
	timeS = systime(0)
	xyouts, 0., -.03,scname+' APM_E,B,Poyn_figure.ps, version '$
	+version+ ',  from APM_Analysis_Poynting,  output at '$
	+ timeS,charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth8
;
ihc = 1			; here 1 is screen, do first
;
plotboth9: print,'at plotboth9, ihc = ',ihc
;
;	This figure belongs in the first paper, but it got to complicated to
;		put it into APM_Analysis,  then anyway I decided on only one paper
;
	!p.background=white
	!p.color=black
	!p.multi=[0,1,2]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_E,B_figure.ps',xsize= 18.,ysize=22.,/color,yoffset=4.
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;	plot the three E waveforms in one panel
	xspmin = float(floor(tsamp1))
	xspmax = float(floor(Tsamp2)) + 1.
	ymax0 = max(1000.*(All_Phys_RTNF[0,N1:N2]-varAxRTNF[0]),min=ymin0)
	ymax1 = max(1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),min=ymin1)
	ymax2 = max(1000.*(All_Phys_RTNF[2,N1:N2]-varCxRTNF[0]),min=ymin2)
	print,'maxs ',ymax0,ymax1,ymax2
	print,'mins ',ymin0,ymin1,ymin2
	YMAX = YMAX0>YMAX1>YMAX2
	ymin = ymin0<ymin1<ymin2
	print,'ymin,max ',N1,N2,ymin,ymax
;	toptitle9 = 'samples '+string(N1f,format='(I2)'+ 'to '+string(N2f,format='(I2)"
	  plot,timeax[N1:N2],1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),$
	  xrange=[xspmin,xspmax],yrange=[ymin,ymax],charsize=1.8,ytitle='E mV/m',$
	  xtitle='sec',title=timesc,/nodata
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),color=red
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[1,N1:N2]-varBxRTNF[0]),color=blue
	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[2,N1:N2]-varCxRTNF[0]),color=green
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .04*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'X',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Y',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Z',color=green
  	  xplace = xplace + xinc
	  xyouts,xplace,yplace,' RTN'
;
;	put rms
;
    varx = moment(1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),sdev=stdx)
    vary = moment(1000.*(All_phys_RTNF[1,N1:N2]-varAxRTNF[0]),sdev=stdy)
    varz = moment(1000.*(All_phys_RTNF[2,N1:N2]-varAxRTNF[0]),sdev=stdz)
    rmstot = sqrt(stdx^2 + stdy^2 + stdz^2)
    print,'E rms tot ',rmstot
    xplace = .3*!x.crange[0] + .7*!x.crange[1]
    xyouts,xplace,yplace,'RMS '+string(rmstot,format='(F6.1)') + ' mV/m'
	Thismax = max(Bevent[0:2,nsel:nsel2],min=thismin)
	print,'BFig min,max ',thismin,thismax
    plot,BtimePL[nsel:nsel2],Bevent[1,nsel:nsel2],$
    yrange=[thismin,thismax],xrange=[xspmin,xspmax],xtitle='sec',$
    title='APM samples '+string(N1f,format='(I4)')+' to'+string(N2f,format='(I4)')$
	,ytitle='B (nT)',charsize=1.8,/nodata
    oplot,BtimePL[nsel:nsel2],Bevent[0,nsel:nsel2],color=red
    oplot,BtimePL[nsel:nsel2],Bevent[1,nsel:nsel2],color=blue
    oplot,BtimePL[nsel:nsel2],Bevent[2,nsel:nsel2],color=green
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .04*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'X',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Y',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Z',color=green
  	  xplace = xplace + xinc
	  xyouts,xplace,yplace,' RTN'
;
;
;	put rms
;
    varBx = moment((Bevent[0,nsel:nsel2]-varBx),sdev=stdBx)
    varBy = moment(Bevent[1,nsel:nsel2],sdev=stdBy)
    varBz = moment(Bevent[2,nsel:nsel2],sdev=stdBz)
    rmstot = sqrt(stdBx^2 + stdBy^2 + stdBz^2)
    print,'B rms tot ',rmstot
    xplace = .3*!x.crange[0] + .7*!x.crange[1]
    yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'RMS '+string(rmstot,format='(F6.2)') + ' nT'
   ;
	  print,'time chk ',BtimePL[0],BtimePL[ NBload-1],BtimePL[sizeB-1]
;
	timeS = systime(0)
	xyouts, 0., -.03,scname+$
	' APM_E,B_figure.ps, from APM_Analysis_Poynting,  version '+version+$
	 ',  output at '+ timeS,charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth9
;
;
ihc = 1			; here 1 is screen, do first
;
plotboth10: print,'at plotboth10, ihc = ',ihc
;
;	This figure belongs in the first paper, but it got to complicated to
;		put it into APM_Analysis
;
	!p.background=white
	!p.color=black
	!p.multi=[0,2,2]
;
if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_E,B_hodo_figure.ps',xsize= 18.,ysize=15.,/color,yoffset=4.
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;
;	plot the three E waveforms in one panel
;	lower left, I want upper right
	xspmin = float(floor(tsamp1))
	xspmax = float(floor(Tsamp2)) + 1.
	ymax0 = max(1000.*(All_Phys_RTNF[0,N1:N2]-varAxRTNF[0]),min=ymin0)
	ymax1 = max(1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),min=ymin1)
	ymax2 = max(1000.*(All_Phys_RTNF[2,N1:N2]-varCxRTNF[0]),min=ymin2)
	print,'maxs ',ymax0,ymax1,ymax2
	print,'mins ',ymin0,ymin1,ymin2
	YMAX = YMAX0>YMAX1>YMAX2
	ymin = ymin0<ymin1<ymin2
	print,'ymin,max ',N1,N2,ymin,ymax
	print,'N1,N2,ymin,max ',N1,N2,ymin,ymax
;	toptitle9 = 'samples '+string(N1,format='(I2)'+ 'to '+string(N2,format='(I2)"
	  plot,timeax[N1:N2],1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),title=' ',$
	  subtitle='sampled'$
	  +string(tsamp1,format='(f4.1)')+' to'+string(tsamp2,format='(f4.1)')+' s'$
	  ,xrange=[xspmin,xspmax],yrange=[ymin,ymax],charsize=1.4,ytitle='E mV/m',$
	    xtitle='sec',/nodata
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),color=red
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[1,N1:N2]-varBxRTNF[0]),color=blue
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[2,N1:N2]-varCxRTNF[0]),color=green
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .05*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'X',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Y',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Z',color=green
  	  xplace = xplace + xinc
	  xyouts,xplace,yplace,' RTN'
;
;	second hodo
;
	  plot,All_Phys_BF[0,N1:N2],All_Phys_BF[1,N1:N2],xrange=PFLIM0,$
	  yrange=PFLIM1,xtitle='X (mV/m) Bsys',ytitle='Y (mV/m) Bsys',$
	  xticks=2,charsize=1.4
	  oplot,All_Phys_BF[0,N1:N1+nstsp],All_Phys_BF[1,n1:N1+nstsp],color=green,thick=2
	  oplot,All_Phys_BF[0,N2-nstsp:N2],$
	     ALL_Phys_BF[1,N2-nstsp:N2],color=red,thick=2
	  Ycenter = [varBxF[0],varBxF[0]]
	  Xcenter = [varAxF[0],varAxF[0]]
	  oplot,Xcenter,Ycenter,psym=4,color=red
  	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(rotxy,format='(F6.1)')
  	  xplace = .2*!x.crange[0] + .8*!x.crange[1]
	  yplace = .9*!y.crange[0] + .1*!y.crange[1]
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
;	  xyouts,xplace,yplace,'(d)'
;
;	lower right, I want B
	Thismax = max((Bevent[0,nsel:nsel2]-varBx[0]),min=thismin)
	Thismaxy = max((Bevent[1,nsel:nsel2]-varBy[0]),min=thisminy)
	thismax = thismax>thismaxy
	thismin = thismin<thisminy
	Thismaxz = max((Bevent[2,nsel:nsel2]-varBz[0]),min=thisminz)
	thismax = thismax>thismaxz
	thismin = thismin<thisminz
	print,'yz ',thismaxy,thismaxz,thisminy,thisminz
	print,'BFig min,max ',thismin,thismax
    plot,BtimePL[nsel:nsel2],Poynting[0,*],$
    yrange=[thismin,thismax],xrange=[xspmin,xspmax],xtitle='sec',$
	ytitle='B (nT)',xticks=4,charsize=1.4,/nodata
    oplot,BtimePL[nsel:nsel2],Bevent[0,nsel:nsel2]-varBx[0],color=red
    oplot,BtimePL[nsel:nsel2],Bevent[1,nsel:nsel2]-varBy[0],color=blue
    oplot,BtimePL[nsel:nsel2],Bevent[2,nsel:nsel2]-varBz[0],color=green
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .05*(!x.crange[1] - !x.crange[0])
	  xyouts,xplace,yplace,'X',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Y',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'Z',color=green
  	  xplace = xplace + xinc
	  xyouts,xplace,yplace,' RTN'
;
	  print,'time chk ',BtimePL[0],BtimePL[ NBload-1],BtimePL[sizeB-1]
;
	timeS = systime(0)
	xyouts, 0., -.03,scname+$
	' APM_E,B,hodo_figure.ps, from APM_Analysis_Poynting, version '+version+$
	 ',  output at '+ timeS,charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
ihc = 1-ihc
;
if ihc ne 1 then goto,plotboth10
;
;
ihc = 1			; here 1 is screen
;
plotboth11: print,'at plotboth11, ihc = ',ihc
	!p.background=white
	!p.color=black
	!p.multi=[0,3,2]
;
;	Six Panels
;
 print,' B 2nd ',ALL_PHYS_B[*,1024]
 print,'at plotboth11 N1,N2 ',N1,N2
;	stop
;
;	ihc = 1				; necessary to avoid close 15 lines above
;
;	Panel 1
;

if ihc ne 1 then begin
	mydevice = !D.NAME
	set_plot, 'ps',/copy
	device,filename='APM_Hodogram_new.ps',xsize= 23., ysize=14.2,$
	yoffset=25.5,/landscape,/color
endif
;
;	PLOT STUFFF IN HERE, FOR EXAMPLE:
;
;
     xYmax = max(Bphys[N1:N2],min=xYmin)			; for plot range
     print,'Bphys min, max ',xYmin,xYmax
;
;
 	varAxRTNF = moment(All_Phys_RTNF[0,N1:N2])
 	varBxRTNF = moment(All_Phys_RTNF[1,N1:N2])
 	varCxRTNF = moment(All_Phys_RTNF[2,N1:N2])
;
;														Panel 1
;	plot the three waveforms in one panel
	xspmin = float(floor(tsamp1))
	xspmax = float(floor(Tsamp2)) + 1.
	  plot,timeax[N1:N2],1000.*(All_Phys_RTNF[1,N1:N2]-varBxRTNF[0]),$
	  title='samples '+string(N1,format='(I5)')+' to'+string(N2,format='(I5)')$
	  ,xrange=[xspmin,xspmax],charsize=2.0,ytitle='mV',$
	    xtitle='sec',/nodata
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),color=red
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[1,N1:N2]-varBxRTNF[0]),color=blue
	  	  oplot,timeax[N1:N2],1000.*(All_phys_RTNF[2,N1:N2]-varCxRTNF[0]),color=green
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .07*(!x.crange[1] - !x.crange[0])
	  yinc =  .07*(!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'R',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'T',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'N',color=green
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(a)'
	;
;	put rms
;
    varx = moment(1000.*(All_phys_RTNF[0,N1:N2]-varAxRTNF[0]),sdev=stdx)
    vary = moment(1000.*(All_phys_RTNF[1,N1:N2]-varAxRTNF[0]),sdev=stdy)
    varz = moment(1000.*(All_phys_RTNF[2,N1:N2]-varAxRTNF[0]),sdev=stdz)
    rmstot = sqrt(stdx^2 + stdy^2 + stdz^2)
    print,'E rms tot ',rmstot
    xplace = .9*!x.crange[0] + .1*!x.crange[1]
    yplace = yplace - yinc
        yplace = yplace - yinc
    xyouts,xplace,yplace,'RMS '+string(rmstot,format='(F4.1)') + ' mV/m',$
      charsize=.8
;
;														End of Panel 1
;														Panel 2, Poynting,V
;    Poynting[*,nPoY-nsel+1] = PoYAvr[*]
	ymax = max(Poynting[*,*],min=ymin)
;	print,'poynting 0 ',max(poynting[0,*],min=ymin)
;	ymax1 = ymax
;	print,'poynting 1 ',max(poynting[1,*],min=ymin1)>ymax1
;	ymin = ymin<ymin1
;	print,'poynting 2 ',max(poynting[2,*],min=ymin1)
;	ymin=ymin<ymin1
    plot,BtimePL[nsel:nsel2],Poynting[0,*],xstyle=17,$
    yrange=[ymin,ymax],xrange=[xspmin,xspmax],xtitle='sec',$
	ytitle='Poynting V',title='        '+scname+' '+timesc,subtitle=' ',$
	charsize=2.,xticks=4,thick=5,/nodata
    oplot,BtimePL[nsel:nsel2],Poynting[0,*],color=red,thick=2
    oplot,BtimePL[nsel:nsel2],Poynting[1,*],color=blue,thick=2
    oplot,BtimePL[nsel:nsel2],Poynting[2,*],color=green,thick=2
	var = moment(poynting[0,*],sdev=vrms)
	varpl = fltarr(nsel2-nsel+1)
	varpl = [var[0],var[0]]
	oplot,[BtimePL[nsel],BtimePL[nsel2]],varpl,color=red
	xplace = .7*!x.crange[0] + .3*!x.crange[1]
	yplace = .1*!y.crange[0] + .9*!y.crange[1]
	xyouts,xplace,yplace,'avr Pr '+string(var[0],format='(F8.4)')
;
;														End Panel 2 Poynting V
;														Start Panel 3
;	first hodo, upper right
;
	  plot,All_Phys_BF[0,N1:N2],All_phys_BF[2,N1:N2],xtitle='X (mV/m) Bsys',$
	  ytitle='E par B (mV/m)',charsize=2.,xrange=PLIM0,yrange=PLIM2,$
	  title = ' ',xticks=3,xstyle=17,/ynozero,/nodata
	   oplot,All_Phys_BF[0,N1:N2],All_phys_BF[2,N1:N2],thick=2
	   oplot,All_phys_BF[0,N1:N1+nstsp],All_phys_BF[2,N1:N1+nstsp],$
	   color=green,thick=2
	   oplot,All_phys_BF[0,N2-nstsp:N2],$
	     All_phys_BF[2,N2-nstsp:N2],color=red,thick=2
	  Xcenter = [varAxF[0],varAxF[0]]
	  Ycenter = [varCxF[0],varCxF[0]]
;	  print,'Ycenter ',Ycenter
;	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,color=red
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(Angsum/findex,format='(F6.1)')
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(e)'
;														End panel 3 hodo
;														Start panel 4, B
;
;	plot B instead of spectrum
;
;	for the files from CDAWeb
;		already loaded,  now select corresponding samples
	print,' mag,APM time ',scetmag,scet,format='(A18,2F14.6)'
	print,'N1,N2 ',N1,N2
	print,'tsamp1,2 ',tsamp1,tsamp2
;
;		find mag samples corresponding to selection
;
;	Tsamp1 is in sec = 1./86400 of aday
;		Tseries is in days
;		Btime is in days
;	1 msec = 1./8.64D07 of a day
    print,'Tseries[0],Btime[0] ',Tseries[0],Btime[0]
	scetmag = 0.d00
    scett = Tseries[0] + Tsamp1/8.6400D04			; day
    nsel = 0
    print,' '
	while SCETMAG LT scett do begin
	  scetmag = Btime[nsel]
;	  print,'going to get Btime, scetmag,scett ',scetmag,scett
;	  print,'tchk ',hr,mmin,sec1,usec
;	  print,'mag chk ',scetmag-scet,mmin,denst,tideg,Bnowt
	  nsel = nsel+1
	ENDwhile
	print,'at 1471,search Btime select, nsel= ',nsel
	print,'scet,mag,diff ',scett,scetmag,scetmag-scett,format='(A16,3F16.8)'
;	load B data, assumed to be 8 samples per sec so 256 samples
	Nsel2 = Nsel + (N2-N1)/8
;
 print,'nsel,Btime[0],Btime[nsel] ',nsel,Btime[0],Btime[nsel]
;  N1 and N2 are the first and last APM sample numbers for E selection
;  nsel and nsel2 are same for B samples in the selection
 print,'going to plot times ',Btimepl[nsel],Btimepl[nsel2]
 for npl = nsel-1,nsel2 do begin
   Btimepl[npl] = (Btime[npl]-Tseries[0])*8.64D04 			; days to sec
;   print,'pl ',npl,Bevent[0:2,npl],format='(A10,I5,3F14.8)'
 endfor
 print,'BtinePL ',Btimepl[nsel-1],Btimepl[nsel2]
 print,'got to 1483,nsel,nsel2 ', nsel,nsel2
;	NBload is the number of B samples in the full event, 32 sec
;   BtimePL
;
 	varAxBev = moment(Bevent[0, 0:NBload-1])
 	varBxBev = moment(Bevent[1, 0:NBload-1])
 	varCxBev = moment(Bevent[2, 0:NBload-1])
 	varAxBsm = moment(Bevent[0,nsel:nsel2])
 	varBxBsm = moment(Bevent[1,nsel:nsel2])
 	varCxBsm = moment(Bevent[2,nsel:nsel2])
	Thismax = max(Bevent[0,Nsel:nsel2]-VarAxBsm[0],min=thismin)
	thismin0 = thismin
	Thismax = max(Bevent[1,Nsel:nsel2]-VarBxBsm[0],min=thismin1)>thismin
	thismin = thismin<thismin1
	Thismax = max(Bevent[2,Nsel:nsel2]-VarCxBsm[0],min=thismin2)>thismin
	thismin = thismin2<thismin
	print,'thismin,max ',thismin,thismax
	plot,BtimePL[nsel:nsel2],BEVENT[1,nsel:nsel2]-VarBxBev[0],$
	yrange=[thismin,thismax],xrange=[xspmin,xspmax],xtitle='sec',$
	ytitle='B (nT)',charsize=2,/nodata
	oplot,BtimePL[nsel:nsel2],BEVENT[0,nsel:nsel2]-VarAxBsm[0],color=red
	oplot,BtimePL[nsel:nsel2],BEVENT[1,nsel:nsel2]-VarBxBsm[0],color=blue
	oplot,BtimePL[nsel:nsel2],BEVENT[2,nsel:nsel2]-VarCxBsm[0],color=green
	xplace = .15*!x.crange[0] + .85*!x.crange[1]
	yplace = .15*!y.crange[0] + .85*!y.crange[1]
	xyouts,xplace,yplace,'(b)'
;
;		put direction designation
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xinc   = .07*(!x.crange[1] - !x.crange[0])
	  yinc =  .07*(!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'R',color=red
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'T',color=blue
	  xplace = xplace + xinc
	  xyouts,xplace,yplace,'N',color=green
;
;	put rms
;
;    varBx = moment((Bevent[0,nsel:nsel2]-varBx),sdev=stdBx)
;    varBy = moment(Bevent[1,nsel:nsel2],sdev=stdBy)
;    varBz = moment(Bevent[2,nsel:nsel2],sdev=stdBz)
    rmstot = sqrt(stdBx^2 + stdBy^2 + stdBz^2)
    print,'B rms tot ',rmstot
    xplace = .85*!x.crange[0] + .15*!x.crange[1]
    yplace = yplace - yinc
    yplace = yplace - yinc
    xyouts,xplace,yplace,'RMS '+string(rmstot,format='(F4.1)') + ' nT',$
    charsize=.8
;
;
;														end panel 4, B
;														Start Panel 5. Spectrum
 print,'start panel 5 at 2820 '
;
	  nfilt = N2-N1+1
      fltfreq = fltarr(nfilt/2)
      fltfundfr = 1./(delt*nfilt)
	  fltfreq = fltfundfr*indgen(nfilt/2)
;
;	  determine observed frequency
;

;      fltfreq = fltarr(nfilt/2)
;
; print,'Awave, N1,N2 ',N1f,N2f
; print,'nfilt ',nfilt
; print,'Asize ',n_elements(Awave)
; print,Awave[*]
;print,'Awave[1:20]',4096.*abs(Awave[1:20])^2
;print,'Bwave[1:20]',4096.*abs(Bwave[1:20])^2
;print,'Cwave[1:20]',4096.*abs(Cwave[1:20])^2
;print,'xyratio ',xyratio
;print,'Ewave[1:20]',4096.*abs(Ewave[1:20])^2
	  Xmax = max(4096.*abs(Awave[1:20])^2,indX)
  	  Ymax = max(4096.*abs(Bwave[1:20])^2,indY)
  	  Zmax = max(4096.*abs(Cwave[1:20])^2,indZ)
  	  Emax = max(4096.*abs(Ewave[1:20])^2,indE)
  	  findex = indX
 print,'Xmax,Y,Z',xmax,Ymax,Zmax,Emax
;
;Mwave[0:nfilt-1] = Awave[0:nfilt-1]
;  	  TpowerM = moment(abs(Aphys[N1f:N2f]))
;  	  Tpower  = (N2f-N1f+1)*TpowerM[1]
  	  if (Ymax gt Xmax) and (IndY ge Indx) then begin
  	  	findex = indY
 	    Mwave = fft(All_phys_RTN[1,N1f:N2f])
 	  endif
  	  if (Zmax gt Xmax) and (Zmax gt Ymax) and (IndZ ge findex) then begin
  	    findex = indZ
  	    Mwave = fft(Cphys[N1f:N2f])
  	  endif
;  	  print,'Mwave[*]',4096.*abs(Mwave[*])^2
  	  TpowerM = moment(abs(Mwave[1:(N2f-N1f-1)]))
  	  Tpower = (N2f-N1f+1)*TpowerM[1]
;  	  print,'Tpower from moment ',Tpower,TpowerM[1]
  	  findex = findex+1
  	  TpowerC = 0.					;
  	  for ncl = 1,nfilt-1 do begin
  	    TpowerC = TpowerC + abs(Mwave[ncl])^2
  	  endfor
;  	  print,'endchk ',abs(Mwave[nfilt-1])^2
;  	  print,'TpowerC ',TpowerC
; 	  print,'N1,N2 ',N1f,N2f
  	  print,'i ',indx,indy,indz,inde
  	  print, 'choice',findex
;		calculate signal power and Bckgpower
	  Sigpower = (Abs(Mwave[findex-1]))^2 + (Abs(Mwave[findex]))^2 + $
	     (Abs(Mwave[findex+1]))^2
;	  Sigpower = 2.*Sigpower
	  Bckgpower = TpowerC - Sigpower
	  print,'Total',TpowerC
	  print,'Signl',Sigpower
	  print,'Bckg ',Bckgpower
	  Bckgpower = Bckgpower>0.
  	  print,'ctr_freq ',ctr_freq
	  plot,fltfreq[1:nfilt/2-1],4096.*abs(Awave[1:nfilt/2-1])^2,xtitle='freq Hz',$
	  title='frequency spectrum ',xrange=[(1./32.),32.],xstyle=17,$
	  charsize=1.8,thick=2,/xlog
	  oplot,fltfreq[1:nfilt/2-1],4096.*abs(Awave[1:nfilt/2-1])^2,psym=4,symsize=.5
	  print,'nfilt ',nfilt
	  print,'Awave ',abs(Awave[1]),abs(Awave[2]),abs(Awave[nfilt/2-1])
;	stop
	  xplace = .7*!x.crange[0] + .3*!x.crange[1]
	  yplace = .2*!y.crange[0] + .8*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,10.^(xplace),yplace,'freq ='+string(ctr_freq,format='(F7.3)')+' Hz'
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xplace = exp(xplace)
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xyouts,xplace,yplace,'(d)'
;  	  stop
;												End Panel 5
;
;	calculate powers in peak and rest
;	  stop
;
;			Calculate rotation for X,Y, unfiltered data
;			VarAx is moment for unfiltered data, N1 to N2,
 	varAxB = moment(All_Phys_B[0,N1:N2])
 	varBxB = moment(All_Phys_B[1,N1:N2])
 	varCxB = moment(All_Phys_B[2,N1:N2])
;
;		for plot limits
	  xBmax0 = max(All_Phys_B[0,N1:N2],min=XBmin0)
	  xBmax1 = max(All_Phys_B[1,N1:N2],min=XBmin1)
	  xBmax2 = max(All_Phys_B[2,N1:N2],min=XBmin2)
	  print,'Xbmin,XBmax ',XBmin0,XBmax0,XBmin1,XBmax1,XBmin2,XBmax2
;
;	make square plots, all with same range, Hrange
;
  Hrange = (XBmax0-XBmin0)>(XBmax1-XBmin1)>(XBmax2-XBmin2)
  Dhr    = .5*Hrange
  Plim0 = [.5*(XBmax0+xBmin0)-Dhr,.5*(XBmax0+xBmin0)+Dhr]
  Plim1 = [.5*(XBmax1+xBmin1)-Dhr,.5*(XBmax1+xBmin1)+Dhr]
  Plim2 = [.5*(XBmax2+xBmin2)-Dhr,.5*(XBmax2+xBmin2)+Dhr]
;
	  plot,All_Phys_BF[0,N1:N2],All_Phys_BF[1,N1:N2],xtitle='X Bsys',ytitle='Y Bsys',$
	    charsize=1.8,xrange=Plim0,yrange=Plim1,xticks=3,xstyle=17,/ynozero,/nodata
	  oplot,All_Phys_BF[0,N1:N2],All_Phys_BF[1,N1:N2],thick=2
	  oplot,All_phys_BF[0,N1:N1+nstsp],All_phys_BF[1,N1:N1+nstsp],$
	    color=green,thick=2
	  oplot,All_phys_BF[0,N2-nstsp:N2],$
	     All_phys_BF[1,N2-nstsp:N2],color=red,thick=2

	  Ycenter = [varBxB[0],varBxB[0]]
	  Xcenter = [varCxB[0],varCxB[0]]
;	  print,'Ycenter ',Ycenter
;	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,thick=2,color=red
;
	  xplace = .9*!x.crange[0] + .1*!x.crange[1]
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  yinc   = (!y.crange[1] - !y.crange[0])
	  xyouts,xplace,yplace,'rot ='+string(rotxy,format='(F4.1)')
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(f)'
;
  	  xplace = .2*!x.crange[0] + .8*!x.crange[1]
	  yplace = .9*!y.crange[0] + .1*!y.crange[1]
	  Binout = Binoutsv
	  if Binout eq 1 then xyouts,xplace,yplace,'B in'
	  if Binout eq 0 then xyouts,xplace,yplace,'B out'
;
;
;check plot
;
;	   band pass filter
;
;	The Nyquist frequency is .5/delt
;   peak freq in terms of nyquist
;   ctr/nyquist = Ctr_freq/fnyquist
;    filter from 1/2 center freq to 5*center freq
;    or .5 * findex/nfilt to 5 * findex/nfilt
     flow = .5*ctr_freq/fnyquist
	 fhigh = 5.*ctr_freq/fnyquist
 Ifilt = 1
 if Ifilt ne 0 then begin
   order = 8
   bpfilter = Digital_Filter(flow, fhigh, 40., order)
; the IDL routine does not normalize????.  So do it
; yes it does with ??
   bpsize = n_elements(bpfilter)
;   print,'bpsize ',bpsize
   Sumbp = 0.
   for nbp = 0,2*order do begin
     sumbp = sumbp + bpfilter[nbp]
   endfor
; print,'sumbp ',sumbp
   bpfilter[*] = bpfilter[*]/sumbp
	print,'f ',ctr_freq,flow,fhigh
;   stop
;
;	example   AFphys = Convol(Aphys,hpfilter)
;	print,'bpfilter ',bpfilter
   All_phys_BF = fltarr(3,size)
   Atemp1 = fltarr(size)
   Atemp2 = fltarr(size)
;
;	make All_Phys_BF filtered
;
   Atemp2[*] = All_Phys_B[0,0:size-1]
   print,'_B 4th',All_Phys_B[*,1024]
   Atemp1[0:size-1] =  Convol(Atemp2,bpfilter)
   print,'BP filter ',bpfilter
   All_Phys_BF[0,0:size-1] = Atemp1[0:size-1]
;   print,'At 1313 ',ATemp1[*]
;
   ;stop
   print,'_BF ',All_Phys_BF[0,0],All_Phys_BF[0,1],All_Phys_BF[0,size-1]
   Atemp2[*] = All_Phys_B[1,0:size-1]
   ATemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_BF[1,0:size-1] = Atemp1[0:size-1]
;
   Atemp2[*] = All_Phys_B[2,0:size-1]
   Atemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_BF[2,0:size-1] = Atemp1[0:size-1]
; endif
;
;	Make RTN filtered
;
 print,'check for panel s '
 print,'all_phys_RTNF ',All_Phys_RTNF[*,1]
 print,'rot ',rotxy
   All_phys_RTNF = fltarr(3,size)
;
   Atemp2[*] = All_Phys_RTN[0,0:size-1]
;   print,'_B 4th',All_Phys_B[*,1024]
   Atemp1[0:size-1] =  Convol(Atemp2,bpfilter)
   All_Phys_RTNF[0,0:size-1] = Atemp1[0:size-1]
;   print,'At ',ATemp1[*]
;
   ;stop
;   print,'_BF ',All_Phys_BF[0,0],All_Phys_BF[0,1],All_Phys_BF[0,size-1]
   Atemp2[*] = All_Phys_RTN[1,0:size-1]
   ATemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_RTNF[1,0:size-1] = Atemp1[0:size-1]
;
   Atemp2[*] = All_Phys_RTN[2,0:size-1]
   Atemp1[0:size-1] = Convol(Atemp2,bpfilter)
   All_Phys_RTNF[2,0:size-1] = Atemp1[0:size-1]
;stop

;			Calculate rotation for Z,X, filtered data
      Del_Ang 	= 0.
      Angsum	= 0.
;      print,'chk ',VarAx[0],varBx[0],varCx[0]
;		it seems than IDL atan(y,x) is arc tan (y/x)
      Angsv = atan((All_Phys_BF[2,N1]-varCxF[0]),(All_phys_BF[0,N1]-varAxF[0]))
      Ysv = (All_Phys_BF[2,N1]-VarBxF[0])
      for n = order+1,nfilt+order do begin
        Angxy = atan((All_PHys_BF[2,n]-varCxF[0]),(ALL_phys_BF[0,n]-varAxF[0]))
        Del_ang = Angxy - Angsv
        ynew = (All_Phys_BF[2,n]-varCxF[0])
;        if (All_phys_BF[1,n]-VarBxF[0])*(All_phys_BF[0,n-1]-VarAxF[0]) gt 0. then $
		 if (ysv*ynew gt 0.) and ((All_Phys_BF[0,n]-varAxF[0]) gt 0.) then $
        	Angsum = Angsum+Del_Ang
;	    print,'Dx,Angxy2Del,Sv,tot ',n,Aphys[n]-VarAx[0],Bphys[n]-VarBx[0],Angxy,$
;	    Del_ang,Angsv,angsum,format='(A18,I5,6F9.5)'
	    Angsv = Angxy
	    ysv = ynew
      endfor
;
;      print,'ZY rotation ',Angsum/findex, ' radians'
	  rotzy =  Angsum/findex
	  Xcenter = [varAxF[0],varAxF[0]]
	  Ycenter = [varBxF[0],varBxF[0]]
	  print,'Ycenter ',Ycenter
	  print,Xcenter
	  oplot,Xcenter,Ycenter,psym=4,thick=2,color=red
	  yplace = .1*!y.crange[0] + .9*!y.crange[1]
	  xplace = .15*!x.crange[0] + .85*!x.crange[1]
	  xyouts,xplace,yplace,'(f)'
;	*********** to fill squares
;
	endif
; *****************
;
	timeS = systime(0)
	xyouts,0.,0.,scname+' APM_Hodogram_new.ps, from APM_Analysis_Poynting_sv, version '$
	+version+ ',  output at '+ timeS,charsize=.7,/normal
;
if ihc ne 1 then begin
	device, /close
	set_plot, mydevice
endif
;
;**************
; ihc = 1
; goto, plotboth4
;**************
print,' end plotboth11, ihc= ',ihc
ihc = 1-ihc
;
;if ihcA eq 1 then read,dispose
if ihc ne 1 then goto,plotboth11
;
 Eenergy = 0.
 Benergy = 0.
 Duration = (Tseries[N2]-Tseries[N1])
 for nen = N1,N2 do begin
   Eenergy = Eenergy + (All_phys_RTNF[0,nen]-varAxRTNF[0])^2
   Eenergy = Eenergy + (All_phys_RTNF[1,nen]-varBxRTNF[0])^2
   Eenergy = Eenergy + (All_phys_RTNF[2,nen]-varCxRTNF[0])^2
 endfor
 Eenergy = .5*eps0*Eenergy
 for nen = nsel,nsel2 do begin
   Benergy = Benergy + (BEVENT[0,nen]-VarAxBsm[0])^2
   Benergy = Benergy + (BEVENT[1,nen]-VarBxBsm[0])^2
   Benergy = Benergy + (BEVENT[2,nen]-VarCxBsm[0])^2
 endfor
; E is V/m but B is nT
 Benergy = .5*1.e-18*Benergy/mu0
 print,'E.B,energy ',Eenergy,Benergy
	err = TM_UR8_to_Ydoyh(Scet,yyyy,doy,hh)
	err = TM_UR8_to_YMD(Scet,yyyy,mmon,dd,hh,min,sec,msec)
	yyyymmddp = 10000*yyyy + 100*mmon + dd
;	print,'yyyy chk ',yyyymmddp
	hhmmssp = 10000*HH + 100*MIN + SEC
	SCETI4 = [yyyymmddp,hhmmssp]
;
 print,'min,denst ',nmin,denst
;
close,3
close,8
close,66
;close,4
stop
end
