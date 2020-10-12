;Make plots of the burst 1 data availability.

rbsp_efw_init
tplot_options,'title','rbsp_burst1_availability_plots.pro'
tplot_options,'xmargin',[20.,16.]
tplot_options,'ymargin',[3,9]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',2
tplot_options,'ythick',2
tplot_options,'labflag',-1


d0 = time_double('2012-09-01')
d1 = time_double('2019-11-01')
timespan,d0,d1-d0,/sec

kyoto_load_dst

;get B2 times
b2ta = rbsp_get_burst2_times_list('a')
b2tb = rbsp_get_burst2_times_list('b')



;get B1 times and rates from this routine
b1ta = rbsp_get_burst_times_rates_list('a')
b1tb = rbsp_get_burst_times_rates_list('b')

goo16a = where((b1ta.samplerate gt 15000) and (b1ta.samplerate lt 18000))
goo8a = where((b1ta.samplerate gt 6000) and (b1ta.samplerate lt 11000))
goo4a = where((b1ta.samplerate gt 3500) and (b1ta.samplerate lt 5500))
goo2a = where((b1ta.samplerate gt 1700) and (b1ta.samplerate lt 3000))
goo1a = where((b1ta.samplerate gt 900) and (b1ta.samplerate lt 1600))
goo5a = where((b1ta.samplerate gt 0) and (b1ta.samplerate lt 800))
goo16b = where((b1tb.samplerate gt 15000) and (b1tb.samplerate lt 18000))
goo8b = where((b1tb.samplerate gt 6000) and (b1tb.samplerate lt 11000))
goo4b = where((b1tb.samplerate gt 3500) and (b1tb.samplerate lt 5500))
goo2b = where((b1tb.samplerate gt 1700) and (b1tb.samplerate lt 3000))
goo1b = where((b1tb.samplerate gt 900) and (b1tb.samplerate lt 1600))
goo5b = where((b1tb.samplerate gt 0) and (b1tb.samplerate lt 800))


;Total # collections of each type
;help,goo16a,goo8a,goo4a,goo2a,goo1a,goo5a
;GOO16a           LONG      = Array[3655]
;GOO8a            LONG      = Array[260]
;GOO4a            LONG      = Array[472]
;GOO2a            LONG      = Array[72]
;GOO1a            LONG      = Array[125]
;GOO5a            LONG      = Array[687]

;Total # hours of each type
dur16a = total(b1ta.duration[goo16a])/3600.
dur8a = total(b1ta.duration[goo8a])/3600.
dur4a = total(b1ta.duration[goo4a])/3600.
dur2a = total(b1ta.duration[goo2a])/3600.
dur1a = total(b1ta.duration[goo1a])/3600.
dur5a = total(b1ta.duration[goo5a])/3600.
dur16b = total(b1tb.duration[goo16b])/3600.
dur8b = total(b1tb.duration[goo8b])/3600.
dur4b = total(b1tb.duration[goo4b])/3600.
dur2b = total(b1tb.duration[goo2b])/3600.
dur1b = total(b1tb.duration[goo1b])/3600.
dur5b = total(b1tb.duration[goo5b])/3600.
;help,dur16a,dur8a,dur4a,dur2a,dur1a,dur5a
;DUR16a           FLOAT     =       622.864
;DUR8a            FLOAT     =       36.8708
;DUR4a            FLOAT     =       354.897
;DUR2a            FLOAT     =       45.8794
;DUR1a            FLOAT     =       191.093
;DUR5a            FLOAT     =       1439.94
;print,dur16a+dur8a+dur4a+dur2a+dur1a+dur5a
;total hrs =       2691.54

;help,dur16b,dur8b,dur4b,dur2b,dur1b,dur5b
;DUR16B          FLOAT     =       501.490
;DUR8B           FLOAT     =       23.2925
;DUR4B           FLOAT     =       616.875
;DUR2B           FLOAT     =       43.3994
;DUR1B           FLOAT     =       531.998
;DUR5B           FLOAT     =       1052.67
;print,dur16b+dur8b+dur4b+dur2b+dur1b+dur5b
;total hrs =       2769.72




;-----------------------------------------------------------
;Create tplot variables
;-----------------------------------------------------------

;create a time base
cadence = 1.*60.
nelem = (d1 - d0)/cadence
times = d0 + dindgen(nelem)*cadence

b1a_flag = fltarr(nelem) & b1b_flag = fltarr(nelem)
b2a_flag = fltarr(nelem) & b2b_flag = fltarr(nelem)

for q=0,n_elements(b1ta.startb1)-1 do begin $
  goodtimes = where((times ge b1ta.startb1[q]) and (times le b1ta.endb1[q])) & $
  if goodtimes[0] ne -1 then b1a_flag[goodtimes] = b1ta.samplerate[q]
endfor
for q=0,n_elements(b1tb.startb1)-1 do begin $
  goodtimes = where((times ge b1tb.startb1[q]) and (times le b1tb.endb1[q])) & $
  if goodtimes[0] ne -1 then b1b_flag[goodtimes] = b1tb.samplerate[q]
endfor

;for q=0,n_elements(b2ta.startb2)-1 do begin $
;  goodtimes = where((times ge b2ta.startb2[q]) and (times le b2ta.endb2[q])) & $
;  if goodtimes[0] ne -1 then b2a_flag[goodtimes] = 16384.
;for q=0,n_elements(b2tb.startb2)-1 do begin $
;  goodtimes = where((times ge b2tb.startb2[q]) and (times le b2tb.endb2[q])) & $
;  if goodtimes[0] ne -1 then b2b_flag[goodtimes] = 16384.



store_data,'b1a',times,b1a_flag
store_data,'b1b',times,b1b_flag
store_data,'b2a',times,b2a_flag
store_data,'b2b',times,b2b_flag
ylim,'b1?',400,20000,1
options,'b1?','psym',-2

;create grid lines
store_data,'grid5',times,replicate(512.,nelem)
store_data,'grid1',times,replicate(1024.,nelem)
store_data,'grid2',times,replicate(2048.,nelem)
store_data,'grid4',times,replicate(4096.,nelem)
store_data,'grid8',times,replicate(8192.,nelem)
store_data,'grid16',times,replicate(16384.,nelem)

store_data,'b1acomb',data=['b1a','grid5','grid1','grid2','grid4','grid8','grid16']
store_data,'b1bcomb',data=['b1b','grid5','grid1','grid2','grid4','grid8','grid16']
options,'b1acomb','ytitle','B1 data/datarate!Ctimeline!CRBSPa!C(Samples/sec)'
options,'b1bcomb','ytitle','B1 data/datarate!Ctimeline!CRBSPb!C(Samples/sec)'

ylim,'b1?comb',400,20000,1
;tplot,'b1?comb'


;-----------------------------------------------------------
;Running integration of telemetered burst data
;***OVERPLOT VARIOUS CAMPAIGNS ON HERE
;-----------------------------------------------------------



totsa = fltarr(n_elements(times))  ;running total of # samples
totsb = fltarr(n_elements(times))  ;running total of # samples
totsrun = 0.
for i=0d,n_elements(times)-1 do begin $
  totsrun += b1a_flag[i]*cadence & $
  totsa[i] = totsrun
endfor
store_data,'running_samplesa',times,totsa

totsrun = 0.
for i=0d,n_elements(times)-1 do begin $
  totsrun += b1b_flag[i]*cadence & $
  totsb[i] = totsrun
endfor
store_data,'running_samplesb',times,totsb


;Running integration of B2 files
tots = 0.
totsaccuma = fltarr(n_elements(b2ta.duration))
totsaccumb = fltarr(n_elements(b2tb.duration))
for j=0.,n_elements(b2ta.duration)-1 do begin $
  tots += b2ta.duration[j]*16384. & $
  totsaccuma[j] = tots
endfor
store_data,'running_samplesb2a',b2ta.startb2,totsaccuma

tots = 0.
for j=0.,n_elements(b2tb.duration)-1 do begin $
  tots += b2tb.duration[j]*16384. & $
  totsaccumb[j] = tots
endfor
store_data,'running_samplesb2b',b2tb.startb2,totsaccumb


;----------------------------------
;@@@@UNDER CONSTRUCTION@@@@@@
;Conversion from Samples to MB
;Pre  2013-10-14 rate is 4096 Samples/block
;Post 2014-10-14 rate is 5461 Samples/block

;;2018-01-28 MSCB1 RBa total size is 2.7 GB
;2018-01-28/15:56:19      869.000      16384.0
;2018-01-28/16:10:49      926.000      16384.0
;2018-01-28/16:26:16      1573.00      16384.0
;2018-01-28/16:52:30      1280.00      16384.0
;2018-01-28/17:13:51      3425.00      16384.0

;durtmp = 8073 ;sec
;sizee = 2700. ;MB

;@16K have 0.334448 MB/sec for 3 channels of mscb1
;I would think that the Vb1 should be twice this b/c there are 6 channels.
;Unfortunately, the Vb1 files is only 3.5 GB in size. Compression??


;@@@@UNDER CONSTRUCTION@@@@@@
;----------------------------------



;Data rate is once/min.
deriv_data,'running_samplesa',nsmooth=1440.*30.,newname='running_samplesa_monthly_rate'
deriv_data,'running_samplesb',nsmooth=1440.*30.,newname='running_samplesb_monthly_rate'
deriv_data,'running_samplesb2a',nsmooth=1440.*30.,newname='running_samplesb2a_monthly_rate'
deriv_data,'running_samplesb2b',nsmooth=1440.*30.,newname='running_samplesb2b_monthly_rate'


;deriv_data,'running_samplesa',nsmooth=1440.*360.,newname='running_samplesa_yearly_rate'
;deriv_data,'running_samplesa',nsmooth=1440.*7.,newname='running_samplesa_weekly_rate'


;BARREL campaigns
;https://earthweb.ess.washington.edu/mccarthy/BARREL/SWD2015/welcome.html
b13 = ['2013-01-01','2013-01-30']
b14 = ['2013-12-27','2014-02-11']
b_kiruna1 = ['2015-08-10','2015-08-26']
b_kiruna2 = ['2016-08-13','2016-08-31']
b_kiruna3 = ['2018-06-25','2018-06-27']
b18 = ['2018-12-09','2019-02-21']
;FIREBIRD campaigns
fb1 = ['2015-02-01','2015-02-21']
fb2 = ['2015-03-21','2015-04-19']
fb3 = ['2015-05-16','2015-06-15']
fb4 = ['2015-07-03','2015-08-04']
fb5 = ['2015-08-08','2015-09-04']
fb6 = ['2015-11-15','2015-12-15']
fb7 = ['2016-01-15','2016-02-03']
fb8 = ['2016-05-20','2016-06-20']
fb9 = ['2016-08-12','2016-09-07']
fb10 = ['2016-12-21','2017-01-04']
fb11 = ['2017-05-01','2017-05-21']
fb12 = ['2017-07-01','2017-07-21']
fb13 = ['2017-11-19','2017-12-14']
fb14 = ['2018-02-27','2018-03-28']
fb15 = ['2018-04-20','2018-05-13']
fb16 = ['2018-06-25','2018-07-18']
fb17 = ['2018-07-31','2018-08-20']
fb18 = ['2018-09-21','2018-10-13']
fb19 = ['2018-12-16','2019-01-11']
fb20 = ['2019-01-24','2019-02-02']
fb21 = ['2019-03-16','2019-04-05']
fb22 = ['2019-05-05','2019-05-17']
fb23 = ['2019-07-04','2019-07-17']
fb24 = ['2019-09-10','2019-09-20']

;Schumann resonance campaign
schum = ['2019-03-17','2019-05-08']


;***NEED ALL 2018/2019 CAMPAIGNS

;store_data,'running_samples_monthly_rate_comb',data=['running_samplesa_monthly_rate','running_samplesb_monthly_rate']
;store_data,'running_samples_monthly_rateb2_comb',data=['running_samplesb2a_monthly_rate','running_samplesb2b_monthly_rate']
store_data,'running_samples_monthly_rate_comb',data=['running_samplesa_monthly_rate','running_samplesb_monthly_rate','running_samplesb2a_monthly_rate','running_samplesb2b_monthly_rate']
options,'running_samples_monthly_rate_comb','colors',[0,50,200,250]
options,'running_samples_monthly_rate_comb','ytitle','Avg monthly burst!Cdata rate!C(Samples/sec)!Cblack=B1a,blue=B1b!COrange=B2a,Red=B2b'

;options,'running_samples_monthly_rate_comb','colors',[50,250]
;options,'running_samples_monthly_rateb2_comb','colors',[50,250]
;store_data,'running_samples_comb',data=['running_samplesa','running_samplesb']
;store_data,'running_samplesb2_comb',data=['running_samplesb2a','running_samplesb2b']
store_data,'running_samples_comb',data=['running_samplesa','running_samplesb','running_samplesb2a','running_samplesb2b']
options,'running_samples_comb','ytitle','Accumulated data volume!C(Samples)!Cblack=B1a,blue=B1b!COrange=B2a,Red=B2b'
;options,'running_samples_comb','colors',[50,250]
options,'running_samples_comb','colors',[0,50,200,250]
options,'running_samplesb2_comb','colors',[50,250]


;timebar,[b13,b14,b18,b_kiruna1,b_kiruna2,b_kiruna3],linestyle=2
;timebar,[fb1,fb2,fb3,fb4,fb5,fb6,fb7,fb8,fb9,fb10,fb11,fb12,fb13,fb14,fb15,fb16,fb17,fb18,fb19,fb20,fb21,fb22,fb23,fb24],color=250,linestyle=2
;timebar,schum,color=200


store_data,'b13',time_double(b13),replicate(1.2,2) & ylim,'b13',0,2 & options,'b13','psym',0 & options,'b13','thick',4
store_data,'b14',time_double(b14),replicate(1.2,2) & ylim,'b14',0,2 & options,'b14','psym',0 & options,'b14','thick',4
store_data,'b_kiruna1',time_double(b_kiruna1),replicate(1.2,2) & ylim,'b_kiruna1',0,2 & options,'b_kiruna1','psym',0 & options,'b_kiruna1','thick',4
store_data,'b_kiruna2',time_double(b_kiruna2),replicate(1.2,2) & ylim,'b_kiruna2',0,2 & options,'b_kiruna2','psym',0 & options,'b_kiruna2','thick',4
store_data,'b_kiruna3',time_double(b_kiruna3),replicate(1.2,2) & ylim,'b_kiruna3',0,2 & options,'b_kiruna3','psym',0 & options,'b_kiruna3','thick',4
store_data,'fb1',time_double(fb1),replicate(1,2) & ylim,'fb1',0,2 & options,'fb1','psym',0 & options,'fb1','thick',4 & options,'fb1','colors',250
store_data,'fb2',time_double(fb2),replicate(1,2) & ylim,'fb2',0,2 & options,'fb2','psym',0 & options,'fb2','thick',4 & options,'fb2','colors',250
store_data,'fb3',time_double(fb3),replicate(1,2) & ylim,'fb3',0,2 & options,'fb3','psym',0 & options,'fb3','thick',4 & options,'fb3','colors',250
store_data,'fb4',time_double(fb4),replicate(1,2) & ylim,'fb4',0,2 & options,'fb4','psym',0 & options,'fb4','thick',4 & options,'fb4','colors',250
store_data,'fb5',time_double(fb5),replicate(1,2) & ylim,'fb5',0,2 & options,'fb5','psym',0 & options,'fb5','thick',4 & options,'fb5','colors',250
store_data,'fb6',time_double(fb6),replicate(1,2) & ylim,'fb6',0,2 & options,'fb6','psym',0 & options,'fb6','thick',4 & options,'fb6','colors',250
store_data,'fb7',time_double(fb7),replicate(1,2) & ylim,'fb7',0,2 & options,'fb7','psym',0 & options,'fb7','thick',4 & options,'fb7','colors',250
store_data,'fb8',time_double(fb8),replicate(1,2) & ylim,'fb8',0,2 & options,'fb8','psym',0 & options,'fb8','thick',4 & options,'fb8','colors',250
store_data,'fb9',time_double(fb9),replicate(1,2) & ylim,'fb9',0,2 & options,'fb9','psym',0 & options,'fb9','thick',4 & options,'fb9','colors',250
store_data,'fb10',time_double(fb10),replicate(1,2) & ylim,'fb10',0,2 & options,'fb10','psym',0 & options,'fb10','thick',4 & options,'fb10','colors',250
store_data,'fb11',time_double(fb11),replicate(1,2) & ylim,'fb11',0,2 & options,'fb11','psym',0 & options,'fb11','thick',4 & options,'fb11','colors',250
store_data,'fb12',time_double(fb12),replicate(1,2) & ylim,'fb12',0,2 & options,'fb12','psym',0 & options,'fb12','thick',4 & options,'fb12','colors',250
store_data,'fb13',time_double(fb13),replicate(1,2) & ylim,'fb13',0,2 & options,'fb13','psym',0 & options,'fb13','thick',4 & options,'fb13','colors',250
store_data,'fb14',time_double(fb14),replicate(1,2) & ylim,'fb14',0,2 & options,'fb14','psym',0 & options,'fb14','thick',4 & options,'fb14','colors',250
store_data,'fb15',time_double(fb15),replicate(1,2) & ylim,'fb15',0,2 & options,'fb15','psym',0 & options,'fb15','thick',4 & options,'fb15','colors',250
store_data,'fb16',time_double(fb16),replicate(1,2) & ylim,'fb16',0,2 & options,'fb16','psym',0 & options,'fb16','thick',4 & options,'fb16','colors',250
store_data,'fb17',time_double(fb17),replicate(1,2) & ylim,'fb17',0,2 & options,'fb17','psym',0 & options,'fb17','thick',4 & options,'fb17','colors',250
store_data,'fb18',time_double(fb18),replicate(1,2) & ylim,'fb18',0,2 & options,'fb18','psym',0 & options,'fb18','thick',4 & options,'fb18','colors',250
store_data,'fb19',time_double(fb19),replicate(1,2) & ylim,'fb19',0,2 & options,'fb19','psym',0 & options,'fb19','thick',4 & options,'fb19','colors',250
store_data,'fb20',time_double(fb20),replicate(1,2) & ylim,'fb20',0,2 & options,'fb20','psym',0 & options,'fb20','thick',4 & options,'fb20','colors',250
store_data,'fb21',time_double(fb21),replicate(1,2) & ylim,'fb21',0,2 & options,'fb21','psym',0 & options,'fb21','thick',4 & options,'fb21','colors',250
store_data,'fb22',time_double(fb22),replicate(1,2) & ylim,'fb22',0,2 & options,'fb22','psym',0 & options,'fb22','thick',4 & options,'fb22','colors',250
store_data,'fb23',time_double(fb23),replicate(1,2) & ylim,'fb23',0,2 & options,'fb23','psym',0 & options,'fb23','thick',4 & options,'fb23','colors',250
store_data,'fb24',time_double(fb24),replicate(1,2) & ylim,'fb24',0,2 & options,'fb24','psym',0 & options,'fb24','thick',4 & options,'fb24','colors',250
store_data,'schum',time_double(schum),replicate(0.8,2) & ylim,'schum',0,2 & options,'schum','psym',0 & options,'schum','thick',4 & options,'schum','colors',50

store_data,'campaigns',data=['b13','b14','b_kiruna1','b_kiruna2','b_kiruna3','fb1','fb2','fb3','fb4','fb5','fb6','fb7','fb8','fb9','fb10','fb11','fb12','fb13','fb14','fb15','fb16','fb17','fb18','fb19','fb20','fb21','fb22','fb23','fb24','schum'] & ylim,'campaigns',0,2
options,'campaigns','panel_size',0.8
options,'campaigns','ytitle','Burst collection!Ccampaigns!CBlack=BARREL!CRed=FIREBIRD!CBlue=Schumann'
ylim,'running_samples_monthly_rate_comb',0,600,0

tplot,['b1?comb','running_samples_monthly_rate_comb','running_samples_comb','campaigns']

;tplot,['running_samples_comb','campaigns']
;tplot,['running_samples_monthly_rate_comb','running_samples_comb']

stop

end
