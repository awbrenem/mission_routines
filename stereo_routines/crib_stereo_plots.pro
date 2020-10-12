;Plots for Cindy paper in 2019

;STEREO A, March 24th @12:00 through Mar-26, 2017
;-magnetic field
;-density
;-velocity
;-specific entropy
;-STE (downstream, ~1-1000 eV). Few different PA. Lance knows these, apparently (come as bins 0-79)
;-spectral waveform w/ zlim modified to see upstream waves clearly
;-Insert ~3 TDS captures/hodograms at bottom
;-TDSMax


;Get data from
;https://stereo-ssc.nascom.nasa.gov/data/ins_data/



timespan,'2017-03-24',3,/days
probe='a' ;

rbsp_efw_init
charsz_plot = 0.8             ;character size for plots
charsz_win = 1.2
!p.charsize = charsz_win
tplot_options,'xmargin',[20.,16.] & tplot_options,'ymargin',[3,9]
tplot_options,'xticklen',0.08 & tplot_options,'yticklen',0.02
tplot_options,'xthick',2 & tplot_options,'ythick',2
tplot_options,'labflag',-1

st_mag_load,probe=probe,coords='RTN'



;***Below routines basically don't work
;st_part_moments,probe='a'
;st_plastic_load,probe='a'
;st_swaves_load,probe='a'


;Load PLASTIC CDF files for moments
path = '/Users/aaronbreneman/Desktop/Research/OTHER/Stuff_for_other_people/Cattell_Cindy/2019_paper_plots/'
fn = 'STA_L2_PLA_1DMax_1min_20170324_V11.cdf' & cdf2tplot,path+fn
get_data,'proton_number_density',t1,d1
get_data,'proton_temperature',t1,temp1
get_data,'proton_Vr_RTN',t21,vr1 & get_data,'proton_Vt_RTN',t21,vt1 & get_data,'proton_Vn_RTN',t21,vn1
fn = 'STA_L2_PLA_1DMax_1min_20170325_V11.cdf' & cdf2tplot,path+fn
get_data,'proton_number_density',t2,d2
get_data,'proton_temperature',t2,temp2
get_data,'proton_Vr_RTN',t22,vr2 & get_data,'proton_Vt_RTN',t22,vt2 & get_data,'proton_Vn_RTN',t22,vn2
fn = 'STA_L2_PLA_1DMax_1min_20170326_V11.cdf' & cdf2tplot,path+fn
get_data,'proton_number_density',t3,d3
get_data,'proton_temperature',t3,temp3
get_data,'proton_Vr_RTN',t23,vr3 & get_data,'proton_Vt_RTN',t23,vt3 & get_data,'proton_Vn_RTN',t23,vn3

store_data,'temp',[t1,t2,t3],[temp1,temp2,temp3]
store_data,'density',[t1,t2,t3],[d1,d2,d3]
store_data,'Vrtn',[t21,t22,t23],[[vr1,vr2,vr3],[vt1,vt2,vt3],[vn1,vn2,vn3]]




;---------------------------------------------------
;Load SWEA and Mag data
path = '/Users/aaronbreneman/Desktop/Research/OTHER/Stuff_for_other_people/Cattell_Cindy/2019_paper_plots/'
tplot_restore,filename=path+'sta_20170324pads.tplot'
copy_data,'sta_pad_58eV','sta_pad_58eV_0324'
copy_data,'sta_pad_93eV','sta_pad_93eV_0324'
copy_data,'sta_pad_246eV','sta_pad_246eV_0324'
copy_data,'sta_pad_650eV','sta_pad_650eV_0324'
copy_data,'sta_pad_172.5deg','sta_pad_172.5deg_0324'
tplot_restore,filename=path+'sta_20170325pads.tplot'
copy_data,'sta_pad_58eV','sta_pad_58eV_0325'
copy_data,'sta_pad_93eV','sta_pad_93eV_0325'
copy_data,'sta_pad_246eV','sta_pad_246eV_0325'
copy_data,'sta_pad_650eV','sta_pad_650eV_0325'
copy_data,'sta_pad_172.5deg','sta_pad_172.5deg_0325'
tplot_restore,filename=path+'sta_20170326pads.tplot'
copy_data,'sta_pad_58eV','sta_pad_58eV_0326'
copy_data,'sta_pad_93eV','sta_pad_93eV_0326'
copy_data,'sta_pad_246eV','sta_pad_246eV_0326'
copy_data,'sta_pad_650eV','sta_pad_650eV_0326'
copy_data,'sta_pad_172.5deg','sta_pad_172.5deg_0326'


get_data,'sta_pad_58eV_0324',data=d1
get_data,'sta_pad_58eV_0325',data=d2
get_data,'sta_pad_58eV_0326',data=d3
store_data,'sta_pad_58eV',[d1.x,d2.x,d3.x],[d1.y,d2.y,d3.y],d1.v

get_data,'sta_pad_93eV_0324',data=d1
get_data,'sta_pad_93eV_0325',data=d2
get_data,'sta_pad_93eV_0326',data=d3
store_data,'sta_pad_93eV',[d1.x,d2.x,d3.x],[d1.y,d2.y,d3.y],d1.v

get_data,'sta_pad_246eV_0324',data=d1
get_data,'sta_pad_246eV_0325',data=d2
get_data,'sta_pad_246eV_0326',data=d3
store_data,'sta_pad_246eV',[d1.x,d2.x,d3.x],[d1.y,d2.y,d3.y],d1.v

get_data,'sta_pad_650eV_0324',data=d1
get_data,'sta_pad_650eV_0325',data=d2
get_data,'sta_pad_650eV_0326',data=d3
store_data,'sta_pad_650eV',[d1.x,d2.x,d3.x],[d1.y,d2.y,d3.y],d1.v

get_data,'sta_pad_172.5deg_0324',data=d1
get_data,'sta_pad_172.5deg_0325',data=d2
get_data,'sta_pad_172.5deg_0326',data=d3
store_data,'sta_pad_172.5deg',[d1.x,d2.x,d3.x],[d1.y,d2.y,d3.y],d1.v




;-------------------------------------------------
;Load STE CDF files
;STE data:    Array[8640, 8, 32]
;Eight channels are: STE-U 0    STE-U 1    STE-U 2    STE-U 3    STE-D 0    STE-D 1    STE-D 2    STE-D 3


st_ste_load,probe='a'
get_data,'sta_ste_D0',data=d
options,'sta_ste_D0','spec',1
;STE energy bins (keV)
; 1.94900      2.43921      2.95401      3.44438      3.95935      4.44990      4.96503      5.45574      5.97106
; 6.46193      6.97741      7.46846      8.05768      8.69606      9.43266      10.2675      11.2006      12.2811
; 13.5090      14.9334      16.6035      18.5683      20.9017      23.6773      27.0671      31.2432      36.4757
; 43.1088      51.7322      63.2059      79.0045      100.554

options,'sta_ste_D0','ytitle','STA STE flux D0!Cenergy [keV]'
tplot,'sta_ste_D0'


;Load CDF file directly:
;***NOTE: doesn't work properly. The structure isn't correct (missing "v" part)

fn = 'STA_L1_STE_20170324_V01.cdf'
pathtmp = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/stereo_routines/'
get_data,'STE_spectra',data=d1,dlim=dlim,lim=lim
;RBSP_EFW> help,d1,/st
;** Structure <168ba18>, 2 tags, length=8916480, data length=8916480, refs=1:
;   X               DOUBLE    Array[8640]
;   Y               FLOAT     Array[8640, 8, 32]

RBSP_EFW> help,dlim.cdf.vatt,/st
** Structure <168fdf8>, 14 tags, length=192, data length=188, refs=2:
   CATDESC         STRING    'STE Electron Spectra'
   DEPEND_0        STRING    'Epoch'
   DEPEND_1        STRING    'STE_spectra_LABL_1'
   DEPEND_2        STRING    'STE_detector'
   DISPLAY_TYPE    STRING    'time_series'
   FIELDNAM        STRING    'STE Electron Spectra Energy in units of cnts'
   FILLVAL         FLOAT      -1.00000e+31
   FORMAT          STRING    'E12.2'
   LABL_PTR_1      STRING    'STE_spectra_LABL_1'
   LABL_PTR_2      STRING    'STE_spectra_LABL_2'
   UNITS           STRING    'cnts'
   VALIDMIN        FLOAT           0.00000
   VALIDMAX        FLOAT       1.00000e+20
   VAR_TYPE        STRING    'data'




;fn = 'STA_L1_STE_20170324_V01.cdf' & cdf2tplot,path+fn
;get_data,'STE_spectra',data=d1
;fn = 'STA_L1_STE_20170325_V01.cdf' & cdf2tplot,path+fn
;get_data,'STE_spectra',data=d2
;fn = 'STA_L1_STE_20170326_V01.cdf' & cdf2tplot,path+fn
;get_data,'STE_spectra',data=d3

;store_data,'STE_spectra_D0',[d1.x,d2.x,d3.x],[reform(d1.y[*,4,*]),reform(d2.y[*,4,*]),reform(d3.y[*,4,*])]
;store_data,'STE_spectra_D1',[d1.x,d2.x,d3.x],[reform(d1.y[*,5,*]),reform(d2.y[*,5,*]),reform(d3.y[*,5,*])]
;store_data,'STE_spectra_D2',[d1.x,d2.x,d3.x],[reform(d1.y[*,6,*]),reform(d2.y[*,6,*]),reform(d3.y[*,6,*])]
;store_data,'STE_spectra_D3',[d1.x,d2.x,d3.x],[reform(d1.y[*,7,*]),reform(d2.y[*,7,*]),reform(d3.y[*,7,*])]
;options,'STE_spectra_D?','spec',1
;tplot,'STE_spectra_D?'














options,'sta_B_RTN','colors',[0,50,250]
options,'Vrtn','colors',[0,50,250]
zlim,'sta_ste_D?',1d-2,1d2,1



get_data,'sta_B_RTN',x2,brtn
bmag = sqrt(brtn[*,0]^2 + brtn[*,1]^2 + brtn[*,2]^2)
store_data,'bfield',x2,[[bmag],[reform(brtn[*,0])],[reform(brtn[*,1])],[reform(brtn[*,2])]]


tinterpol_mxn,'bfield','density',newname='bfield_interp'
get_data,'bfield_interp',data=dd
options,'bfield_interp','colors',[0,50,200,250]

store_data,'zeroline',dd.x,replicate(0.,n_elements(dd.x))
store_data,'bcomb',data=['bfield_interp','zeroline']


;-------------------------------
;THINGS TO ADD: note that the units should be the same as those in plastic files: n=cm-3, Tp=K, Bt=nT
;entropy=alog10((tp^1.5)/np) ) - 6*alog10(10);
;beta=np*10^(-5)*(469895.8 + 4.048*tp) / bt^2
;-------------------------------

get_data,'temp',x,tp
get_data,'bfield_interp',x2,b
get_data,'density',x3,np


bt = b[*,0]

v1 = (tp^1.5)/np
entropy=alog10(v1) - 6*alog10(10)
store_data,'entropy',x,entropy
options,'entropy','ytitle','entropy'

v2 = np*(10d^(-5.))
beta=v2*(469895.8 + 4.048*tp) / (bt^2)
store_data,'beta',x,beta
options,'beta','ytitle','beta'
ylim,'beta',0.1,100,1

get_data,'Vrtn',data=v
store_data,'Vrtn',v.x,[[v.y[*,0]],[v.y[*,1]],[v.y[*,2]]]

vmag_scaled = sqrt(v.y[*,0]^2 + v.y[*,1]^2 + v.y[*,2]^2)/10.
store_data,'Vmag_scaled',v.x,vmag_scaled
store_data,'Vmag_dens_comb',data=['Vmag_scaled','density']
options,'Vmag_dens_comb','colors',[0,50]


options,'bfield_interp','colors',[0,50,110,250]
options,'Vrtn','colors',[0,50,250]


options,'Vmag_dens_comb','ytitle','|V|/10 (km/s)!Cblack!CDensity (cm-3)!Cblue'
options,'bfield_interp','ytitle','Bfield [nT]!Cmag=black!CBr=darkblue!CBt=lightblue!CBn=red'
options,'bcomb','ytitle','Bfield [nT]!Cmag=black!CBr=darkblue!CBt=lightblue!CBn=red'
options,'Vrtn','ytitle','Velocity RTN!CVr=black!CVt=blue!CVn=red'
options,'temp','ytitle','Proton Temp!CK'
options,'beta','ytitle','Plasma beta'




options,'*','panel_size',1
options,'sta_pad_172.5deg','spec',1


ylim,'STE_spectra_D?',1,10,0
zlim,'STE_spectra_D?',1,100,1
ylim,'beta',0.1,100,1
ylim,'sta_pad_172.5deg',50,1000,1
ylim,'sta_ste_D0',2,6,0

zlim,'sta_pad_58eV',10,1000,1
zlim,'sta_pad_93eV',10,1000,1
zlim,'sta_pad_246eV',1,100,1
zlim,'sta_pad_650eV',0.01,1,1
zlim,'sta_pad_172.5deg',1d0,1d3,1

tplot,['bcomb',$
'Vmag_dens_comb',$
'entropy',$
'beta',$
'sta_pad_93eV','sta_pad_246eV','sta_pad_650eV',$
'sta_pad_172.5deg',$
'sta_ste_D0']


;'STE_spectra_D0']



stop
;get_data,'density[Plastic]',x3,np
;'SWEA_dist_0',$
;'SWEA_dist_79',$
;'sta_SWEA_pad-1-7:7',$
;'sta_SWEA_pad-1-4:4',$
;'sta_SWEA_pad-1-2:2',$
;'sta_SWEA_pad-2-7:7',$


;sta_B_RTN
;sta_ste_D?









;Load the TDSMax data
;lrshk2tplot,sc='a',date='2017-03-24'



end
