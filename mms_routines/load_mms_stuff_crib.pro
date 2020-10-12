;Plot interesting MMS data. 
;Aaron's custom cribsheet, based primarily on the MMS crib sheets. 


;Notes:
; srvy = survey rate data 
; fast = fast rate data
; brst = burst rate data
; "e" = electron, e.g. "des"
; "i" = ion, e.g. "dis"



; FPI = fast plasma investigation (10 eV to 30 keV e- and ions)
; see https://lasp.colorado.edu/mms/sdc/public/datasets/fpi/FPI_dpg.pdf
    ;des = dual electron spectrometers
    ;fast = 4.5 sec


;Make pretty tplot plots
  !p.charsize = 1.0
  tplot_options,'xmargin',[20.,16.]
  tplot_options,'ymargin',[3,9]
  tplot_options,'xticklen',0.08
  tplot_options,'yticklen',0.02
  tplot_options,'xthick',2
  tplot_options,'ythick',2
  tplot_options,'labflag',-1







probe = '2'
level = 'l2'

date = '2017-12-10'
;t0 = time_string(date + '/' + '02:48:13')
;t1 = time_string(date + '/' + '02:58')
t0 = time_string(date + '/' + '02:30')
t1 = time_string(date + '/' + '02:58')
trange = [t0,t1]


t0d = time_double(t0)
t1d = time_double(t1)
timespan,t0d,t1d-t0d,/seconds



;---------------------------------------------------------
;Load DC magnetic field data
;---------------------------------------------------------

mms_load_fgm,probe=probe

tplot,['mms2_fgm_b_gse_srvy_l2_bvec',$
        'mms2_fgm_b_gse_srvy_l2_btot']

mms_load_fgm,probe=probe,data_rate='brst'

tplot,['mms2_fgm_b_gse_brst_l2_bvec',$
        'mms2_fgm_b_gse_brst_l2_btot']





;---------------------------------------------------------
;Load fast particle data
;   From mms_load_fpi_crib.pro
;---------------------------------------------------------


; DES/DIS moments file (contains moments, as well as spectra and pitch angle distributions)


mms_load_fpi,probes=probe,datatype=['des-moms', 'dis-moms'],$
    level='l2',data_rate='fast',min_version='2.2.0'
;mms_load_fpi,probes=probe,datatype='dis-moms',$
;    level='l2',data_rate='fast',min_version='2.2.0'


prefix = 'mms'+strcompress(string(probe), /rem)

tplot,prefix + '_' + ['des_pitchangdist_avg',$
'dis_energyspectr_omni_fast',$
'dis_numberdensity_fast',$
'des_numberdensity_fast']
tlimit,trange


name =  'mms'+probe+'_des_dist_fast'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_des_bulkv_gse_fast'     ;name of bulk velocity vector
dist_fast_e = mms_get_dist(name, trange=trange,data_rate='fast')

name =  'mms'+probe+'_dis_dist_fast'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_dis_bulkv_gse_fast'     ;name of bulk velocity vector
dist_fast_i = mms_get_dist(name, trange=trange,data_rate='fast')



;--------------------------------
;Determine energy ranges of plots
;--------------------------------

get_data,'mms2_des_pitchangdist_lowen_fast',dlim=dlim
print,dlim.cdf.vatt.var_notes
;low energy bin: 0 eV - 200 eV. pitch-angle bin size: 6 deg. 
get_data,'mms2_des_pitchangdist_miden_fast',dlim=dlim
print,dlim.cdf.vatt.var_notes
;mid energy bin: 200 eV - 2 keV.  pitch-angle bin size: 6 deg. 
get_data,'mms2_des_pitchangdist_highen_fast',dlim=dlim
print,dlim.cdf.vatt.var_notes
;high energy bin: 2 keV - 30 keV.  pitch-angle bin size: 6 deg. 

options,'mms2_des_pitchangdist_lowen_fast','ytitle','mm2 des dist!Cfast [0-200 eV]!Cs^3/cm^6'
options,'mms2_des_pitchangdist_miden_fast','ytitle','mm2 des dist!Cfast [200-2000 eV]!Cs^3/cm^6'
options,'mms2_des_pitchangdist_highen_fast','ytitle','mm2 des dist!Cfast [2000-30000 eV]!Cs^3/cm^6'

;---------------------------------------------------------
;Load burst particle data
;---------------------------------------------------------

level='l2'
data_rate='brst'

;electron burst
mms_load_fpi, data_rate=data_rate, level=level, datatype='des-dist', $
              probe=probe, trange=trange, min_version='2.2.0'

name =  'mms'+probe+'_des_dist_'+data_rate
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_des_bulkv_gse_'+data_rate     ;name of bulk velocity vector
dist_burst_e = mms_get_dist(name, trange=trange,data_rate='brst')


;ion burst
name =  'mms'+probe+'_dis_dist_'+data_rate
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_dis_bulkv_gse_'+data_rate     ;name of bulk velocity vector
mms_load_fpi, data_rate=data_rate, level=level, datatype='dis-dist', $
              probe=probe, trange=trange, min_version='2.2.0'

dist_burst_i = mms_get_dist(name, trange=trange,data_rate='brst')



;load burst velocity moment
mms_load_fpi, data_rate=data_rate, level=level, datatype=['des-moms','dis-moms'], $
              probe=probe, trange=trange, min_version='2.2.0'
;---------------------------------------------




options,'mms2_des_pitchangdist_lowen_brst','ytitle','mm2 des dist!Cbrst [0-200 eV]!Cs^3/cm^6'
options,'mms2_des_pitchangdist_miden_brst','ytitle','mm2 des dist!Cbrst [200-2000 eV]!Cs^3/cm^6'
options,'mms2_des_pitchangdist_highen_brst','ytitle','mm2 des dist!Cbrst [2000-30000 eV]!Cs^3/cm^6'


;--------------------------------------------
;Load EDP (electric field) burst data 
;--------------------------------------------
;mms_load_edp

mms_load_edp, data_rate='brst', probes=probe, level='l2'

; Display colors for parallel E (black) and error (pink)
; Large error bars signifies possible presence of cold plasma
; or spacecraft charging, which can make axial electric field
; measurements difficult. Please always use error bars on e-parallel!!
options, 'mms?_edp_dce_par_epar_fast_l2', colors = [1, 0]
options, 'mms?_edp_dce_par_epar_fast_l2', labels = ['Error', 'E!D||!N']

; Since the electric field is often close to zero in multiple components, label spacing tends to get bunched
; together
options, '*', 'labflag', -1
split_vec,'mms2_edp_dce_gse_brst_l2'

;tplot, ['mms2_edp_dce_gse_brst_l2','mms2_edp_dce_dsl_brst_l2']

;dif_data,'mms2_edp_dce_gse_brst_l2','mms2_edp_dce_dsl_brst_l2',newname='ediff'
;tplot,'ediff'

;--------------------------------------------
;Load SCM burst data 
;--------------------------------------------

;; Select data rate ('srvy' or 'burst')
scm_data_rate = 'brst';'brst';'srvy'

;; Select mode ('scsrvy' for survey data rate (both slow and fast have 32 S/s), 
;                'scb' (8192 S/s) or 'schb' (16384 S/s) for burst data rate)
scm_datatype = 'scb';'scb';'scsrvy'

scm_name = 'mms'+probe+'_scm_acb_gse_'+scm_datatype+'_'+scm_data_rate+'_l2'

mms_load_scm, trange=trange, probes=probe, level='l2', data_rate=scm_data_rate, datatype=scm_datatype, tplotnames=tplotnames



options, scm_name, colors=[2, 4, 6]
options, scm_name, labels=['X', 'Y', 'Z']
options, scm_name, labflag=-1


window, 0, ysize=650
tplot_options, 'xmargin', [15, 15]
tplot_options,title= 'MMS'+probe+' '+ scm_data_rate+' period, '+scm_datatype +' SCM data in GSE frame'

; plot the SCM data
;tplot, scm_name
;tlimit,trange


; calculate the dynamic power spectra without overlapping nshiftpoints=nboxpoints
if scm_datatype eq 'scb' then nboxpoints_input = 8192 else nboxpoints_input = 512

tdpwrspc, scm_name, nboxpoints=nboxpoints_input,nshiftpoints=nboxpoints_input,bin=1

if scm_datatype eq 'scsrvy' then Fmin = 0.5
if scm_datatype eq 'scsrvy' then Fmax = 16.
if scm_datatype eq 'scb'    then Fmin = 1.
if scm_datatype eq 'scb'    then Fmax = 4096.
if scm_datatype eq 'schb'   then Fmin = 32.
if scm_datatype eq 'schb'   then Fmax = 8192.

options, scm_name+'_?_dpwrspc', 'ytitle', 'MMS'+probe+' '+scm_datatype
options, scm_name+'_x_dpwrspc', 'ysubtitle', 'dynamic power!CX!C[Hz]'
options, scm_name+'_y_dpwrspc', 'ysubtitle', 'dynamic power!CY!C[Hz]'
options, scm_name+'_z_dpwrspc', 'ysubtitle', 'dynamic power!CZ!C[Hz]'
options, scm_name+'_?_dpwrspc', 'ztitle', '[nT!U2!N/Hz]

ylim, scm_name+'_?_dpwrspc',Fmin,Fmax,1
 
tplot, [scm_name, scm_name+'_?_dpwrspc']




;---------------------------------------------------------
;Load fast wave data from EDI (electron drift instrument)
;---------------------------------------------------------


;; load the E-field data
;mms_load_edi, probes=probe, data_rate='brst', datatype='efield', level='l2'
;mms_load_edi, probes=probe, data_rate='fast', datatype='efield', level='l2'

;; set the colors
;options, 'mms'+probe+'_edi_*_gsm_srvy_l2', colors=[2, 4, 6]
;
;; plot the data
;tplot, 'mms'+probe+['_edi_e_gsm_srvy_l2',$
;                    '_edi_vdrift_gsm_srvy_l2'] ; ExB drift velocity


;----------------------------------------
;Load fast/burst particle data 
;----------------------------------------
;From mms_slice2d_fpi_crib.pro
;Field-aligned slice part

;         datatype:     valid datatypes are:
;                         Quicklook: ['des', 'dis'] 
;                         SITL: '' (none; loads both electron and ion data from single CDF)
;                         L1b/L2: ['des-dist', 'dis-dist', 'dis-moms', 'des-moms']
;         data_rate:    instrument data rates for MMS FPI include 'fast', 'brst'. 




;--------------------------------------------
;Load heavier ion data
;--------------------------------------------

stop

mms_load_hpca, probe=probe, data_rate='srvy', level='l2', datatype='ion', versions=hpca_versions
mms_load_hpca, probe=probe, data_rate='srvy', level='l2', datatype='moments', versions=hpca_versions
;mms_load_hpca, probe=probe, data_rate='brst', level='l2', datatype='moments', versions=hpca_versions

; sum the HPCA spectra over the full field of view
mms_hpca_calc_anodes, fov=[0, 360], probe=probe


ylim,['mms2_hpca_hplus_flux','mms2_hpca_heplus_flux','mms2_hpca_heplusplus_flux','mms2_hpca_oplus_flux'],1,10,1
tplot,['mms2_hpca_hplus_flux','mms2_hpca_heplus_flux','mms2_hpca_heplusplus_flux','mms2_hpca_oplus_flux']



;-------------------
;Heavy ion PADs

level = 'l2'
species = 'oplus'
data_rate = 'srvy'

name = 'mms'+probe+'_hpca_'+species+'_phase_space_density'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec'             ;name of bfield vector
vname = 'mms'+probe+'_hpca_'+species+'_ion_bulk_velocity' ;name of bulk velocity vector

;timespan, '2015-10-16/13:06:00', 1, /min  ;time range to load
;trange = timerange()
;time = trange[0]  ;slice time 

mms_load_hpca, probes=probe, trange=trange[0], data_rate=data_rate, level=level, datatype='ion'

dist = mms_get_dist(name)

;load B field data
mms_load_fgm, probe=probe, trange=timerange(), level='l2'

;load velocity moment
mms_load_hpca, probes=probe, trange=timerange(), data_rate=data_rate, level=level, $
               datatype='moments', varformat='*_'+species+'_ion_bulk_velocity'

slice = spd_slice2d(dist, time=trange[0], window=window, $
                    rotation='bv', mag_data=bname, vel_data=vname)

;plot
spd_slice2d_plot, slice

;Plot slice
; SLICE_NORM='tvar1'

;-----------------------------------------------
;Run cribsheet crib_master_v5.pro 
;-----------------------------------------------

;mms1_fgm_b_gsm_srvy_l2_btot
;mms1_fgm_b_gsm_srvy_clipped
;mms1_dis_energyspectr_omni_fast
;mms1_des_energyspectr_omni_fast
;mms1_dis_numberdensity_fast
;mms1_dis_bulkv_gse_fast
;mms1_hpca_hplus_flux_elev_0-360
;mms1_hpca_heplusplus_flux_elev_0-360
;mms1_hpca_oplus_flux_elev_0-360
;mms1_hpca_hplus_number_density
;mms1_hpca_hplus_ion_bulk_velocity_GSM
;mms1_hpca_oplus_ion_bulk_velocity_GSM
;mms1_hpca_hplus_scalar_temperature
;mms1_edp_dce_gse_fast_l2
;mms1_dsp_epsd_omni
;mms1_dsp_bpsd_omni_fast_l2    ;'Omni-directional magnetic power spectral density: square root of the sum o'...
;mms1_epd_feeps_srvy_l2_electron_intensity_omni   ;'Unidirectional differential flux per spin sector Bottom Sens12'





;--------------------------------------------
;Various plots
;--------------------------------------------

store_data,'mms2_des_temps_fast_comb',data=['mms2_des_temppara_fast','mms2_des_tempperp_fast']


;Survey plot electrons
tplot,['mms2_edp_dce_gse_brst_l2_y',$
        'mms2_scm_acb_gse_scb_brst_l2_x',$
        'mms2_fgm_b_gse_srvy_l2_bvec',$
        'mms2_fgm_b_gse_srvy_l2_btot',$
        'mms2_des_numberdensity_fast',$
        'mms2_des_pitchangdist_lowen_fast',$
        'mms2_des_pitchangdist_miden_fast',$
        'mms2_des_pitchangdist_highen_fast',$
        'mms2_des_energyspectr_omni_fast']

;Survey plot ions
tplot,['mms2_edp_dce_gse_brst_l2_y',$
'mms2_scm_acb_gse_scb_brst_l2_x',$
'mms2_fgm_b_gse_srvy_l2_bvec',$
        'mms2_fgm_b_gse_srvy_l2_btot',$
        'mms2_dis_energyspectr_px_fast',$
        'mms2_dis_energyspectr_mx_fast',$
        'mms2_dis_energyspectr_py_fast',$
        'mms2_dis_energyspectr_my_fast',$
        'mms2_dis_energyspectr_pz_fast',$
        'mms2_dis_energyspectr_mz_fast',$
        'mms2_dis_energyspectr_omni_fast',$
'mms2_hpca_hplus_flux_elev_0-360',$
'mms2_hpca_heplus_flux_elev_0-360',$
'mms2_hpca_heplusplus_flux_elev_0-360',$
'mms2_hpca_oplus_flux_elev_0-360']

;        'mms2_hpca_oplus_tperp','mms2_hpca_oplos_tparallel']

        



ylim,'mms2_des_energyspectr_omni_brst',6,1000,1


;Burst plot 
tplot,['mms2_edp_dce_gse_brst_l2_y',$
'mms2_scm_acb_gse_scb_brst_l2_x',$
scm_name+'_x_dpwrspc',$
'mms2_des_dist_brst',$
'mms2_fgm_b_gse_brst_l2_bvec',$
'mms2_fgm_b_gse_brst_l2_btot',$
'mms2_des_numberdensity_brst',$
'mms2_des_pitchangdist_lowen_brst',$
'mms2_des_pitchangdist_miden_brst',$
'mms2_des_pitchangdist_highen_brst',$
'mms2_des_energyspectr_omni_brst']


tplot,['mms2_edp_dce_gse_brst_l2_y',$
'mms2_scm_acb_gse_scb_brst_l2_x',$
'mms2_des_pitchangdist_lowen_brst',$
'mms2_des_pitchangdist_miden_brst',$
'mms2_des_pitchangdist_highen_brst',$
'mms2_des_energyspectr_omni_brst']


;burst plot ions
tplot,['mms2_edp_dce_gse_brst_l2_y',$
        'mms2_scm_acb_gse_scb_brst_l2_x',$
        'mms2_fgm_b_gse_srvy_l2_bvec',$
        'mms2_fgm_b_gse_srvy_l2_btot',$
        'mms2_dis_energyspectr_px_brst',$
        'mms2_dis_energyspectr_mx_brst',$
        'mms2_dis_energyspectr_py_brst',$
        'mms2_dis_energyspectr_my_brst',$
        'mms2_dis_energyspectr_pz_brst',$
        'mms2_dis_energyspectr_mz_brst']



;--------------------------------------------
;Plot energy or velocity contours for interesting times
;--------------------------------------------

;Ion fast distributions
time = '2017-12-10/02:51:41'  ;ES wave and ion beam
time = '2017-12-10/02:52:12'  ;

name =  'mms'+probe+'_dis_dist_fast'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_dis_bulkv_gse_fast'     ;name of bulk velocity vector


window = 60.
erange = [1,30000.]  ;eV
slice = spd_slice2d(dist_fast_i, time=time, window=window, $
                    rotation='perp', mag_data=bname, vel_data=vname,$
                    /energy,erange=erange,/subtract_bulk)

spd_slice2d_plot, slice,/sundir,/plotbfield


;Create a series of png plots for trange

window = 5.
nelem = (time_double(t1) - time_double(t0))/window
times = window*dindgen(nelem) + time_double(t0)
;print,time_string(times)



for i=0, n_elements(times)-1 do begin ;$
  slice = spd_slice2d(dist_fast_i, time=times[i], window=window, rotation='bv', mag_data=bname, vel_data=vname,/subtract_bulk,/energy,erange=erange); & $
  filename = '~/Desktop/mms'+probe+'_i_'+time_string(times[i],format=2); & $ 
  spd_slice2d_plot, slice, export=filename ;,/eps
endfor

;velocity
for i=0, n_elements(times)-1 do begin ;$
  slice = spd_slice2d(dist_fast_i,time=times[i],window=window,rotation='bv',mag_data=bname,vel_data=vname,/subtract_bulk) ;& $
  filename = '~/Desktop/mms'+probe+'_i_'+time_string(times[i],format=2) ;& $ 
  spd_slice2d_plot, slice, export=filename ;,/eps
endfor




;Electron fast distributions
time = '2017-12-10/02:51:41'  ;ES wave and ion beam
time = '2017-12-10/02:52:12'  ;

name =  'mms'+probe+'_des_dist_fast'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_des_bulkv_gse_fast'     ;name of bulk velocity vector


window = 60.
erange = [1,30000.]  ;eV
slice = spd_slice2d(dist_fast_e, time=time, window=window, $
                    rotation='bv', mag_data=bname, vel_data=vname,$
                    /energy,erange=erange,/subtract_bulk)

spd_slice2d_plot, slice



;Create a series of png plots for trange
window = 5.
nelem = (time_double(t1) - time_double(t0))/window
times = window*dindgen(nelem) + time_double(t0)
;print,time_string(times)



for i=0, n_elements(times)-1 do begin ;$
  slice = spd_slice2d(dist_fast_e, time=times[i], window=window, rotation='bv', mag_data=bname, vel_data=vname,/subtract_bulk,/energy,erange=erange) ;& $
  filename = '~/Desktop/mms'+probe+'_e_'+time_string(times[i],format=2) ;& $ 
  spd_slice2d_plot, slice, export=filename ;,/eps
endfor 

;velocity
for i=0, n_elements(times)-1 do begin ;$
  slice = spd_slice2d(dist_fast_e,time=times[i],window=window,rotation='bv',mag_data=bname,vel_data=vname,/subtract_bulk) ;& $
  filename = '~/Desktop/mms'+probe+'_e_'+time_string(times[i],format=2) ;& $ 
  spd_slice2d_plot, slice, export=filename,xrange=[-2d4,2d4],yrange=[-2d4,2d4] ;,/eps
endfor



















;Ion burst distributions

;time = '2017-12-10/02:49:15'
time = '2017-12-10/02:52:12'  ;

name =  'mms'+probe+'_dis_dist_brst'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_dis_bulkv_gse_brst'     ;name of bulk velocity vector

window = 0.5
erange = [1,30000.]  ;eV
slice = spd_slice2d(dist_burst_i, time=time, window=window, $
                    rotation='perp', mag_data=bname, vel_data=vname,$
                    /subtract_bulk,/energy,erange=erange)
;slice = spd_slice2d(dist_burst, time=time, window=window, $
;                    rotation='bv', mag_data=bname, vel_data=vname,$
;                    /subtract_bulk)

spd_slice2d_plot, slice,/plotbulk


;t0z = '2017-12-10/02:51:38'  ;
;t1z = '2017-12-10/02:51:50'  ;
t0z = '2017-12-10/02:51:00'  ;
t1z = '2017-12-10/02:52:30'  ;


window = 0.5
nelem = (time_double(t1z) - time_double(t0z))/window
times = window*dindgen(nelem) + time_double(t0z)

for i=0, n_elements(times)-1 do begin $
  slice = spd_slice2d(dist_burst_i, time=times[i], window=window,rotation='bv',mag_data=bname,vel_data=vname,/subtract_bulk,/energy,erange=erange) & $
  filename = '~/Desktop/mms'+probe+'_i_'+time_string(times[i],format=2) & $
  spd_slice2d_plot, slice, export=filename,/plotbulk


;  endfor





;Electron burst distributions

;time = '2017-12-10/02:49:15'
time = '2017-12-10/02:54:10'  ;

name =  'mms'+probe+'_des_dist_brst'
bname = 'mms'+probe+'_fgm_b_gse_srvy_l2_bvec' ;name of bfield vector
vname = 'mms'+probe+'_des_bulkv_gse_brst'     ;name of bulk velocity vector

window = 0.5
erange = [1,3000.]  ;eV
slice = spd_slice2d(dist_burst_e, time=time, window=window, $
                    rotation='perp', mag_data=bname, vel_data=vname,$
                    /subtract_bulk,/energy,erange=erange)
;slice = spd_slice2d(dist_burst_e, time=time, window=window, $
;                    rotation='bv', mag_data=bname, vel_data=vname,$
;                    /subtract_bulk)

spd_slice2d_plot, slice,xrange=[-2d4,2d4],yrange=[-2d4,2d4],/sundir,/plotbfield


;t0z = '2017-12-10/02:51:38'  ;
;t1z = '2017-12-10/02:51:50'  ;
t0z = '2017-12-10/02:54:00'  ;
t1z = '2017-12-10/02:54:30'  ;

;t0z = '2017-12-10/02:49:20'
;t1z = '2017-12-10/02:49:30'


window = 1.
nelem = (time_double(t1z) - time_double(t0z))/window
times = window*dindgen(nelem) + time_double(t0z)
erange = [1,3000.]  ;eV

for i=0, n_elements(times)-1 do begin ;$
  slice = spd_slice2d(dist_burst_e, time=times[i], window=window,rotation='bv',mag_data=bname,vel_data=vname,/subtract_bulk,/energy,erange=erange) ;& $
  filename = '~/Desktop/mms'+probe+'_e_'+time_string(times[i],format=2) ;& $
  spd_slice2d_plot, slice, export=filename ;,/eps
endfor


;--------------------------------------------
;Find angle b/t Bo and Vsw 
;--------------------------------------------

tplot,['mms2_des_bulkv_gse_fast','mms2_fgm_b_gse_srvy_l2_bvec']

tinterpol_mxn,'mms2_fgm_b_gse_srvy_l2_bvec','mms2_des_bulkv_gse_fast',newname='mms2_fgm_b_gse_srvy_l2_bvec_interp'

get_data,'mms2_des_bulkv_gse_fast',data=vsw
get_data,'mms2_fgm_b_gse_srvy_l2_bvec_interp',data=bo

vmag = sqrt(vsw.y[*,0]^2 + vsw.y[*,1]^2 + vsw.y[*,2]^2)
bmag = sqrt(bo.y[*,0]^2 + bo.y[*,1]^2 + bo.y[*,2]^2)

angle = fltarr(n_elements(vmag))

for i=0,n_elements(angle)-1 do angle[i] = acos(total(vsw.y[i,*]*bo.y[i,*])/vmag[i]/bmag[i])/!dtor

store_data,'angle_vsw_bo',vsw.x,angle
tplot,['mms2_des_bulkv_gse_fast','mms2_fgm_b_gse_srvy_l2_bvec_interp','angle_vsw_bo']


;--------------------------------------
;Different kind of slices to plot
;--------------------------------------


;get single distribution
;  -3d/2d interpolation show smooth contours
;  -3d interpolates entire volume
;  -2d interpolates projection of a subset of data near the slice plane 
;  -geometric interpolation is slow but shows bin boundaries
;---------------------------------------------
;slice = spd_slice2d(dist, time=time) ;3D interpolation
;slice = spd_slice2d(dist, time=time, /two) ;2D interpolation
;slice = spd_slice2d(dist, time=time, /geo) ;geometric interpolation

;average all data in specified time window
;slice = spd_slice2d(dist, time=time, /geo, window=0.5)  ; window (sec) starts at TIME  
;slice = spd_slice2d(dist, time=time, /geo, window=0.5, /center_time)  ; window centered on TIME
 
;average specific number of distributions (uses N closest to specified time)
;slice = spd_slice2d(dist, time=time, /geo, samples=3)

stop


;======================================================================
; Export time series
;======================================================================

name =  'mms'+probe+'_d'+species+'s_dist_'+data_rate
trange = [t0,t1] ;trange=['2015-10-16/13:06', '2015-10-16/13:07']


;produce a plot of 0.5 seconds of data every 10 seconds for 1 minute
times = time_double(trange[0]) + 10 * findgen(7)
window = 60.

for i=0, n_elements(times)-1 do begin 
  slice = spd_slice2d(dist, time=times[i], window=window) 
  filename = 'mms'+probe+'_'+species+'_'+time_string(times[i],format=2)
  spd_slice2d_plot, slice, export=filename ;,/eps
endfor

end



;--------------------------------------------
;Various things that can be plotted
;--------------------------------------------

;;Determine energies for particle plots 
;ylim,['mms2_des_energy_fast','mms2_des_energy_brst'],5,3d4,1
;tplot,['mms2_des_energy_fast','mms2_des_energy_brst']
;tplot,['mms2_des_dist_fast','mms2_des_energy_fast']
;;4.5 sec electron spectrometers low, mid, and high energy PA spectra
;tplot,['mms2_des_pitchangdist_lowen_fast',$
;'mms2_des_pitchangdist_miden_fast',$
;'mms2_des_pitchangdist_highen_fast',$
;'mms2_des_pitchangdist_avg']
;tplot,['mms2_des_energyspectr_par_fast',$
;'mms2_des_energyspectr_anti_fast',$
;'mms2_des_energyspectr_perp_fast',$
;'mms2_des_energyspectr_omni_fast']
;tplot,['mms2_des_numberdensity_fast',$
;'mms2_des_densityextrapolation_low_fast',$
;'mms2_des_densityextrapolation_high_fast']
;tplot,['mms2_des_bulkv_gse_fast','mms2_fgm_b_gse_srvy_l2_bvec']
;tplot,['mms2_des_heatq_gse_fast',$
;'mms2_des_temppara_fast',$
;'mms2_des_tempperp_fast']
;tplot,['mms2_des_phi_brst',$
;'mms2_des_dist_brst',$
;'mms2_des_energy_brst',$
;'mms2_des_energy_delta_brst']
;tplot,['mms2_des_pitchangdist_lowen_brst',$
;'mms2_des_pitchangdist_miden_brst',$
;'mms2_des_pitchangdist_highen_brst']
;tplot,['mms2_des_energyspectr_par_brst',$
;'mms2_des_energyspectr_anti_brst',$
;'mms2_des_energyspectr_perp_brst',$
;'mms2_des_energyspectr_omni_brst']
;tplot,['mms2_des_numberdensity_brst',$
;'mms2_des_numberdensity_err_brst',$
;'mms2_des_densityextrapolation_low_brst',$
;'mms2_des_densityextrapolation_high_brst']
;tplot,['mms2_des_prestensor_gse_brst',$
;'mms2_des_temptensor_gse_brst',$
;'mms2_des_heatq_gse_brst']


