;Load and plot FAST data 

;http://sprg.ssl.berkeley.edu/fast/scienceops/fast_fields_help.html




rbsp_efw_init

;Make tplot plots looks pretty
charsz_plot = 0.8               ;character size for plots
charsz_win = 1.2
!p.charsize = charsz_win
tplot_options,'xmargin',[20.,15.]
tplot_options,'ymargin',[3,6]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',2
tplot_options,'ythick',2
tplot_options,'labflag',-1



;--------------------------------------------------------------------
;Determine what orbit number you're interested in (accurate to one orbit)

;orbitNo = fa_orbit_time('1997-01-19/15:00:00')
orbitNo = fa_orbit_time('1997-02-19/15:00:00')
print, orbitNo
  ;1635.77




;--------------------------------------------------------------------
;load e4k data

;Here is a list of variables in the FAST L2 despun e4k files. These
;include the E_NEAR_B and E_ALONG_V files, B phase and Sun phase, the
;spin plane field components E5-8 and E14-58, which are used to
;calculate E_NEAR_B and E_ALONG_V. All of the voltage differences and
;voltages avaliable from SDT are included.

'fa_e_near_b_4k': 'FAST 4k Burst mode electric field, in the direction corresponding to zero degrees magnetic phase, See:http://sprg.ssl.berkeley.edu/fast/scienceops/fast_fields_help.html'
'fa_e_along_v_4k': 'FAST 4k Burst mode electric field, in the direction corresponding to ninety degrees magnetic phase'
'fa_e1458_4k': 'FAST E12, probe 14-58 field, calibrated, not despun'
'fa_e58_4k': 'FAST E58, probe 5-8 field, calibrated, not despun'
'fa_bphase_4k': 'FAST B field phase, used for despin'
'fa_sphase_4k': 'FAST Sun phase, the projection of the sun on the spin plane, SPHASE = 0 corresponds to the X-axis of DSC coordinates'
'fa_v1_v2_4k': 'FAST probe 1-2 Voltage difference, direct SDT output'
'fa_v1_v4_4k': 'FAST probe 1-4 Voltage difference, direct SDT output'
'fa_v2_v4_4k': 'FAST probe 2-4 Voltage difference, direct SDT output'
'fa_v5_v8_4k': 'FAST probe 5-8 Voltage difference, direct SDT output'
'fa_v5_v6_4k': 'FAST probe 5-6 Voltage difference, direct SDT output'
'fa_v5_v7_4k': 'FAST probe 5-7 Voltage difference, direct SDT output'
'fa_v6_v8_4k': 'FAST probe 6-8 Voltage difference, direct SDT output'
'fa_v7_v8_4k': 'FAST probe 7-8 Voltage difference, direct SDT output'
'fa_v2_4k': 'FAST probe 2 Voltage, direct SDT output'
'fa_v6_4k': 'FAST probe 6 Voltage, direct SDT output'
'fa_v7_4k': 'FAST probe 7 Voltage, direct SDT output'
'fa_v9_4k': 'FAST probe 9 Voltage, direct SDT output'
'fa_v10_4k': 'FAST probe 10 Voltage, direct SDT output'

;--------------------------------------------------------------------


fa_despun_e_load, orbit = floor(orbitNo), datatype = 'e4k'

tplot_names
  ;1 fa_e_near_b_4k
  ;2 fa_e_along_v_4k
  ;3 fa_e1458_4k
  ;4 fa_e58_4k
  ;5 fa_bphase_4k
  ;6 fa_sphase_4k
  ;7 fa_v1_v4_4k
  ;8 fa_v5_v8_4k
  ;9 fa_v10_4k

get_data,'fa_e_near_b_4k',data=d




;--------
;Check sample rate and its consistency
print, 1/(d.x[1] - d.x[0])
;Sample rate = 8192.0000
df = 1/(d.x - shift(d.x,1))
df[0] = 0
store_data,'sampleRateTest',data={x:d.x,y:df}
options,'sampleRateTest','psym',4
ylim,'sampleRateTest',8150,8200
tplot,'sampleRateTest'
;--------


tplot,['fa_e_near_b_4k']



;See how V4 performs relative to other booms for AC








;load 16k
fa_despun_e_load, orbit = 1635, datatype = 'e16k'
;DSP data
fa_dsp_load, orbit = 1635
;SFA data
fa_sfa_load, orbit = 1635
tplot, ['fa_e_along_v_4k','fa_e_near_b_4k', 'fa_sfaave_mag3ac', 'fa_dspadc_mag3ac']

You can also pass in a time range:
fa_despun_e_load, orbit = 1635, datatype= 'e4k', trange = ['1996-12-11/00:00', '1996-12-11/04:00']
tplot,  ['fa_e_near_b_4k', 'fa_e_along_v_4k']




