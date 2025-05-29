;Determine if I need to flip sign in one of the Endurance probes 
;

.compile /Users/abrenema/Desktop/code/Aaron/github/spectral-coherence-phaselag-analysis/aaron_chaston_polarization.pro


rbsp_efw_init
timespan,'1970-01-01/00:00',20,/minutes

;------------------------
;Load skins
path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_skins/'
fn = '47001_TM1_LFDSP_S5Skins_V1SD2SD3SD4SD_cal.sav'
restore,path+fn

v1 = DV1S_VOLTS
v2 = DV2S_VOLTS
v3 = DV3S_VOLTS
v4 = DV4S_VOLTS

store_data,'v1',tv1s,v1
store_data,'v2',tv2s,v2
store_data,'v3',tv3s,v3
store_data,'v4',tv4s,v4

dif_data,'v1','v2'
dif_data,'v3','v4'

;V1 and V2 out of phase (true for large bumps and smallest wiggles)
tplot,['v1','v2']

;HOWEVER, V3 and V4 in phase (true for large bumps and smallest wiggles)
tplot,['v3','v4']




;------------------------
;Load differential DC

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_DC/'
fn = '47001_TM1_LFDSP_S5DCE_DCES5_calibrated.sav'
restore,path+fn

w12 = DV12_MVM
w34 = DV34_MVM
timesDC = times
store_data,'w12',data={x:time_double('1970-01-01/00:00:00') + timesDC,y:w12}
store_data,'w34',data={x:time_double('1970-01-01/00:00:00') + timesDC,y:w34}
;rbsp_detrend,'w12',0.001
;rbsp_detrend,'w34',0.001
;rbsp_detrend,'w12_smoothed',3.
;rbsp_detrend,'w34_smoothed',3.


tplot,['v1','v2','v3','v4','v1-v2','w12','v3-v4','w34']


tplot,['w12','w34']

