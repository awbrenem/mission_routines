;Routines for producing the SWEA PAD and MOMENTS from the
;full distributions.
;****OBSOLETE. BENOIT LAVRUD TOLD ME NOT TO USE THESE. INSTEAD, USE
;THE OFFICIAL PADs. SEE READ_SWEA_PAD_FILES.PRO
;*****


timespan,'2017-03-24',3,/days
probe='a'

rbsp_efw_init
tplot_options,'xmargin',[20.,16.] & tplot_options,'ymargin',[3,9]
tplot_options,'xticklen',0.08 & tplot_options,'yticklen',0.02
tplot_options,'xthick',2 & tplot_options,'ythick',2
tplot_options,'labflag',-1


; load position and mag data
st_position_load,probe=probe
st_mag_load,probe=probe,coords='SC'


;This line is critical. If we don't rename the Bo variable then ST_SWEA_LOAD.pro
;(which calls st_swea_mag_load.pro) can't find the Bo data and can't define
;the pitch angle directions.
copy_data,'st'+probe+'_B_SC','st'+probe+'_l1_mag_sc'


; load SWEA 3D distribution (full distributions)
st_swea_load,probes=probe

; get the PAD
st_part_moments,probe=probe,/get_pad;,/get_moments

; reduce_pads usage:
; reduce_pads,'name',x,y,z
;   x=1 --> PAD across energy range specified by y,z.  y=z gives single energy channel
;   x=2 --> Energy distibution across PA range specified by y,z. "" single PA.

; break out PADs for each of the 16 energy channels
reduce_pads,'st'+probe+'_SWEA_pad',1,0,0    ; 1720 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,1,1    ; 1060 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,2,2    ;  650 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,3,3    ;  400 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,4,4    ;  250 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,5,5    ;  152 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,6,6    ;   93 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,7,7    ;   58 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,8,8    ;   35 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,9,9    ;   22 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,10,10  ;   13.4 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,11,11  ;    8.3 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,12,12  ;    5.1 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,13,13  ;    3.1 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,14,14  ;    1.93 eV
reduce_pads,'st'+probe+'_SWEA_pad',1,15,15  ;    1.19 eV


; break out energy spectrum for each of the 8 PA bins
reduce_pads,'st'+probe+'_SWEA_pad',2,0,0    ;   15.6 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,1,1    ;   35 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,2,2    ;   57 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,3,3    ;   79 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,4,4    ;  101 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,5,5    ;  123 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,6,6    ;  145 deg
reduce_pads,'st'+probe+'_SWEA_pad',2,7,7    ;  164 deg


; set up colors
device,decomposed=0
loadct,39

options,'*pad-2*','spec',1  ;  1 for spectral plots, 0 for line plots
ylim,'*pad-2*',1,2000,1  ; set up yrange

options,'*pad-*','no_interp',1  ; turn off interpolation

; plot PAD at for each energy bin
tplot,'*pad-1*'

; plot energy specs for each PA bin
;Ignore SWEA energy bins below ~50 eV as they can't be trusted
ylim,'*pad-2*',50,2000,1
zlim,'*pad-2*',1d2,1d9,1
tplot,'*pad-2*'
