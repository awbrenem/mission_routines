;------------------------------
;Rotate the dbm_scm UVW data into SC coordinates. 
;Returns tplot variable with '_SC' appended to name

;Rotation array from Bowen20, eqn 8, which rotates UVW into MAG coordinates. 
;They note that this is equivalent to SC coordinates. 

;This has been tested by direct comparison to the merged product
;psp_fld_l3_merged_scam_wf_2018110300_v01.cdf
;which comes in SC coord. 
;After proper filtering, the results are very close 
;---see below for the test. 


pro psp_rotate_scm_uvw2sc,var


    rot_mat_scm = [[0.81654,-0.40827,-0.40827],$
                    [0.,-0.70715,0.70715],$
                    [-0.57729,-0.57729,-0.57729]]



    get_data,var,data=dd
    Bsc = reform(rot_mat_scm ## dd.y)
    store_data,var + '_SC',dd.x,Bsc




  ;*********************************************************************************
  ;*********************************************************************************
  ;*********************************************************************************
  ;Test SCM data in SC coord against merged data product that comes in SC coord
  if keyword_set(test_compare) then begin 

  ;Load the 146 S/sec merged data
  cdf2tplot,'/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/mission_routines/solar_probe_routines/psp_fld_l3_merged_scam_wf_2018110300_v01.cdf'
  get_data,'psp_fld_l3_merged_scam_wf_SC',data=bw


  ;Let's bandpass the burst data to put it in a similar range to the lower cadence merged data
  lf = 15.  ;Hz
  hf = 80.  ;Hz

  ;First bandpass the merged data
  srt = 1/(bw.x[1]-bw.x[0])  ;sample rate
  x = rbsp_vector_bandpass(bw.y,srt,lf,hf)  ;automatically applies Hanning window
  store_data,'psp_fld_l3_merged_scam_wf_SC_bp',bw.x,x

  ;Now bandpass the burst data
;  get_data,'burst_psp_fld_l2_dfb_dbm_scm_SC_018-11-03/00:12:01.392820',data=bwb
  get_data,'burst_psp_fld_l2_dfb_dbm_scm_2018-11-03/00:12:01.392820_SC',data=bwb
  srt = 1/(bwb.x[1]-bwb.x[0])  ;sample rate
  xb = rbsp_vector_bandpass(bwb.y,srt,lf,hf)  ;automatically applies Hanning window
  store_data,'burst_psp_fld_l2_dfb_dbm_scm_018-11-03/00:12:01.392820_SC_bp',bwb.x,xb

  ;The two compare very closely 
  tplot,['burst_psp_fld_l2_dfb_dbm_scm_018-11-03/00:12:01.392820_SC_bp',$
        'psp_fld_l3_merged_scam_wf_SC_bp']
  ;**********************************************************************************
  endif






end




