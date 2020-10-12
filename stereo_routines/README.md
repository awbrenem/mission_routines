# stereo_wave_id
Routines for grabbing and identifying STEREO burst captures


;  EXAMPLES:    This program is meant to be used in conjunction with
;				get_burst_groups.pro and get_group_elements.pro
;
;				sc = 'STA'
;				t0 = '2007-01-01/00:00:00'
;				t1 = '2007-12-31/24:00:00'
;				file2 = 'STEREO_TDS_STA_q_waveclass_2007_complete.txt'
;				tds = identify_waves(t0,t1,file='~/Desktop/code/Aaron/datafiles/stereo/wave_quality/'+file2)
;				groups = get_burst_groups('STA',tds,channel='Ex')
;				elements = get_group_elements('STA',tds,groups,channel='Ex')
;
