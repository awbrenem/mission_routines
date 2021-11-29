;for testing rbsp_efw_create_esvy_uvw_from_vsvy.pro

;testing = 0.
;date='2015-03-17'  ;test event in monthly reports (won't work for method 2 b/c only V2 and V4 are not misbehaving)
date='2015-11-05'  ;Jinxing and Bortnik's event  5:25-8:40 UT
;date='2015-11-13'  ;Xiaojia's event 13:01-13:51 UT.
;date='2015-10-12'
;date='2013-12-30' ;little-no wave activity
;date='2013-12-18'  ;some wave activity
;date='2013-12-08'
;date='2015-12-31' ;little-no wave activity
;date='2015-12-21' ;Big 20 mV/m spike



timespan,date
probe = 'a'
bad_probe = 1


tplot_options,'xmargin',[20.,16.]
tplot_options,'ymargin',[3,9]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',2
tplot_options,'ythick',2
tplot_options,'labflag',-1


;First load "good" data
rbsp_load_efw_waveform, probe=probe, datatype='vsvy', coord = 'uvw',/noclean
get_data,'rbsp'+probe+'_efw_vsvy',times,v
efield = 10*[[v[*,0]-v[*,1]],[v[*,2]-v[*,3]],[v[*,4]-v[*,5]]]
store_data,'rbsp'+probe+'_efw_esvy',times,efield
copy_data,'rbsp'+probe+'_efw_esvy','rbsp'+probe+'_efw_esvy_orig'
split_vec,'rbsp'+probe+'_efw_esvy_orig',suffix='_'+['x','y','z']


;Now load adjusted data
rbsp_efw_create_esvy_uvw_from_vsvy,date,probe,bad_probe,method=1
copy_data,'rbsp'+probe+'_efw_esvy','rbsp'+probe+'_efw_esvy_adj'
split_vec,'rbsp'+probe+'_efw_esvy_adj',suffix='_'+['x','y','z']

;Now load adjusted data with method 2
rbsp_efw_create_esvy_uvw_from_vsvy,date,probe,bad_probe,method=2
copy_data,'rbsp'+probe+'_efw_esvy','rbsp'+probe+'_efw_esvy_adj2'
split_vec,'rbsp'+probe+'_efw_esvy_adj2',suffix='_'+['x','y','z']





store_data,'Eucomb',data=['rbsp'+probe+'_efw_esvy_orig_x','rbsp'+probe+'_efw_esvy_adj_x']
store_data,'Evcomb',data=['rbsp'+probe+'_efw_esvy_orig_y','rbsp'+probe+'_efw_esvy_adj_y']
store_data,'Eucomb2',data=['rbsp'+probe+'_efw_esvy_orig_x','rbsp'+probe+'_efw_esvy_adj2_x']
store_data,'Evcomb2',data=['rbsp'+probe+'_efw_esvy_orig_y','rbsp'+probe+'_efw_esvy_adj2_y']
store_data,'Eucomb3',data=['rbsp'+probe+'_efw_esvy_adj_x','rbsp'+probe+'_efw_esvy_adj2_x']
store_data,'Evcomb3',data=['rbsp'+probe+'_efw_esvy_adj_y','rbsp'+probe+'_efw_esvy_adj2_y']
options,'E?comb','colors',[0,250]
options,'E?comb?','colors',[0,250]


tplot,['Eucomb','Evcomb']
tplot,['Eucomb2','Evcomb2']
tplot,['Eucomb3','Evcomb3']

stop



rbsp_uvw_to_mgse,'a','rbsp'+probe+'_efw_esvy_orig'
split_vec,'rbsp'+probe+'_efw_esvy_orig_mgse'

rbsp_uvw_to_mgse,'a','rbsp'+probe+'_efw_esvy_adj'
split_vec,'rbsp'+probe+'_efw_esvy_adj_mgse'

rbsp_uvw_to_mgse,'a','rbsp'+probe+'_efw_esvy_adj2'
split_vec,'rbsp'+probe+'_efw_esvy_adj2_mgse'


store_data,'emgse_ycomb',data=['rbsp'+probe+'_efw_esvy_orig_mgse_y','rbsp'+probe+'_efw_esvy_adj_mgse_y']
store_data,'emgse_zcomb',data=['rbsp'+probe+'_efw_esvy_orig_mgse_z','rbsp'+probe+'_efw_esvy_adj_mgse_z']
store_data,'emgse_ycomb2',data=['rbsp'+probe+'_efw_esvy_orig_mgse_y','rbsp'+probe+'_efw_esvy_adj2_mgse_y']
store_data,'emgse_zcomb2',data=['rbsp'+probe+'_efw_esvy_orig_mgse_z','rbsp'+probe+'_efw_esvy_adj2_mgse_z']
store_data,'emgse_ycomb3',data=['rbsp'+probe+'_efw_esvy_adj_mgse_y','rbsp'+probe+'_efw_esvy_adj2_mgse_y']
store_data,'emgse_zcomb3',data=['rbsp'+probe+'_efw_esvy_adj_mgse_z','rbsp'+probe+'_efw_esvy_adj2_mgse_z']
options,'emgse_?comb','colors',[0,250]
options,'emgse_?comb?','colors',[0,250]

tplot,['emgse_ycomb','emgse_zcomb']
tplot,['emgse_ycomb2','emgse_zcomb2']
tplot,['emgse_ycomb3','emgse_zcomb3']


beep
beep
beep


;plot Ey_orig vs Ey_adj, etc.
t0z = date + '/09:00'
t1z = date + '/10:00'

time_clip,'rbsp'+probe+'_efw_esvy_adj_mgse_y',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_adj_mgse_y_clip'
time_clip,'rbsp'+probe+'_efw_esvy_adj_mgse_z',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_adj_mgse_z_clip'
time_clip,'rbsp'+probe+'_efw_esvy_adj2_mgse_y',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_adj2_mgse_y_clip'
time_clip,'rbsp'+probe+'_efw_esvy_adj2_mgse_z',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_adj2_mgse_z_clip'
time_clip,'rbsp'+probe+'_efw_esvy_orig_mgse_y',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_orig_mgse_y_clip'
time_clip,'rbsp'+probe+'_efw_esvy_orig_mgse_z',t0z,t1z,newname='rbsp'+probe+'_efw_esvy_orig_mgse_z_clip'


get_data,'rbsp'+probe+'_efw_esvy_orig_mgse_y_clip',data=y0c
get_data,'rbsp'+probe+'_efw_esvy_orig_mgse_z_clip',data=z0c
get_data,'rbsp'+probe+'_efw_esvy_adj_mgse_y_clip',data=y1c
get_data,'rbsp'+probe+'_efw_esvy_adj_mgse_z_clip',data=z1c
get_data,'rbsp'+probe+'_efw_esvy_adj2_mgse_y_clip',data=y2c
get_data,'rbsp'+probe+'_efw_esvy_adj2_mgse_z_clip',data=z2c


;Compare original data with adjusted data
!p.multi = [0,2,2]
;method1
plot,y0c.y,y1c.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method1',xrange=[-10,10],yrange=[-10,10]
oplot,[-10,10],[-10,10],color=250
plot,z0c.y,z1c.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method1',xrange=[-10,10],yrange=[-10,10]
oplot,[-10,10],[-10,10],color=250
;method2
plot,y0c.y,y2c.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method2',xrange=[-10,10],yrange=[-10,10]
oplot,[-10,10],[-10,10],color=250
plot,z0c.y,z2c.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method2',xrange=[-10,10],yrange=[-10,10]
oplot,[-10,10],[-10,10],color=250


;Plot detrended versions to ignore DC offset

tval = 2  ;min
rbsp_detrend,'rbsp'+probe+'_efw_esvy_orig_mgse_y_clip',60.*tval
rbsp_detrend,'rbsp'+probe+'_efw_esvy_orig_mgse_z_clip',60.*tval
rbsp_detrend,'rbsp'+probe+'_efw_esvy_adj_mgse_y_clip',60.*tval
rbsp_detrend,'rbsp'+probe+'_efw_esvy_adj_mgse_z_clip',60.*tval
rbsp_detrend,'rbsp'+probe+'_efw_esvy_adj2_mgse_y_clip',60.*tval
rbsp_detrend,'rbsp'+probe+'_efw_esvy_adj2_mgse_z_clip',60.*tval

get_data,'rbsp'+probe+'_efw_esvy_orig_mgse_y_clip_detrend',data=y0cd
get_data,'rbsp'+probe+'_efw_esvy_orig_mgse_z_clip_detrend',data=z0cd
get_data,'rbsp'+probe+'_efw_esvy_adj_mgse_y_clip_detrend',data=y1cd
get_data,'rbsp'+probe+'_efw_esvy_adj_mgse_z_clip_detrend',data=z1cd
get_data,'rbsp'+probe+'_efw_esvy_adj2_mgse_y_clip_detrend',data=y2cd
get_data,'rbsp'+probe+'_efw_esvy_adj2_mgse_z_clip_detrend',data=z2cd


;Compare original detrended data with adjusted data
!p.multi = [0,2,2]
;method1
plot,y0cd.y,y1cd.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z0cd.y,z1cd.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
;method2
plot,y0cd.y,y2cd.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z0cd.y,z2cd.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250

stop


;Compare method 1 and method 2
!p.multi = [0,2,2]
;method1
plot,y1c.y,y2c.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z1c.y,z2c.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
;method2
plot,y1c.y,y2c.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z1c.y,z2c.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250


;Compare detrended method 1 and method 2
!p.multi = [0,2,2]
;method1
plot,y1cd.y,y2cd.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z1cd.y,z2cd.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method1';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
;method2
plot,y1cd.y,y2cd.y,psym=4,xtitle='mV/m (MGSEy)',ytitle='mV/m (MGSEy fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250
plot,z1cd.y,z2cd.y,psym=4,xtitle='mV/m (MGSEz)',ytitle='mV/m (MGSEz fixed)',title='method2';,xrange=[-3,3],yrange=[-3,3]
oplot,[-10,10],[-10,10],color=250

stop


end
