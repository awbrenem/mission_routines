;Testing cribsheet

rbsp_efw_init

path = '/Users/aaronbreneman/Desktop/code/Aaron/RBSP/TDAS_trunk_svn/general/missions/rbsp/efw/examples/'
tplot_restore,filenames=path+'chastontest.tplot'
restore,path + 'chastontest.idl'



get_data,varE_s4+'_tmp',data=d

;interpolate time steps to constant timegrid.

ttmp0 = d.x[0]
ttmp1 = d.x[n_elements(d.x)-1]

;find average sample rate.
sr = rbsp_sample_rate(d.x,OUT_MED_AVG=medavg)
sr = medavg[0]
stepsz = 1./sr

;smoothtimes = stepsz*dindgen(n_elements(d.x))/(n_elements(d.x)-1) + ttmp0
smoothtimes = stepsz*dindgen(n_elements(d.x)) + ttmp0

;Called from rbsp_efw_burst_fa_rotate_crib
;twavpol_modified,varE_s4+'_tmp',prefix='tmp',nopfft=nopfft,steplength=steplength



;From within twavpol_modified
t = smoothtimes
d1 = d.y[*,0]
d2 = d.y[*,1]
d3 = d.y[*,2]

  wavpol, t, d1,d2,d3, timeline, freqline,$
   powspec, degpol, waveangle, elliptict, helict, pspec3, nopfft=nopfft,steplength=steplength;,bin_freq=1;bin_freq

print,total(powspec,/nan)
print,total(elliptict,/nan)



store_data,'powtmp',data={x:timeline,y:powspec,v:freqline}
options,'powtmp','spec',1
zlim,'powtmp',1d-9,1d-6,1
tplot,['powtmp',varE_s4+'_tmp']


store_data,'elltmp',data={x:timeline,y:elliptict,v:freqline}
options,'elltmp','spec',1
zlim,'elltmp',-0.5,0.5,0
tplot,['powtmp','elltmp',varE_s4+'_tmp']



nopoints = n_elements(d.x)
nosteps=(nopoints-nopfft)/steplength
