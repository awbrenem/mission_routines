;Read GOES csv mag files from Sam Califf.
;He threw these together for GOES14 for me on 2019-05-09

;t_datenum: Matlab datenum in days, t = 0 at 0000-00-00, t = 730486 at 2000-01-01 00:00:00
;b_nT_epn: EPN magnetic field from the outboard mag


file = 'goes_14_mag_2014-01-11.csv'
path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events_campaign2/Jan11/'


ft=[5,4,4,4,7,7,7,7,7,7,7,7,7,7,7,4,4,4,7,7,7,4,4,4,7,7,7,4,7,4,4,4,7,7,7,4,4,4]
fn=['FIELD01','FIELD02','FIELD03','FIELD04','FIELD05','FIELD06','FIELD07','FIELD08','FIELD09','FIELD10','FIELD11','FIELD12','FIELD13','FIELD14','FIELD15','FIELD16','FIELD17','FIELD18','FIELD19','FIELD20','FIELD21','FIELD22','FIELD23','FIELD24','FIELD25','FIELD26','FIELD27','FIELD28','FIELD29','FIELD30','FIELD31','FIELD32','FIELD33','FIELD34','FIELD35','FIELD36','FIELD37','FIELD38']
fl = [0,18,35,52,69,73,77,81,85,89,93,97,101,105,109,113,128,142,157,161,165,169,183,198,212,216,220,224,239,243,257,271,285,289,293,297,311,325]
fg = indgen(38)

s = {version: 1.,$
    datastart: 2L,$
    delimiter:44B,$
    MISSINGVALUE:!values.f_nan,$
    COMMENTSYMBOL:'',$
    fieldcount:38L,$
    fieldtypes:ft,$
    fieldnames:fn,$
    fieldlocations:fl,$
    fieldgroups:fg}


x = read_ascii(path+file,template=s)


t_datenumT = x.FIELD01
b1t = x.FIELD16
b2t = x.FIELD17
b3t = x.FIELD18


;t = 730486 at 2000-01-01 00:00:00
t0 = 730486d

;G14 seconds relative to first listed time
tsec = (t_datenumT mod 1d) *86400d

;Absolute first time listed
nsecs = (t_datenumT[0] - t0)*86400d
tinit = time_double('2000-01-01/00:00:00') + nsecs
print,time_string(tinit[0])


newtimes = tinit + tsec

bmag = sqrt(b1t^2 + b2t^2 + b3t^2)

store_data,'G14_bfield_EPN',newtimes,[[b1t],[b2t],[b3t]]
store_data,'G14_bmag',newtimes,bmag
tplot,['G14_bfield_EPN','G14_bmag']


tplot_save,['G14_bfield_EPN','G14_bmag'],filename='~/Desktop/G14_bfield'

end
