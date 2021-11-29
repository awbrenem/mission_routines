;simple loop for creating L3, L4 files


d0 = time_double('2017-02-14')
d1 = time_double('2017-04-01')

d0 = time_double('2017-05-15')
d1 = time_double('2017-07-01')

d0 = time_double('2017-08-13')
d1 = time_double('2017-10-01')

d0 = time_double('2017-11-14')
d1 = time_double('2018-01-01')


probe = 'b'
ndays = (d1-d0)/86400
for i=0,ndays-1 do begin $
    dtmp = time_string(d0 + 86400*i) & $
    dtmp = strmid(dtmp,0,10) & $
    rbsp_efw_make_l3,probe,dtmp,boom_pair='34',/testing & $
    store_data,'*',/del


