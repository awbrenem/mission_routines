;Darius is trying to find the angular separation between two footpoints defined by L, MLT. 
;SDK is telling him the separation is 12.2970 deg and he's calculating a separation of 8.0233 deg. 
;L1 = 1.1252
;L2 = 1.2570 
;MLT = 7.0970




;Take the lat, long, alt files from Xinlin's group and map them to magnetic coordinates for conjunction comparisons. 

rbsp_efw_init



ft = [3,7,3,7,4,4,4]
fn = ['day','month','year','time','lat','lon','alt']
fl = [1,3,7,12,31,42,54]
fg = indgen(7)

template = {version: 1.,$
            datastart:6L,$
            delimiter:32B,$
            missingvalue:!values.f_nan,$
            commentsymbol:'',$
            fieldcount:7L,$
            fieldtypes:ft,$
            fieldnames:fn,$
            fieldlocations:fl,$
            fieldgroups:fg}

path = '/Users/aaronbreneman/Desktop/Research/OTHER/proposals/2020_CubeSat_Xinlin/'
;plds = ['C1','C2','C3','C4','G1','G2','G3']
plds = ['C1','G1']
suffix = '_data.txt'

for i=0,n_elements(plds)-1 do begin 
    vals = read_ascii(path+plds[i]+suffix,template=template)

    monthf = replicate('',n_elements(vals.day))
    goo = where(vals.month eq 'Jan') & if goo[0] ne -1 then monthf[goo] = '01'
    goo = where(vals.month eq 'Feb') & if goo[0] ne -1 then monthf[goo] = '02'
    goo = where(vals.month eq 'Mar') & if goo[0] ne -1 then monthf[goo] = '03'
    goo = where(vals.month eq 'Apr') & if goo[0] ne -1 then monthf[goo] = '04'
    goo = where(vals.month eq 'May') & if goo[0] ne -1 then monthf[goo] = '05'
    goo = where(vals.month eq 'Jun') & if goo[0] ne -1 then monthf[goo] = '06'
    goo = where(vals.month eq 'Jul') & if goo[0] ne -1 then monthf[goo] = '07'
    goo = where(vals.month eq 'Aug') & if goo[0] ne -1 then monthf[goo] = '08'
    goo = where(vals.month eq 'Sep') & if goo[0] ne -1 then monthf[goo] = '08'
    goo = where(vals.month eq 'Oct') & if goo[0] ne -1 then monthf[goo] = '10'
    goo = where(vals.month eq 'Nov') & if goo[0] ne -1 then monthf[goo] = '11'
    goo = where(vals.month eq 'Dec') & if goo[0] ne -1 then monthf[goo] = '12'

    dayf = replicate('',n_elements(vals.day))
    goo = where(vals.day lt 10.)
    if goo[0] ne -1 then dayf[goo] = '0' + strtrim(vals.day[goo],2)
    goo = where(vals.day ge 10.)
    if goo[0] ne -1 then dayf[goo] = strtrim(vals.day[goo],2)

    yearf = strtrim(vals.year,2)

    dtime = yearf + '-' + monthf + '-' + dayf + '/' + vals.time
    dtimed = time_double(dtime)


    store_data,plds[i]+'_geolat',dtimed,vals.lat 
    store_data,plds[i]+'_geolon',dtimed,vals.lon
    store_data,plds[i]+'_alt',dtimed,vals.alt + 6370.


    ;tplot,['geolat','geolon','alt']


    ;;--------------------------------------------
    ;;Compare plot to what Darius sent me (looks very close)
    ;t0z = time_double('2020-08-16')
    ;t1z = time_double('2020-08-16/03:00')
    ;geolat = tsample('G3_geolat',[t0z,t1z],times=tms)
    ;geolon = tsample('G3_geolon',[t0z,t1z],times=tms)
    ;plot_geo_coord_on_map,geolon,geolat
    ;;--------------------------------------------



;    get_data,plds[i]+'_geolat',data=geolat
 ;   get_data,plds[i]+'_geolon',data=geolon
 ;   get_data,plds[i]+'_alt',data=geoalt
;
;    times = geolat.x

    ;Scale down the times
 ;   times2 = times 
    alts6370 = vals.alt + 6370.
 ;   lats2 = geolat.y 
 ;   lons2 = geolon.y


    xgeo = alts6370*cos(!dtor*vals.lats)*cos(!dtor*vals.lons)
    ygeo = alts6370*cos(!dtor*vals.lats)*sin(!dtor*vals.lons)
    zgeo = alts6370*sin(!dtor*vals.lats)
    store_data,plds[i]+'_geo',dtimed,[[xgeo],[ygeo],[zgeo]]


    ;Calculate MLT
    cotrans,plds[i]+'_geo',plds[i]+'_gei',/geo2gei
    cotrans,plds[i]+'_gei',plds[i]+'_gse',/gei2gse
    cotrans,plds[i]+'_gse',plds[i]+'_gsm',/gse2gsm
    cotrans,plds[i]+'_gsm',plds[i]+'_sm',/gsm2sm


endfor





;Now that data is loaded, let's do the mapping
for i=0,n_elements(plds)-1 do begin 

    ;Reduce to times of interest
    t0z = time_double('2020-08-16/00:00')
    t1z = time_double('2020-08-18/00:00')

    get_data,plds[i]+'_gsm',data=dd
    yv = tsample(plds[i]+'_gsm',[t0z,t1z],times=tms)
    store_data,plds[i]+'_gsm2',tms,yv



    ;Make sure R0 is ABOVE the Earth's surface. Takes too long to map to the surface
    ;Note that R0=50 above surface and rlim=10 lead to glitches
    R0 = 6370. + 25.
    rlim = 20.*6370.
    duration = t1z - t0z
    start_time = time_string(t0z)
    kp = 2.

    model = 't89'
;    model = 'none'  ;IGRF
    aaron_map_with_tsy,model,start_time,duration,plds[i],plds[i]+'_gsm2',Kp,R0=R0,rlim=rlim


    get_data,plds[i]+'_out_iono_foot_north_glat_glon',data=tmp
    geolat = tmp.y[*,0] & geolon = tmp.y[*,1]
    store_data,plds[i]+'_out_iono_foot_north_geolat',tmp.x,geolat
    store_data,plds[i]+'_out_iono_foot_north_geolon',tmp.x,geolon

    get_data,plds[i]+'_out_iono_foot_south_glat_glon',data=tmp
    geolat = tmp.y[*,0] & geolon = tmp.y[*,1]
    store_data,plds[i]+'_out_iono_foot_south_geolat',tmp.x,geolat
    store_data,plds[i]+'_out_iono_foot_south_geolon',tmp.x,geolon


    ;Mapped quantities have a lower cadence. Interpolate the input times to the 
    ;mapped times. 

    tinterpol_mxn,plds[i]+'_geolat',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_geolon',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_alt',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_gse',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_geo',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_gei',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_gsm',tmp.x,/overwrite,/quadratic
    tinterpol_mxn,plds[i]+'_sm',tmp.x,/overwrite,/quadratic





    ;For the C1-C4 spacecraft, remove the N or S mapped points based on which hemisphere the sc is in. 
    
    if plds[i] eq 'C1' or plds[i] eq 'C2' or plds[i] eq 'C3' or plds[i] eq 'C4' then begin
        get_data,plds[i]+'_geolat',data=tmp 
        times = tmp.x
        cxgeolat = tmp.y
        get_data,plds[i]+'_geolon',data=tmp
        cxgeolon = tmp.y

        get_data,plds[i]+'_out_iono_foot_north_geolat',data=tmp
        cxgeolatN = tmp.y
        get_data,plds[i]+'_out_iono_foot_north_geolon',data=tmp
        cxgeolonN = tmp.y
        get_data,plds[i]+'_out_iono_foot_south_geolat',data=tmp
        cxgeolatS = tmp.y
        get_data,plds[i]+'_out_iono_foot_south_geolon',data=tmp
        cxgeolonS = tmp.y


        goo = where(cxgeolat lt 0.)
        if goo[0] ne -1 then cxgeolatN[goo] = !values.f_nan
        if goo[0] ne -1 then cxgeolonN[goo] = !values.f_nan
        goo = where(cxgeolat ge 0.)
        if goo[0] ne -1 then cxgeolatS[goo] = !values.f_nan
        if goo[0] ne -1 then cxgeolonS[goo] = !values.f_nan

        store_data,plds[i]+'_out_iono_foot_north_geolat',times,cxgeolatN
        store_data,plds[i]+'_out_iono_foot_north_geolon',times,cxgeolonN
        store_data,plds[i]+'_out_iono_foot_south_geolat',times,cxgeolatS
        store_data,plds[i]+'_out_iono_foot_south_geolon',times,cxgeolonS

    endif
endfor


stop


;-----------------------------------------------------------
;Remove mapped MLT values for low alt sats when Lshell gets beyond a certain threshold. 
;This is clearly leading to anomalous MLT values. 
;-----------------------------------------------------------

lmax = 15.
get_data,'C1!CL-shell-t89',data=dd 
goo = where(dd.y ge lmax)
if goo[0] ne -1 then dd.y[goo] = !values.f_nan
store_data,'C1!CL-shell-t89',data=dd 

get_data,'C1!Cnorth-foot-MLT!Ct89',data=dd
if goo[0] ne -1 then dd.y[goo] = !values.f_nan
store_data,'C1!Cnorth-foot-MLT!Ct89',data=dd 
get_data,'C1!Csouth-foot-MLT!Ct89',data=dd
if goo[0] ne -1 then dd.y[goo] = !values.f_nan
store_data,'C1!Csouth-foot-MLT!Ct89',data=dd 






;--------------------------------------------------------
;TEST MLT VALUES **********************
;---------------------------------------------------------

get_data,'C1_gse',data=tmpgse
postimes_gse=tmpgse.x


; MLT of SC location (not foot point)
angle_tmp = atan(tmpgse.y[*,1],tmpgse.y[*,0])/!dtor
goo = where(angle_tmp lt 0.)
if goo[0] ne -1 then angle_tmp[goo] = 360. - abs(angle_tmp[goo])
mloctime = angle_tmp * 12/180. + 12.
goo = where(mloctime ge 24.)
if goo[0] ne -1 then mloctime[goo] = mloctime[goo] - 24
store_data,'C1_MLT',data={x:postimes_gse,y:mloctime}
tinterpol_mxn,'C1_MLT',times,/overwrite,/quadratic

tplot,['C1!CL-shell-t89','C1_MLT','C1!Cnorth-foot-MLT!Ct89','C1!Csouth-foot-MLT!Ct89']



;---------------------------------------------------------------------------
;Intelligently take differences (MLT values, mostly)
;---------------------------------------------------------------------------



dif_data,'C1!CL-shell-t89','G1!CL-shell-t89',newname='C1G1_Ldiff'
dif_data,'C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89',newname='C1G1N_MLTdiff'
dif_data,'C1!Csouth-foot-MLT!Ct89','G1!Csouth-foot-MLT!Ct89',newname='C1G1S_MLTdiff'
dif_data,'C1_MLT','G1!Cnorth-foot-MLT!Ct89',newname='C1G1NE_MLTdiff'
dif_data,'C1_MLT','G1!Csouth-foot-MLT!Ct89',newname='C1G1SE_MLTdiff'
dif_data,'C1_out_iono_foot_north_gse','G1_out_iono_foot_north_gse',newname='C1G1N_GSEdiff'
dif_data,'C1_out_iono_foot_south_gse','G1_out_iono_foot_south_gse',newname='C1G1S_GSEdiff'



store_data,'C1G1_Lcomb',data=['C1!CL-shell-t89','G1!CL-shell-t89']
options,'C1G1_Lcomb','colors',[0,250] & options,'C1G1_Lcomb','psym',-4

get_data,'C1G1_Ldiff',data=tmpL
store_data,'C1G1_Ldiff_abs',tmpL.x,abs(tmpL.y)


;-------------------------
;Create absolute difference values (ignores switch at 24 MLT)
;NOTE: delta-MLT values should range from -12 to 12. 



get_data,'C1G1N_MLTdiff',data=tmp
goo = where(tmp.y ge 12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] - 24.
goo = where(tmp.y lt -12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] + 24.
store_data,'C1G1N_MLTdiff',data=tmp
store_data,'C1G1N_MLTdiff_abs',tmp.x,abs(tmp.y)

;tplot,['C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89','C1G1N_MLTdiff']

;TESTING
;    options,['C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89','C1G1N_MLTdiff','C1G1N_MLTdiff_tmp'],'psym',-4
;    store_data,'combtmp',data=['C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89']
;    options,'combtmp','colors',[0,250]
;tplot,['C1!CL-shell-t89','combtmp','C1G1N_MLTdiff_tmp','C1G1N_MLTdiff']

get_data,'C1G1S_MLTdiff',data=tmp
goo = where(tmp.y ge 12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] - 24.
goo = where(tmp.y lt -12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] + 24.
store_data,'C1G1S_MLTdiff',data=tmp
store_data,'C1G1S_MLTdiff_abs',tmp.x,abs(tmp.y)

get_data,'C1G1NE_MLTdiff',data=tmp
goo = where(tmp.y ge 12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] - 24.
goo = where(tmp.y lt -12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] + 24.
store_data,'C1G1NE_MLTdiff',data=tmp
store_data,'C1G1NE_MLTdiff_abs',tmp.x,abs(tmp.y)

;   TESTING (DATA ARE GOOD!!!)
;    options,['C1_MLT','G1!Cnorth-foot-MLT!Ct89','C1G1NE_MLTdiff','C1G1NE_MLTdiff_tmp'],'psym',-4
;    store_data,'combtmp',data=['C1_MLT','G1!Cnorth-foot-MLT!Ct89']
;    options,'combtmp','colors',[0,250]
;    tplot,['combtmp','C1G1NE_MLTdiff_tmp','C1G1NE_MLTdiff']

get_data,'C1G1SE_MLTdiff',data=tmp
goo = where(tmp.y ge 12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] - 24.
goo = where(tmp.y lt -12.)
if goo[0] ne -1 then tmp.y[goo] = tmp.y[goo] + 24.
store_data,'C1G1SE_MLTdiff',data=tmp
store_data,'C1G1SE_MLTdiff_abs',tmp.x,abs(tmp.y)


tplot,['C1G1N_MLTdiff','C1G1S_MLTdiff','C1G1NE_MLTdiff','C1G1SE_MLTdiff']
tplot,['C1G1N_MLTdiff','C1G1S_MLTdiff','C1G1NE_MLTdiff','C1G1SE_MLTdiff']+'_abs'





;-------------------------------------------------------------------------------
;FIND VARIOUS TYPES OF CONJUNCTIONS
;--ABSOLUTE CONJUNCTION 
;--DRIFT CONJUNCTION, ETC.
;-------------------------------------------------------------------------------



    ;------------------------------------------------------------------------------
    ;Find absolute L, MLT conjunction (not drift conjunction)
    ;------------------------------------------------------------------------------

    deltaMLTmax = 1.
    deltaLmax = 1.


    ;let's reduce to the deltaMLTmax and deltaLmax times 
    get_data,'C1G1N_MLTdiff_abs',data=tmpMLT
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_geolat',data=tmpgeo

    deltaL_abs_line = replicate(deltaLmax,n_elements(tmpL.x))
    store_data,'deltaL_abs_line',tmpL.x,deltaL_abs_line & options,'deltaL_abs_line','thick',1 & options,'deltaL_abs_line','color',50
    store_data,'C1G1_Ldiff_abs_comb',data=['C1G1_Ldiff_abs','deltaL_abs_line']

    deltaMLT_abs_line = replicate(deltaMLTmax,n_elements(tmpL.x)) & options,'deltaMLT_abs_line','psym',-4
    store_data,'deltaMLT_abs_line',tmpL.x,deltaL_abs_line & options,'deltaMLT_abs_line','thick',1 & options,'deltaMLT_abs_line','color',50
    store_data,'C1G1N_MLTdiff_abs_comb',data=['C1G1N_MLTdiff_abs','deltaMLT_abs_line']



    goo1 = where(tmpMLT.y le deltaMLTmax)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(tmpL.y le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpgeo.y ge 0.)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    comb = tmpy1 + tmpy2 + tmpy3
    goo = where(comb eq 3)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1N_absolute_conjunction',tmpL.x,val & ylim,'C1G1N_absolute_conjunction',0,2
    options,'C1G1N_absolute_conjunction','symsize',2
    options,'C1G1N_absolute_conjunction','thick',2
    options,'C1G1N_absolute_conjunction','psym',-4
    options,'C1G1N_absolute_conjunction','panel_size',0.5
    options,['C1_geolat','C1G1N_absolute_conjunction','C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89','C1G1N_MLTdiff_abs','C1G1_Ldiff_abs_comb'],'psym',-4
    tplot,['C1_geolat','C1G1N_absolute_conjunction','C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89','C1G1N_MLTdiff_abs','C1G1_Ldiff_abs_comb']



    ;-----
    get_data,'C1G1S_MLTdiff_abs',data=tmpMLT
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_geolat',data=tmpgeo


    store_data,'C1G1S_MLTdiff_abs_comb',data=['C1G1S_MLTdiff_abs','deltaMLT_abs_line']

    goo1 = where(tmpMLT.y le deltaMLTmax)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(tmpL.y le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpgeo.y lt 0.)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    comb = tmpy1 + tmpy2 + tmpy3
    goo = where(comb eq 3)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1S_absolute_conjunction',tmpL.x,val & ylim,'C1G1S_absolute_conjunction',0,2
    options,'C1G1S_absolute_conjunction','symsize',2
    options,'C1G1S_absolute_conjunction','thick',2
    options,'C1G1S_absolute_conjunction','psym',-4
    options,'C1G1S_absolute_conjunction','panel_size',0.5

    options,['C1_geolat','C1G1S_absolute_conjunction','C1!Csouth-foot-MLT!Ct89','G1!Csouth-foot-MLT!Ct89','C1G1S_MLTdiff_abs','C1G1_Ldiff_abs_comb'],'psym',-4
    tplot,['C1_geolat','C1G1S_absolute_conjunction','C1!Csouth-foot-MLT!Ct89','G1!Csouth-foot-MLT!Ct89','C1G1S_MLTdiff_abs','C1G1_Ldiff_abs_comb']

    ;------
    ;SAME CALCULATION AS ABOVE BUT USING C1_MLT INSTEAD OF THE FOOTPOINT VERSION

    get_data,'C1G1NE_MLTdiff_abs',data=tmpMLT
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_geolat',data=tmpgeo

    store_data,'C1G1NE_MLTdiff_abs_comb',data=['C1G1NE_MLTdiff_abs','deltaMLT_abs_line']

    goo1 = where(tmpMLT.y le deltaMLTmax)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(tmpL.y le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpgeo.y ge 0.)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    comb = tmpy1 + tmpy2 + tmpy3
    goo = where(comb eq 3)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1NE_absolute_conjunction',tmpL.x,val & ylim,'C1G1NE_absolute_conjunction',0,2
    options,'C1G1NE_absolute_conjunction','symsize',2
    options,'C1G1NE_absolute_conjunction','thick',2
    options,'C1G1NE_absolute_conjunction','psym',-4
    options,'C1G1NE_absolute_conjunction','panel_size',0.5

    tplot,['C1_geolat','C1G1NE_absolute_conjunction','C1_MLT','G1!Cnorth-foot-MLT!Ct89','C1G1NE_MLTdiff_abs','C1G1_Ldiff_abs_comb']


    ;-------
    get_data,'C1G1SE_MLTdiff_abs',data=tmpMLT
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_geolat',data=tmpgeo

    store_data,'C1G1SE_MLTdiff_abs_comb',data=['C1G1SE_MLTdiff_abs','deltaMLT_abs_line']

    goo1 = where(tmpMLT.y le deltaMLTmax)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(tmpL.y le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpgeo.y lt 0.)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    comb = tmpy1 + tmpy2 + tmpy3
    goo = where(comb eq 3)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1SE_absolute_conjunction',tmpL.x,val & ylim,'C1G1SE_absolute_conjunction',0,2
    options,'C1G1SE_absolute_conjunction','symsize',2
    options,'C1G1SE_absolute_conjunction','thick',2
    options,'C1G1SE_absolute_conjunction','psym',-4
    options,'C1G1SE_absolute_conjunction','panel_size',0.5
    tplot,['C1_geolat','C1G1SE_absolute_conjunction','C1_MLT','G1!Csouth-foot-MLT!Ct89','C1G1SE_MLTdiff_abs','C1G1_Ldiff_abs_comb']

stop





    ;--------------------------------------------------------------------
    ;Find drift conjunctions (C payload required to be East of G payload)
    ;--------------------------------------------------------------------

    deltaMLTmax_drift = 4.

    deltaMLT_drift_line = replicate(deltaMLTmax_drift,n_elements(tmpL.x))
    store_data,'deltaMLT_drift_line',times,deltaMLT_drift_line
    options,'deltaMLT_drift_line','thick',2

    ;Now determine where C1 is East of G1
    driftNgood = replicate(1.,n_elements(tmp.x))
    get_data,'C1G1N_MLTdiff',data=tmp
    goo = where(tmp.y lt 0.)  ;C East of G (bad)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    goo = where(tmp.y gt 12.) ;C East of G (bad)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    goo = where(finite(tmp.y) eq 0.)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    store_data,'C1_east_of_G1N',tmp.x,driftNgood & ylim,'C1_east_of_G1N',0,2
    options,'C1_east_of_G1N','psym',-4 & options,'C1_east_of_G1N','panel_size',0.5
    store_data,'C1G1N_MLT_comb',data=['C1!Cnorth-foot-MLT!Ct89','G1!Cnorth-foot-MLT!Ct89']
    options,'C1G1N_MLT_comb','colors',[0,250]
    tplot,['C1_east_of_G1N','C1G1N_MLT_comb','C1G1N_MLTdiff']


    driftSgood = replicate(1.,n_elements(tmp.x))
    get_data,'C1G1S_MLTdiff',data=tmp
    goo = where(tmp.y lt 0.)  ;C East of G (bad)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    goo = where(tmp.y gt 12.) ;C East of G (bad)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    goo = where(finite(tmp.y) eq 0.)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    store_data,'C1_east_of_G1S',tmp.x,driftSgood & ylim,'C1_east_of_G1S',0,2
    options,'C1_east_of_G1S','psym',-4 & options,'C1_east_of_G1S','panel_size',0.5
    store_data,'C1G1S_MLT_comb',data=['C1!Csouth-foot-MLT!Ct89','G1!Csouth-foot-MLT!Ct89']
    options,'C1G1S_MLT_comb','colors',[0,250]
    tplot,['C1_east_of_G1S','C1G1S_MLT_comb','C1G1S_MLTdiff']


    driftNgood = replicate(1.,n_elements(tmp.x))
    get_data,'C1G1NE_MLTdiff',data=tmp
    goo = where(tmp.y lt 0.)  ;C East of G (bad)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    goo = where(tmp.y gt 12.) ;C East of G (bad)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    goo = where(finite(tmp.y) eq 0.)
    if goo[0] ne -1 then driftNgood[goo] = !values.f_nan
    store_data,'C1_east_of_G1NE',tmp.x,driftNgood & ylim,'C1_east_of_G1NE',0,2
    options,'C1_east_of_G1NE','psym',-4 & options,'C1_east_of_G1NE','panel_size',0.5
    store_data,'C1G1NE_MLT_comb',data=['C1_MLT','G1!Cnorth-foot-MLT!Ct89']
    options,'C1G1NE_MLT_comb','colors',[0,250]
    tplot,['C1_east_of_G1NE','C1G1NE_MLT_comb','C1G1NE_MLTdiff']


    driftSgood = replicate(1.,n_elements(tmp.x))
    get_data,'C1G1SE_MLTdiff',data=tmp
    goo = where(tmp.y lt 0.)  ;C East of G (bad)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    goo = where(tmp.y gt 12.) ;C East of G (bad)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    goo = where(finite(tmp.y) eq 0.)
    if goo[0] ne -1 then driftSgood[goo] = !values.f_nan
    store_data,'C1_east_of_G1SE',tmp.x,driftSgood & ylim,'C1_east_of_G1SE',0,2
    options,'C1_east_of_G1SE','psym',-4 & options,'C1_east_of_G1SE','panel_size',0.5
    store_data,'C1G1SE_MLT_comb',data=['C1_MLT','G1!Csouth-foot-MLT!Ct89']
    options,'C1G1SE_MLT_comb','colors',[0,250]
    tplot,['C1_east_of_G1SE','C1G1SE_MLT_comb','C1G1SE_MLTdiff']







    ;Now that we have all the times when the C payload is East of the G payload, let's reduce to the 
    ;deltaMLTmax and deltaLmax times 
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_east_of_G1N',data=tmpEast
    get_data,'C1G1N_MLTdiff',data=tmpMLT
    get_data,'C1_geolat',data=tmpgeo

    store_data,'C1G1N_MLTdiff_comb',data=['C1G1N_MLTdiff','deltaMLT_drift_line']
    options,'C1G1N_MLTdiff_comb','thick',2
    options,'C1G1N_MLTdiff_comb','colors',[0,50]

    goo1 = where(tmpMLT.y le deltaMLTmax_drift)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(abs(tmpL.y) le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpEast.y eq 1)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpEast.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    goo4 = where(tmpgeo.y ge 0.)
    tmpy4 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo4[0] ne -1 then tmpy4[goo4] = 1


    comb = tmpy1 + tmpy2 + tmpy3 + tmpy4
    goo = where(comb eq 4)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1N_drift_conjunction',tmpL.x,val & ylim,'C1G1N_drift_conjunction',0,2
    options,'C1G1N_drift_conjunction','symsize',2
    options,'C1G1N_drift_conjunction','thick',2
    options,'C1G1N_drift_conjunction','psym',-4
    options,'C1G1N_drift_conjunction','panel_size',0.5

    options,['C1_geolat','C1G1N_drift_conjunction','C1G1N_MLT_comb','C1_east_of_G1N','C1G1N_MLTdiff_comb','C1G1_Ldiff_abs_comb'],'psym',-4
    tplot,['C1_geolat','C1G1N_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1N_MLT_comb','C1G1N_MLTdiff_comb','C1_east_of_G1N']





    ;-----------
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_east_of_G1S',data=tmpEast
    get_data,'C1G1S_MLTdiff',data=tmpMLT
    get_data,'C1_geolat',data=tmpgeo

    store_data,'C1G1S_MLTdiff_comb',data=['C1G1S_MLTdiff','deltaMLT_drift_line']
    options,'C1G1S_MLTdiff_comb','thick',2
    options,'C1G1S_MLTdiff_comb','colors',[0,50]


    goo1 = where(tmpMLT.y le deltaMLTmax_drift)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(abs(tmpL.y) le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpEast.y eq 1)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpEast.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    goo4 = where(tmpgeo.y lt 0.)
    tmpy4 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo4[0] ne -1 then tmpy4[goo4] = 1

    comb = tmpy1 + tmpy2 + tmpy3 + tmpy4
    goo = where(comb eq 4)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1S_drift_conjunction',tmpL.x,val & ylim,'C1G1S_drift_conjunction',0,2
    options,'C1G1S_drift_conjunction','symsize',2
    options,'C1G1S_drift_conjunction','thick',2
    options,'C1G1S_drift_conjunction','psym',-4
    options,'C1G1S_drift_conjunction','panel_size',0.5

    options,['C1_geolat','C1G1S_drift_conjunction','C1G1_Ldiff_abs_comb','C1G1S_MLT_comb','C1_east_of_G1S','C1G1S_MLTdiff_comb'],'psym',-4
    tplot,['C1_geolat','C1G1S_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1S_MLT_comb','C1G1S_MLTdiff_comb','C1_east_of_G1S']



    ;---------
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_east_of_G1NE',data=tmpEast
    get_data,'C1G1NE_MLTdiff',data=tmpMLT
    get_data,'C1_geolat',data=tmpgeo

    store_data,'C1G1NE_MLTdiff_comb',data=['C1G1NE_MLTdiff','deltaMLT_drift_line']
    options,'C1G1NE_MLTdiff_comb','thick',2
    options,'C1G1NE_MLTdiff_comb','colors',[0,50]

    goo1 = where(tmpMLT.y le deltaMLTmax_drift)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(abs(tmpL.y) le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpEast.y eq 1)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpEast.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    goo4 = where(tmpgeo.y ge 0.)
    tmpy4 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo4[0] ne -1 then tmpy4[goo4] = 1


    comb = tmpy1 + tmpy2 + tmpy3 + tmpy4
    goo = where(comb eq 4)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1NE_drift_conjunction',tmpL.x,val & ylim,'C1G1NE_drift_conjunction',0,2
    options,'C1G1NE_drift_conjunction','symsize',2
    options,'C1G1NE_drift_conjunction','thick',2
    options,'C1G1NE_drift_conjunction','psym',-4
    options,'C1G1NE_drift_conjunction','panel_size',0.5

    options,['C1_geolat','C1G1NE_drift_conjunction','C1G1NE_MLT_comb','C1_east_of_G1NE','C1G1NE_MLTdiff_comb','C1G1_Ldiff_abs_comb'],'psym',-4
    tplot,['C1_geolat','C1G1NE_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1NE_MLT_comb','C1G1NE_MLTdiff_comb','C1_east_of_G1NE']



    ;------
    get_data,'C1G1_Ldiff_abs',data=tmpL
    get_data,'C1_east_of_G1SE',data=tmpEast
    get_data,'C1G1SE_MLTdiff',data=tmpMLT
    get_data,'C1_geolat',data=tmpgeo
   
    store_data,'C1G1SE_MLTdiff_comb',data=['C1G1SE_MLTdiff','deltaMLT_drift_line']
    options,'C1G1SE_MLTdiff_comb','thick',2
    options,'C1G1SE_MLTdiff_comb','colors',[0,50]


    goo1 = where(tmpMLT.y le deltaMLTmax_drift)
    tmpy1 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo1[0] ne -1 then tmpy1[goo1] = 1

    goo2 = where(abs(tmpL.y) le deltaLmax)
    tmpy2 = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo2[0] ne -1 then tmpy2[goo2] = 1

    goo3 = where(tmpEast.y eq 1)
    tmpy3 = replicate(!values.f_nan,n_elements(tmpEast.x))
    if goo3[0] ne -1 then tmpy3[goo3] = 1

    goo4 = where(tmpgeo.y lt 0.)
    tmpy4 = replicate(!values.f_nan,n_elements(tmpgeo.x))
    if goo4[0] ne -1 then tmpy4[goo4] = 1

    comb = tmpy1 + tmpy2 + tmpy3 + tmpy4
    goo = where(comb eq 4)
    val = replicate(!values.f_nan,n_elements(tmpL.x))
    if goo[0] ne -1 then val[goo] = 1
    store_data,'C1G1SE_drift_conjunction',tmpL.x,val & ylim,'C1G1SE_drift_conjunction',0,2
    options,'C1G1SE_drift_conjunction','symsize',2
    options,'C1G1SE_drift_conjunction','thick',2
    options,'C1G1SE_drift_conjunction','psym',-4
    options,'C1G1SE_drift_conjunction','panel_size',0.5

    options,['C1_geolat','C1G1SE_drift_conjunction','C1G1SE_MLT_comb','C1_east_of_G1SE','C1G1SE_MLTdiff_comb','C1G1_Ldiff_abs_comb'],'psym',-4
    tplot,['C1_geolat','C1G1SE_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1SE_MLT_comb','C1G1SE_MLTdiff_comb','C1_east_of_G1SE']









stop


    ;compare drift vs absolute conjunction times (SEEMS TO BE WORKING)
    tplot,['C1G1N_MLT_comb','C1G1N_drift_conjunction','C1G1N_absolute_conjunction','C1G1S_MLT_comb','C1G1S_drift_conjunction','C1G1S_absolute_conjunction']
    tplot,['C1G1NE_MLT_comb','C1G1NE_drift_conjunction','C1G1NE_absolute_conjunction','C1G1SE_MLT_comb','C1G1SE_drift_conjunction','C1G1SE_absolute_conjunction']




    tplot,['C1G1N_absolute_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1N_MLT_comb','C1G1N_MLTdiff_abs_comb']
    tplot,['C1G1S_absolute_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1S_MLT_comb','C1G1S_MLTdiff_abs_comb']
    tplot,['C1G1NE_absolute_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1NE_MLT_comb','C1G1NE_MLTdiff_abs_comb']
    tplot,['C1G1SE_absolute_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1SE_MLT_comb','C1G1SE_MLTdiff_abs_comb']
    tplot,['C1G1N_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1N_MLT_comb','C1G1N_MLTdiff_comb','C1_east_of_G1N']
    tplot,['C1G1S_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1S_MLT_comb','C1G1S_MLTdiff_comb','C1_east_of_G1S']
    tplot,['C1G1NE_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1NE_MLT_comb','C1G1NE_MLTdiff_comb','C1_east_of_G1NE']
    tplot,['C1G1SE_drift_conjunction','C1G1_Lcomb','C1G1_Ldiff_abs_comb','C1G1SE_MLT_comb','C1G1SE_MLTdiff_comb','C1_east_of_G1SE']


    ;Combine N/S conjunctions
    tplot,['C1G1N_absolute_conjunction','C1G1S_absolute_conjunction']
    tplot,['C1G1N_drift_conjunction','C1G1S_drift_conjunction']

    map_set,/mollweide,0,0,/grid,label=1,title='Title'
    map_continents
    t0z2 = time_double('2020-08-16/00:00')
    t1z2 = time_double('2020-08-17/00:00')

    geolatN = tsample('C1_out_iono_foot_north_geolat',[t0z2,t1z2],times=tms)
    geolonN = tsample('C1_out_iono_foot_north_geolon',[t0z2,t1z2],times=tms)
    geolatS = tsample('C1_out_iono_foot_south_geolat',[t0z2,t1z2],times=tms)
    geolonS = tsample('C1_out_iono_foot_south_geolon',[t0z2,t1z2],times=tms)
    plots,geolonS,geolatS,psym=2,color=250
    plots,geolonN,geolatN,psym=2,color=250,/continue

    geolatN = tsample('G1_out_iono_foot_north_geolat',[t0z2,t1z2],times=tms)
    geolonN = tsample('G1_out_iono_foot_north_geolon',[t0z2,t1z2],times=tms)
    geolatS = tsample('G1_out_iono_foot_south_geolat',[t0z2,t1z2],times=tms)
    geolonS = tsample('G1_out_iono_foot_south_geolon',[t0z2,t1z2],times=tms)
    plots,geolonS,geolatS,psym=2,color=50
    plots,geolonN,geolatN,psym=2,color=50,/continue











    ;---------------
    ;Find GEO separation at ionosphere
    get_data,'C1G1N_GSEdiff',data=tmp 
    store_data,'C1G1N_GSEdiff_magnitude',tmp.x,sqrt(tmp.y[*,0]^2 + tmp.y[*,1]^2 + tmp.y[*,2]^2)
    tplot,['C1_out_iono_foot_north_gse','G1_out_iono_foot_north_gse','C1G1N_GSEdiff_magnitude']

    get_data,'C1G1S_GSEdiff',data=tmp 
    store_data,'C1G1S_GSEdiff_magnitude',tmp.x,sqrt(tmp.y[*,0]^2 + tmp.y[*,1]^2 + tmp.y[*,2]^2)
    tplot,['C1_out_iono_foot_south_gse','G1_out_iono_foot_south_gse','C1G1S_GSEdiff_magnitude']






stop
end 
