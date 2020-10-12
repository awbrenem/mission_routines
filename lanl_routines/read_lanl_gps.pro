;Read LANL GPS data files. 
;These files have 4 minute low resolution MeV particle and magnetic field data. 
;Not super useful most of the time, so this routine is not entirely completed. 


;Morley, S. K., J. P. Sullivan, M. R. Carver, R. M. Kippen, R. H. W.Friedel, G. D. Reeves,
;and M. G. Henderson (2017), Energetic Particle Data From the Global Positioning System Constellation, 
;Space Weather, 15, 283–289, doi:10.1002/2017SW001604.


;path = 'http://www.ngdc.noaa.gov/stp/ space-weather/satellite-data/ satellite-systems/gps/'

;GPS satellites:
;ns41,ns48,ns53,ns54,ns55,ns56,ns57,ns58,ns59,ns60,ns61,ns62,ns63,ns64,ns65,ns66,ns67,ns68,ns69,ns70,ns71,ns72,ns73,




;---------------------------------------------------------
;Download online files 
;****************************
;*****
;***
;---------------------------------------------------------

;remote_data_dir = 'https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/gps/data/'
payload = 'ns63'
;file_loaded = spd_download(remote_path=remote_data_dir+payload+'/',remote_file='*.ascii',/no_download)

;local_data_dir = '/Users/aaronbreneman/data/lanl/'+payload + '/'
;file_loaded = spd_download(remote_path=remote_data_dir+payload+'/',remote_file=file,$
;    local_path=local_data_dir,/last_version)





;fn = '~/Desktop/lanl_gps_sample.txt'
path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events/Jan6/'
;fn = 'ns66_140105_v1.03.ascii.txt'
;fn = 'ns65_140105_v1.03.ascii.txt'
fn = payload+'_140105_v1.03.ascii.txt'

openr,lun,path+fn,/get_lun 

;Skip past header files
jnk = ''
for i=0,231 do readf,lun,jnk


datatmp = ''
doy = '' & year = '' & lshell = '' & mlt = ''
geolat = '' & geolon = '' & rad_re = '' & maglon = ''
bmag = ''
proton_density_fit = '' & electron_density_fit = ''
electron_temperature_fit = '' & proton_characteristic_momentum = ''
dropped_data = ''
proton_integrated_flux_fit = strarr(1,6)
integral_flux_instrument = strarr(1,30) & integral_flux_energy = strarr(1,30)
electron_diff_flux_energy = strarr(1,15) & electron_diff_flux = strarr(1,15)



i=0.
while not eof(lun) do begin 
    readf,lun,datatmp
    data = strsplit(datatmp,' ',/extract)

    doy = [doy,data[0]]
    year = [year,data[22]]
    geolat = [geolat,data[1]]
    geolon = [geolon,data[2]]
    maglon = [maglon,data[28]]
    rad_re = [rad_re,data[3]]
    lshell = [lshell,data[29]]
    mlt = [mlt,data[35]]
    bmag = [bmag,data[37]]  ;Gauss
    proton_density_fit = [proton_density_fit,data[57]]  ;cm^-3
    electron_density_fit = [electron_density_fit,data[59]]  ;cm^-3
    electron_temperature_fit = [electron_temperature_fit,data[58]] ;MeV
    proton_characteristic_momentum = [proton_characteristic_momentum,data[56]] ;MeV/c (proton characteristic momentum)
    dropped_data = [dropped_data,data[25]] ;1= bad data
    proton_integrated_flux_fit = [proton_integrated_flux_fit,reform(data[92:97],1,6)] ;"proton rates above 10,15.85,25.11,30,40,79.43 MeV", "UNITS": "cm^-2sec^-1sr^-1"
    integral_flux_instrument = [integral_flux_instrument,reform(data[135:164],1,30)] ;integral of electron flux fit above integral_flux_energy[i]: cm^-2sec^-1sr^-1
    integral_flux_energy = [integral_flux_energy,reform(data[165:194],1,30)] ;energies for integral_flux_instrument; MeV
    electron_diff_flux_energy = [electron_diff_flux_energy,reform(data[195:209],1,15)] ;energies for the fluxes in electron_diff_flux_energy: MeV
    electron_diff_flux = [electron_diff_flux,reform(data[210:224],1,15)] ;electron flux at energies electron_diff_flux[i]: cm^-2sec^-1sr^-1MeV^-1

    i++ 
endwhile
close,lun & free_lun,lun



;Remove leading blank string 
nelem = n_elements(doy)
doy = doy[1:nelem-1]
year = year[1:nelem-1]
geolat = float(geolat[1:nelem-1])
geolon = float(geolon[1:nelem-1])
maglon = float(maglon[1:nelem-1])
rad_re = float(rad_re[1:nelem-1])
lshell = float(lshell[1:nelem-1])
mlt = float(mlt[1:nelem-1])
bmag = float(bmag[1:nelem-1])
proton_density_fit = float(proton_density_fit[1:nelem-1])
electron_density_fit = float(electron_density_fit[1:nelem-1])
electron_temperature_fit = float(electron_temperature_fit[1:nelem-1])
proton_characteristic_momentum = float(proton_characteristic_momentum[1:nelem-1])
dropped_data = float(dropped_data[1:nelem-1])
proton_integrated_flux_fit = float(proton_integrated_flux_fit[1:nelem-1,*])
integral_flux_instrument = float(integral_flux_instrument[1:nelem-1,*])
integral_flux_energy = float(integral_flux_energy[1:nelem-1,*])
electron_diff_flux_energy = float(electron_diff_flux_energy[1:nelem-1,*])
electron_diff_flux = float(electron_diff_flux[1:nelem-1,*])

;Create time variable 
doy = double(doy)

;add one or two leading zeros for correct doy format
goo = where((doy lt 100.) and (doy ge 10.)) 
doy = strtrim(doy,2)
doy[goo] = '0' + doy[goo]
goo = where(doy lt 10.) 
doy = strtrim(doy,2)
doy[goo] = '00' + doy[goo]


dtst = year + strtrim(doy,2)
dtst = double(dtst)

datetimestr = strarr(nelem-1)
for i=0,nelem-2 do datetimestr[i] = date_conv(dtst[i],'F')

datetime = time_double(datetimestr)


;Extract various energy channels for labeling
ife = reform(integral_flux_energy[0,*])  ;MeV




;Create tplot variables 
store_data,payload+'_geolat',datetime,geolat & ylim,payload+'_geolat',-90,90
store_data,payload+'_geolon',datetime,geolon & ylim,payload+'_geolon',0,360
store_data,payload+'_maglon',datetime,maglon & ylim,payload+'_maglon',0,360
store_data,payload+'_radius_re',datetime,rad_re & ylim,payload+'_radius_re',0,10
store_data,payload+'_lshell',datetime,lshell & ylim,payload+'_lshell',0,20
store_data,payload+'_mlt',datetime,mlt & ylim,payload+'_mlt',0,24
store_data,payload+'_bmag',datetime,bmag*1e5 & ylim,payload+'_bmag',0,1000 ;nT
store_data,payload+'_proton_density',datetime,proton_density_fit & ylim,payload+'_proton_density',0,10
store_data,payload+'_electron_density',datetime,electron_density_fit & ylim,payload+'_electron_density',0,1
store_data,payload+'_electron_temp',datetime,electron_temperature_fit & ylim,payload+'_electron_temp',0,5
store_data,payload+'_proton_characteristic_momentum',datetime,proton_characteristic_momentum
store_data,payload+'_dropped_data',datetime,dropped_data
store_data,payload+'_proton_integrated_flux_fit',datetime,proton_integrated_flux_fit
store_data,payload+'_integral_flux_instrument',datetime,integral_flux_instrument & ylim,payload+'_integral_flux_instrument',10,1d8,1
store_data,payload+'_electron_diff_flux_energy',datetime,electron_diff_flux_energy & ylim,payload+'_electron_diff_flux_energy',10,1d8,1
store_data,payload+'_electron_diff_flux',datetime,electron_diff_flux & ylim,payload+'_electron_diff_flux',10,1d8,1

;options,'integral_flux_instrument','ztitle',string(ife)
;options,'integral_flux_instrument','labels',['lab1','lab2']  ; Labels used for individual plot quantities in a single

rbsp_efw_init
stop
end 



;#     "rate_electron_measured": { "DESCRIPTION":"rate_electron_measured",
;#        "DIMENSION": [11],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 4
;#     },
;#     "rate_proton_measured": { "DESCRIPTION":"rate_proton_measured",
;#        "DIMENSION": [5],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 15
;#     },
;#     "LEP_thresh": { "DESCRIPTION":"LEP_threshold: 0 means low, 1 means high",
;#        "DIMENSION": [1],
;#        "UNITS": "none",
;#        "START_COLUMN": 20
;#     },;

;#     "b_coord_radius": { "DESCRIPTION":"radius from earths dipole axis",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 26
;#     },
;#     "b_coord_height": { "DESCRIPTION":"height above earth's equatorial plane",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 27
;#     },
;#     "b_equator": { "DESCRIPTION":"B field at equator",
;#        "DIMENSION": [1],
;#        "UNITS": "gauss",
;#        "START_COLUMN": 38
;#     },








;#     "proton_activity": { "DESCRIPTION":"proton_activity=1 if proton activity is occurring",
;#        "DIMENSION": [1],
;#        "UNITS": "none",
;#        "START_COLUMN": 55
;#     },
;#     "model_counts_electron_fit_pf": { "DESCRIPTION":"E1-E11 rates due to proton background based on proton flux fit -- not filled",
;#        "DIMENSION": [11],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 60
;#     },
;#     "model_counts_proton_fit_pf": { "DESCRIPTION":"proton rate due from proton flux fit",
;#        "DIMENSION": [5],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 71
;#     },
;#     "model_counts_electron_fit": { "DESCRIPTION":"electron rates from 9-parameter electron flux model",
;#        "DIMENSION": [11],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 76
;#     },

;#     "efitpars": { "DESCRIPTION":"fit parameters for 9 parameter electron fit",
;#        "DIMENSION": [9],
;#        "UNITS": "various",
;#        "START_COLUMN": 225
;#     },
;#     "pfitpars": { "DESCRIPTION":"fit parameters for 4 parameter proton fit",
;#        "DIMENSION": [4],
;#        "UNITS": "various",
;#        "START_COLUMN": 234
;#     }
;#};
;
;
;
;





;#     "electron_background": { "DESCRIPTION":"electron_background estimated",
;#        "DIMENSION": [11],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 39
;#     },
;#     "proton_background": { "DESCRIPTION":"proton_background estimated",
;#        "DIMENSION": [5],
;#        "UNITS": "hertz",
;#        "START_COLUMN": 50
;#     },
;#     "L_LGM_TS04IGRF": { "DESCRIPTION":"L_LGM_TS04IGRF, TS04 External field, IGRF Internal field",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 30
;#     },
;#     "L_LGM_OP77IGRF": { "DESCRIPTION":"L_LGM_OP77IGRF, OP77 External field, IGRF Internal field",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 31
;#     },
;#     "L_LGM_T89CDIP": { "DESCRIPTION":"L_LGM_T89CDIP T89 External Field, Centered Dipole Internal Field",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 32
;#     },
;#     "L_LGM_T89IGRF": { "DESCRIPTION":"L_LGM_T89IGRF T89 External Field, IGRF Internal Field",
;#        "DIMENSION": [1],
;#        "UNITS": "R_E",
;#        "START_COLUMN": 33
;#     },
;#     "bfield_ratio": { "DESCRIPTION":"Bsatellite/Bequator",
;#        "DIMENSION": [1],
;#        "UNITS": "none",
;#        "START_COLUMN": 34
;#     },
;#     "utc_lgm": { "DESCRIPTION":"LGM UTC",
;#        "DIMENSION": [1],
;#        "UNITS": "h",
;#        "START_COLUMN": 36
;#     },;


