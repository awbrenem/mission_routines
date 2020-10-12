
;Useful SPEDAS tutorial on loading GOES data
;http://spedas.org/wiki/index.php?title=GOES#Electron.2C_Proton.2C_Alpha_Detector_.28EPEAD.29


goes_load_data, trange=['2014-01-11', '2014-01-12'], datatype='fgm', probes='12', /avg_1m
tplot, 'g12_H_enp'



goes_load_data, trange=['2014-01-11', '2014-01-12'], datatype='fgm', probes='12', /avg_1m
tplot, 'g12_H_enp'


;Energetic Particle Sensor (EPS)
;Load and plot the 1-min averaged 2.4MeV proton flux as measured by the EPS instrument onboard GOES-12 on March 22, 2008:
goes_load_data, trange=['2014-01-11', '2014-01-12'], datatype='eps', probes='12', /avg_1m
tplot, 'g12_prot_2.4MeV_flux'


;We can also plot the integral flux of electrons at 0.6 MeV:
tplot, 'g12_elec_0.6MeV_iflux'


;X-ray Sensor (XRS)
;See the XRS Readme for more information on the GOES XRS instrument.
;To load and plot the 1-min averaged X-ray flux as measured by the XRS instrument onboard GOES-10 on March, 22, 2008:
goes_load_data, trange=['2014-01-11', '2014-01-12'], datatype='xrs', probes='10', /avg_1m
tplot, 'g10_xrs_avg'


;GOES 13-15 SEM Data
;To check the current status of the GOES 13-15 data, see GOES 13-15 status
;For a full description of the GOES 13-15 instruments and their data products, see the GOES 13-15 databook


;Fluxgate Magnetometer (FGM)
;Load and plot 1-min averaged data for both FGM sensors on GOES-15 for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='fgm', probes='15', /avg_1m
tplot, 'g15_H_enp_*'


;Magnetospheric Electron Detector (MAGED)
;Load and plot 1-min averaged, corrected (for deadtimes) flux of 40keV, 75keV electrons measured by the MAGED instrument on GOES-15 for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='maged', probes='15', /avg_1m
tplot, ['g15_maged_40keV_dtc_cor_flux', 'g15_maged_75keV_dtc_cor_flux']


;To calculate the pitch angles corresponding to each MAGED telescope head (FGM data must be loaded):
goes_lib ; compile the GOES library routines
goes_pitch_angles, 'g15_H_enp_1', 'g15_HT_1', prefix = 'g15'
tplot, 'g15_pitch_angles'



;MAGED data at 40 keV, 75 keV, March 17, 2013
;GOES-15 pitch angles, March 17, 2013


;Magnetospheric Proton Detector (MAGPD)
;Load and plot the 95 keV proton flux from the MAGPD instrument onboard GOES-15 for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='magpd', probes='15', /avg_1m
tplot, 'g15_magpd_95keV_dtc_cor_flux'


;Electron, Proton, Alpha Detector (EPEAD)
;Warning about "east" and "west" look directions
;The GOES 13-15 satellites are capable of undergoing a yaw flip, resulting in the EPEAD pointing directions labeled as 'east' and 'west' being switched. Care must be taken to ensure 'east' and 'west' refer to the same directions when comparing EPEAD data from one time interval to another.
;To read more on this issue, see the note at NGDC: Note on GOES 13‐15 Solar Protons and Yaw Flips



;Load and plot the 2.5MeV, uncorrected proton flux as observed by the GOES-15 EPEAD instrument for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='epead', probes='15', /avg_1m
tplot, 'g15_prot_2.5MeV_uncor_flux'

;To calculate and plot the center pitch angles for the east and west heads of the EPEAD instrument, first load the FGM data:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='fgm', probes='15', /avg_1m

;then compile the GOES library file, goes_lib:
goes_lib ; compile GOES support routines
;and finally use goes_epead_center_pitch_angles to calculate the pitch angles from the magnetic field in spacecraft coordinates:
goes_epead_center_pitch_angles, 'g15_Bsc_1', 'g15_BTSC_1'
tplot, 'goes_epead_center_pitch_angles'



;High Energy Proton and Alpha Detector (HEPAD)
;Load and plot the 375MeV proton flux as observed by the GOES-15 HEPAD instrument for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='hepad', probes='15', /avg_1m
tplot, 'g15_hepadp_375MeV_flux'


;X-ray Sensor (XRS)
;Load and plot the X-rays measured by the GOES-15 XRS instrument for the day of March 17, 2013:
goes_load_data, trange=['2013-03-17', '2013-03-18'], datatype='xrs', probes='15', /avg_1m
tplot, 'g15_xrs_avg'
