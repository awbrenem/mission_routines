"""
end_analysis_wake_shadowing.py 

Determine (vs time) if a wake or shadow can exist, and what type
See Matts Andre + 22 (submitted manuscript) for wake descriptions:

Narrow wake:  mvi^2/2 > kTi; mvi^2/2 >> e|Vsc|
Enhanced wake: eVsc >> mvi^2/2 > kTi   (would only happen in low density plasmas)
Focussing wake: -eVsc > mvi^2/2 > kTi

ExB drift velocity:
vi = 800 - 1000 m/s (upleg/downleg during Bernstein waves). Mostly zonal on upleg and half zonal half meridional on downleg

Rocket velocity:
Vrocket ~ 3500 m/s 

Temperature:
Upleg: 
Downleg (840 s):
    Ti ~ 0.1 eV
    Te ~ 0.2 eV
 
Spacecraft potential:
Vsc = -1 V


"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import numpy as np 
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import plasma_params_get_flhr_freq as dflh


e1eV = 1.6e-19            #joules in 1 eV
me     = 9.1093897e-31    #kg
mp     = 1.6726231e-27
kB = 1.380658e-23         #J/K


#----------------------------------------------
#load data
#----------------------------------------------

skins = EFL('V1SD')
v1sd, times_skin = skins.load_data()
v2sd, times_skin = skins.load_data()
v3sd, times_skin = skins.load_data()
v4sd, times_skin = skins.load_data()


ephem, ephem2 = end_data_loader.load_ephemeris()
iri = end_data_loader.load_iri()
ig = end_data_loader.load_ig()
slp = end_data_loader.load_slp()

#Define time arrays
times_slp = slp['ToF [s]']
times_iri_up = iri['times_upleg']
times_iri_dn = iri['times_downleg']
times_ephem = ephem['Time']



#----------------------------------------------
#Determine rocket skin potential energy
#----------------------------------------------


skinpot = (v1sd + v2sd + v3sd + v4sd)/4

plt.plot(times_skin,skinpot)
plt.ylim(-1,0)
plt.xlim(100,900)

Epot = np.abs(skinpot)


#-------------------------------------------
#Determine flow energy based on relative velocity b/t rocket and convection
#...for now I'm just using the rocket velocity
#-------------------------------------------

vmag_rocket = np.sqrt(ephem['ECEF-X-VEL']**2 + ephem['ECEF-Y-VEL']**2 + ephem['ECEF-Z-VEL']**2)
plt.plot(times_ephem,vmag_rocket)

EHflow = 0.5 * (mp) * vmag_rocket**2 / e1eV
EOflow = 0.5 * (mp*16) * vmag_rocket**2 / e1eV





#----------------------------------------------------------
#Thermal energies (assume same temp for each ion species)
#----------------------------------------------------------

Te_iri = iri['Te(K)'] / 11600
Te_slp = slp['SLP Te [K]'] / 11600  #eV (SLP)


#However, at <350 km IRI suggests it is very close to neutral temp (but < Te)
#Tneut_ig = ig['Temp_neutral'] / 11600
Ti_slp = slp['SLP Ti [K]'] / 11600
Ti_iri = iri['Ti(K)'] / 11600

#Temp comparison in eV
plt.plot(times_iri_up,Ti_iri,'.')
plt.plot(times_iri_dn,Ti_iri,'.')
plt.plot(times_slp,Ti_slp,'.')
plt.plot(times_slp,Te_slp,'.')
plt.plot(times_slp,Te_iri,'.')


#Ion temps not measured by SLP at lower altitudes, but EISCAT suggests 1100 K from 100-600 km on both upleg/downleg
#Ti = 1100 / 11600 #EISCAT
#Ti[times_slp < 400] = 1100/11600


#Thermal velocities from temps (km/s)
VH_iri = np.sqrt(kB*(11600*Ti_iri)/mp) / 1000
VO_iri = np.sqrt(kB*(11600*Ti_iri)/(16*mp)) / 1000

Eth = 0.5 * (mp) * (VH_iri*1000)**2  / e1eV
Eth = 0.5 * (mp*16) * (VO_iri*1000)**2 / e1eV

plt.plot(times_iri_up,VH_iri,'.')
plt.plot(times_iri_dn,VH_iri,'.')
plt.plot(times_iri_up,VO_iri,'.')
plt.plot(times_iri_dn,VO_iri,'.')




#-------------------------------------------------
#Plot energy comparison 
#-------------------------------------------------


#Comparison plot of all relevant energies
plt.plot(times_skin,Epot,'.')       #rocket charging (blue)
plt.plot(times_ephem,EOflow,'.') #rocket energy (orange)
plt.plot(times_iri_up,Eth,'.') #thermal energy (red/green)
plt.plot(times_iri_dn,Eth,'.')
plt.ylim(0,1)
plt.xlim(100,900)




