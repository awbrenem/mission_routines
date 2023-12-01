"""
end_analysis_bernstein_wave_properties.py 

Determine properties of the Bernstein waves based on the interferometry results

Short wavelength defined as k*pi >= 1 (e.g. Kintner+91)

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import numpy as np 
import plot_spectrogram as ps
import matplotlib.pyplot as plt

#-------------------------------------------------------
#Determine values for input into online plasma formulary 
#-------------------------------------------------------

t = 140 #sec 
wavelength = 10 
kmag = 2*np.pi/wavelength


ephem = end_data_loader.load_ephemeris(t=t)
alt = ephem[' Altitde (km)'] 


iri = end_data_loader.load_iri(alt=alt)
ig = end_data_loader.load_ig(alt=alt)
slp = end_data_loader.load_slp(t=140)



#------------------
#magnetic field 
#------------------
magv = EFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]
Boz = Bo[np.where(Bot > t)]
Boz = Boz[0]  #nT



#------------------
#temperatures (assume same temp for each ion species)
#------------------
#Te = iri['Te(K)'] / 11600
Te = slp['SLP Te [K]'] / 11600  #eV (SLP)

#However, at <350 km IRI suggests it is very close to neutral temp (but < Te)
Ti = iri['Ti(K)'] / 11600
Tneut_ig = ig['Temp_neutral'] / 11600
#Ion temps not measured by SLP at lower altitudes 
#Ti = slp['SLP Ti [K]']


#-------------------
#composition
#98% O+ at 270 km (comparison of flh to SLP densities)
#-------------------
fOp = 0.98
fHp = 0.02


#------------------
#Densities
#------------------

nneut = ig['Dens_neutral']
ni = slp['SLP Ni [/m3]'] / 100**3
ne = ni #cm-3

nOp = fOp * ne 
nHp = fHp * ne

#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

#gyroradii (m)
pe = 2.177e-2
pH = 0.647 
pO = 2.588  

#Debye length (m) 
debye_e = 6.882e-3
debye_H = 4.773e-3
debye_O = 4.773e-3

#Thermal velocities
Ve = 193 #km/s
VH = 3.13   
VO = 0.78

inertial_e = 10.6  #m
inertial_H = 458 
inertial_O = 1831

kpe = pe * kmag
kpH = pH * kmag
kpO = pO * kmag

#waves considered "short wavelength" relative to O+ gyroradius
