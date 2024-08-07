"""
end_analysis_bernstein_wave_properties.py 

Determine properties of the Bernstein waves based on the interferometry results
(including those needed for WHAMP input)

Bernstein waves typically "short wavelength", defined as k*pi ~ n, where n is an (e.g. Kintner+91)

WHAMP values for 160 sec 
Te = 0.22 eV (SLP - also consistent with IRI)
Tneut = 0.088 (IG)
Ti = 0.1 (EISCAT - consistent with IRI)
ne = 220058 cm-3
H+frac = 0.04
O+frac = 0.96
Bo = 49375 nT
alt = 345 km

WHAMP values for 840 sec. 
Te = 0.2 eV
Ti = 0.09 eV
ne = 246920 cm-3
H+frac = 0.02
O+frac = 0.98
Bo = 47669 nT
alt = 264

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
#import plasma_params_get_flhr_freq
import plasma_params_get_flhr_freq as dflh

#-------------------------------------------------------
#Determine values for input into online plasma formulary 
#-------------------------------------------------------

#t = 160 #sec
t = 750 #sec
wavelength = 5.7  #m
fo = 6200  #Hz

#t = 140 #sec 
#wavelength = 10 
#kmag = 2*np.pi/wavelength
#fo = 6300 #Hz (freq of interest)

ephem, ephem2 = end_data_loader.load_ephemeris(t=t)
alt = ephem['Altitude'] 
#alt = 264

iri = end_data_loader.load_iri(alt=alt)
ig = end_data_loader.load_ig(alt=alt)
slp = end_data_loader.load_slp(t=t)



#------------------
#magnetic field 
#------------------
magv = EFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]
Boz = Bo[np.where(Bot > t)]
Boz = Boz[0]  #nT
#Boz = 47669 nT

plt.plot(Bot,mag[1],Bot,mag[2])
plt.plot(Bot,Bo)
#------------------
#temperatures (assume same temp for each ion species)
#------------------

#Te = iri['Te(K)'] / 11600
Te = slp['SLP Te [K]'] / 11600  #eV (SLP)


#However, at <350 km IRI suggests it is very close to neutral temp (but < Te)
#Ti = iri['Ti(K)'] / 11600
Tneut_ig = ig['Temp_neutral'] / 11600
#Tneut_ig = 0.088 eV
#Ion temps not measured by SLP at lower altitudes, but EISCAT suggests 1100 K
Ti = 1100 / 11600 #EISCAT
#Ti = 0.09 eV

#Ti = slp['SLP Ti [K]'] / 11600


#-------------------
#composition
#98% O+ at 270 km (comparison of flh to SLP densities)
#-------------------
fOp = 0.99
fHp = 0.01


#------------------
#Densities
#------------------

nneut = ig['Dens_neutral']
nneut = 79279090 
ni = slp['SLP Ni [/m3]'] / 100**3
#ni = 246920 cm-3
ne = ni #cm-3

nOp = fOp * ne 
nHp = fHp * ne


#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

#Gyrofrequencies (Hz)
fce = 28*Boz
fcH = fce / 1634
fcHe = fcH / 4
fcO = fcH / 16
#fce = 1.4e6 
#fcH = 769
#fcO = 48 

#plasma freqs (Hz)
fpe = 4.47e6
fpH = 104219
fpO = 26054

#beta 
beta_e = fpe / fce


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
VO = 0.78  #consistent with EISCAT

#inertial lengths
inertial_e = 10.6  #m
inertial_H = 458 
inertial_O = 1831

#k * p comparisons
#waves considered "short wavelength" relative to O+ gyroradius
#Note that Bernstein waves have k*ps ~ 1, 2, 3, ....
kmag = 1

kpe = pe * kmag
kpH = pH * kmag  #=0.4
kpO = pO * kmag  #=1.6

#phase velocities
vphase = fo*(2*np.pi) / kmag / 1000   #km/s
vphase = 63 #km/s


#Lower hybrid freq (f >> fH)
#flh = dflh.flhr_IonMassFractions([ne,ne], [fce,fce], [fHp,fHp], [fOp,fOp])
flh = dflh.flhr_IonMassFractions(ne, fce, fHp, fOp)
#flh = 8338 @ 160 sec


#Ion-ion crossover freq (H+ - O+)

t1 = (2*np.pi*fcH) * (2*np.pi*fcO) 
t2 = (nHp + 16*nOp) / (nOp + 16*nHp)
fHO = np.sqrt(t1 * t2)/ (2*np.pi)
#fHO = 713 Hz @ 160 sec


#Cross-over freq from RH R-mode to LH polarization 
t1 = (nHp/ne) * (2*np.pi*fcO)**2 
t2 = (nOp/ne) * (2*np.pi*fcH)**2
fcr = np.sqrt(t1 + t2) 
#fcr = 4807 Hz 

#Cutoff freq
t1 = (nHp/ne) * (2*np.pi*fcO)
t2 = (nOp/ne) * (2*np.pi*fcH)
fL2 = t1 + t2 
#fL2 = 4786 Hz








