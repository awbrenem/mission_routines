"""
end_analysis_crossover_freq.py 

Determine cross-over freq below which waves are LH polarized and above are RH polarized

Also determine if it is plausible that a remote source (along Bo) could have a crossover frequency at > 5 kHz
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
import pickle
import pyIGRF


iri = end_data_loader.load_iri(alt=800)
tiri_u = iri['times_upleg']
tiri_d = iri['times_downleg']
iriH = iri['H_ions'] #percentage
iriO = iri['O_ions']
ne = iri['Ne(m-3)']

#nH_ne = (ne*iriH/100) / ne



plt.plot(tiri_u,iriH)



ephem = end_data_loader.load_ephemeris()
alt = ephem[0]['Altitude'] 


plot_alt = np.asarray(ephem[0]['Altitude'])

glat = 78 + (55/60)
glon = 11 + (55/60)
BoIGRF = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])

ephemt = ephem[0]['Time']


fceIGRF = 28 * BoIGRF
fcHIGRF = fceIGRF / 1836
fcOIGRF = fcHIGRF / 16





"""
path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/plasma_composition.pkl'
with open(path, 'rb') as file:
    compdat = pickle.load(file)


fHp = compdat['nH_ne']
fOp = compdat['nO_ne']
ne = compdat['ne_langmuirprobe']
compdat_times = compdat['times']

nOp = fOp * ne 
nHp = fHp * ne


#extrapolate cyclotron values for times of composition determination

fcH_interp = np.zeros(len(fHp))
fcO_interp = np.zeros(len(fHp))

for i in range(len(fHp)):
    fcH_interp[i] = np.interp(compdat_times[i],ephemt,fcHIGRF)
    fcO_interp[i] = np.interp(compdat_times[i],ephemt,fcOIGRF)

#Cross-over freq from RH R-mode to LH polarization 
t1 = (nHp/ne) * (2*np.pi*fcO_interp)**2 
t2 = (nOp/ne) * (2*np.pi*fcH_interp)**2
fcr = np.sqrt(t1 + t2) 


plt.plot(compdat_times,fcr)

"""

#--------------------------------
#Cross-over and cutoff freqs based on IRI
#--------------------------------

t = np.arange(100,150,10)
#t = np.arange(100,900,5)
fcr = np.zeros(len(t))
fco = np.zeros(len(t))

for i in range(len(t)):

    tref = np.where(ephem[0]['Time'] > t[i])[0][0]
    altz = alt[tref]

    iri = end_data_loader.load_iri(alt=altz)

    fcH = np.interp(t[i],ephemt,fcHIGRF)
    fcO = np.interp(t[i],ephemt,fcOIGRF)

    #---cross over freq (below is LH above is RH)
    t1 = iri['H_ions']/100. * (2*np.pi*fcO)**2
    t2 = iri['O_ions']/100. * (2*np.pi*fcH)**2
    fcr[i] = np.sqrt(t1 + t2) / (2*np.pi)
    #---cutoff freq (below crossover freq and the lower limit for LH propagation)
    x1 = iri['H_ions']/100. * (2*np.pi*fcO)
    x2 = iri['O_ions']/100. * (2*np.pi*fcH)
    fco[i] = (x1 + x2) / (2*np.pi)


plt.plot(t,fcr,t,fco)

plt.plot(ephemt,fcHIGRF)
plt.plot(ephemt,fcOIGRF)
plt.plot(t,fcr)
plt.plot(t,fco)

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/crossover_cutoff_freq_iri.pkl'
pickle.dump([t,fcr,fco],open(path,'wb'))






#------------------------------------------


iri_all = end_data_loader.load_iri()



































tiriu = tiri_u.to_numpy()
tirid = tiri_d.to_numpy()



#compdat_timesiriu = np.concatenate([tiriu[0:134],np.flip(tirid)])
compdat_timesiriu = tiriu[0:134]
compdat_timesirid = np.flip(tirid)


fcH_interpiriu = np.zeros(len(iriH))
fcO_interpiriu = np.zeros(len(iriH))

for i in range(len(fHp)):
    fcH_interpiriu[i] = np.interp(compdat_timesiriu[i],ephemt,fcHIGRF)
    fcO_interpiriu[i] = np.interp(compdat_timesiriu[i],ephemt,fcOIGRF)



t1iriu = iriH[0:134]/100. * (2*np.pi*fcO_interpiriu[0:134])**2
t2iriu = iriO[0:134]/100. * (2*np.pi*fcH_interpiriu[0:134])**2
fcririu = np.sqrt(t1iriu + t2iriu) 



plt.plot(compdat_timesiriu,fcririu)



"""
ephem, ephem2 = end_data_loader.load_ephemeris(t=t)
alt = ephem['Altitude'] 
#alt = 264

iri = end_data_loader.load_iri()
#ig = end_data_loader.load_ig(alt=alt)
slp = end_data_loader.load_slp()
"""






















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



#Cutoff freq
t1 = (nHp/ne) * (2*np.pi*fcO)
t2 = (nOp/ne) * (2*np.pi*fcH)
fL2 = t1 + t2 
#fL2 = 4786 Hz









