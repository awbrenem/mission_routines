#Endurance - find density from upper and lower hybrid lines 

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
import end_load_data
import plasma_params_get_density_from_flhr_freq as dflh
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
from scipy import signal

import numpy as np 
import matplotlib.pyplot as plt

import plasmapy
from astropy import units as u  

#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]

"""Enable auto module reloading"""
#%load_ext autoreload
#%autoreload 2





"""Load and plot E-field VLF data"""
evlf = end_load_data.efield_vlf()
fs = evlf.samplerate
freq12, tspec12, power12 = signal.spectrogram(evlf.dvlf12_mvm, fs, nperseg=16384, noverlap=16384/2., return_onesided=True)
ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-60], xr=[100,900],yr=[0,15000],pl=1)




"""Load and plot E-field HF data"""
ehf = end_load_data.efield_hf()

hf12 = ehf.afftpow12
hf34 = ehf.afftpow34
hffreqs = ehf.afreq
hftimes = ehf.atimesfft

ps.plot_spectrogram(hftimes,hffreqs,hf12,pl=0,vr=[-80,-30],yr=[2e3,2.6e6],xr=[100,900])



#First try to identify density from fuh
tv = [400, ]
fuh = [1.706, 1.630]

fpe = [plasmapy.formulary.plasma_frequency(ne_lh, particle='electron', to_hz=True)]






#Find density at t~500

#print(mag.keys())

tref = 200  #reference time (sec since launch)
goo = np.squeeze(np.where(mag.tsec > tref))
print(type(goo))
print(np.shape(goo))

print(mag.tsec[goo[0]])
print(mag.Bmag[goo[0]])


#flh = [8100.] * u.Hz   #@500 km
flh = [7600.] * u.Hz   #@800 km
Bo = [mag.Bmag[goo[0]]] * u.nT
print("Bo = ", Bo)

nH_ne = [0.09] * u.dimensionless_unscaled 
nO_ne = [0.91] * u.dimensionless_unscaled

fce = plasmapy.formulary.gyrofrequency(Bo, particle='electron', to_hz=True)
fcH = plasmapy.formulary.gyrofrequency(Bo, particle='H+', to_hz=True)
fcHe = plasmapy.formulary.gyrofrequency(Bo, particle='He+', to_hz=True)
fcO = plasmapy.formulary.gyrofrequency(Bo, particle='O+', to_hz=True)

print(fce, fcH, fcHe, fcO)

ne_lh = dflh.dens_IonMassFractions(flh, fce, nH_ne, nO_ne)
ne_lh2 = dflh.dens_singleion(flh, Bo, 'O+')
print("ne_lh O+", ne_lh2)


#fuh = [1.5e6, 2.2e6]
fuh = 2e6 * u.Hz
fpe = np.sqrt(fuh**2 - fce[0]**2)
ne_uh = (fpe/8980)**2
ne_uh = ne_uh.value * u.cm**-3


print(type(ne_lh[0]))
print("ne_lh fractional ", ne_lh)

print("ne_uh ", ne_uh)

ne_lh2 = ne_lh[0].value 
ne_lh = ne_lh2 * u.cm**-3

fpe = [plasmapy.formulary.plasma_frequency(ne_lh, particle='electron', to_hz=True)]
#fpe = 8980*np.sqrt(ne_lh)

#fuh = np.sqrt(0.4e6**2 + fce[0]**2)
#print(fpe/1e6, fce[0]/1e6, fuh/1e6)

