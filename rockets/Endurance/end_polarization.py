"""
Determine polarization of various waves from Endurance.
"""



import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
import end_load_data
#import plasma_params_get_density_from_flhr_freq as dflh
import numpy as np 
import matplotlib.pyplot as plt
import plasmapy
from astropy import units as u  
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
from scipy import signal



#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]

"""Enable auto module reloading"""
#%load_ext autoreload
#%autoreload 2



mag = end_load_data.mag_dc()
#Index(['tsec', 'Bx', 'By', 'Bz', 'Bmag'], dtype='object')
evlf = end_load_data.efield_vlf()
#Index(['tsec', 'amp'], dtype='object')
ehf = end_load_data.efield_hf()
#dict_keys(['afftpow12', 'afftpow34', 'afreq', 'atimesfft', 'fftsize', 'overlap', 'weight', 'samplerate', 'nfreq', 'nfftlines', 'in_file'])


"""
fig, ax = plt.subplots(2, sharey=True, sharex=True)
ax[0].plot(evlf.tvlf, evlf.dvlf12_mvm)
ax[1].plot(evlf.tvlf, evlf.dvlf34_mvm)
ax[0].set_xlim(320, 340)
ax[0].set_ylim(-1,1)
"""


#Spectral overview
fs = evlf.samplerate
freq12, tspec12, power12 = signal.spectrogram(evlf.dvlf12_mvm, fs, nperseg=512) #, return_onesided=1)
#freq34, tspec34, power34 = signal.spectrogram(evlf.dvlf34_mvm, fs, nperseg=512) #, return_onesided=1)
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40], xr=[900,1000],yr=[6000,10000], yscale='linear')
ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[3500,14000],xr=[600,1000], yscale='linear')

#vals = plt.ginput(2)

#ps.plot_spectrogram(tspec34,freq34,power34,vr=[-100,-40], xr=[0,900],yr=[0,10000], yscale='linear')


#Isolate waveform for polarization

#tr = [320, 340]
tr = [100, 800]
fr = [4000, 8000]

pmask = np.full_like(power12, 1)

badt = np.squeeze(np.where(tspec12 < tr[0]))
if len(badt) != 0: pmask[:,badt] = 0
badt = np.squeeze(np.where(tspec12 > tr[1]))
if len(badt) != 0: pmask[:,badt] = 0
badf = np.squeeze(np.where((freq12 < fr[0])))
if len(badt) != 0: pmask[badf,:] = 0
badf = np.squeeze(np.where((freq12 > fr[1])))
if len(badt) != 0: pmask[badf,:] = 0

p12 = power12 * pmask

ps.plot_spectrogram(tspec12,freq12,p12,vr=[-100,-40], xr=tr,yr=fr, yscale='linear')


pmask = np.full_like(power34, 1)

badt = np.squeeze(np.where(tspec34 < tr[0]))
pmask[:,badt] = 0
badt = np.squeeze(np.where(tspec34 > tr[1]))
pmask[:,badt] = 0
badf = np.squeeze(np.where((freq34 < fr[0])))
pmask[badf,:] = 0
badf = np.squeeze(np.where((freq34 > fr[1])))
pmask[badf,:] = 0

p34 = power34 * pmask

ps.plot_spectrogram(tspec34,freq34,p34,vr=[-100,-40], xr=tr,yr=fr, yscale='linear')



print("here")




