"""
Test aliased signals seen on VDC12, 34, etc...




"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
#import end_load_data as end
from end_fields_loader import Endurance_Fields_Loader as EFL
#import end_load_gainphase as gainphase

from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle


#---------------------------------------------
#Load corrected data
#---------------------------------------------

"""
pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_DC/'
fnsav = 'Endurance_Analog 1_V12D_10-10000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
tdatDCc = wf_corr_load['tvals']
wfDCc12 = wf_corr_load['wf']
fsDC = np.mean([1/(tdatDCc[i+1]-tdatDCc[i]) for i in range(len(tdatDCc)-1)])

pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'
fnsav = 'Endurance_Analog 1_VLF12D_6-30000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
tdatc = wf_corr_load['tvals']
wfc = wf_corr_load['wf']
fs = np.mean([1/(tdatc[i+1]-tdatc[i]) for i in range(len(tdatc)-1)])
"""

vlf12 = EFL('VLF12D')
v12 = EFL('V12D')
v34 = EFL('V34D')

fs = vlf12.chnspecs['fs']


wfc, tvlf = vlf12.load_data_gainphase_corrected()
wfDCc12, tDC = v12.load_data_gainphase_corrected()
wfDCc34, tDC = v34.load_data_gainphase_corrected()

fsv = 30000
fsDC = 1/(tDC[1]-tDC[0])

#---------------------------------------------
#Load uncorrected data
#---------------------------------------------

"""
wavegoo = end.efield_vlf()
tdat = wavegoo.tvlf  #times
wf = wavegoo.dvlf12_mvm

wavegoo = end.efield_dc()
tdatDC = wavegoo['times']
wfDC12 = wavegoo['dv12_mvm']
wfDC34 = wavegoo['dv34_mvm']
"""

#-------------------------------
#FFT data
#-------------------------------

#fspecDC12, tspecDC12, powerDC12 = signal.spectrogram(wfDC12, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex',scaling='density')
#fspecDC34, tspecDC34, powerDC34 = signal.spectrogram(wfDC34, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex',scaling='density')
fspecDC12, tspecDC12, powercDC12 = signal.spectrogram(wfDCc12, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex',scaling='density')
fspec, tspec, powerc = signal.spectrogram(wfc, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex',scaling='density')
#fspecDC, tspecDC, powercDC = signal.spectrogram(wfDCc, fsDC, nperseg=512,return_onesided=True,mode='complex')
#fspec, tspec, powerc = signal.spectrogram(wfc, fs, nperseg=512,return_onesided=True,mode='complex')


#-------------------------------
#Context plot
#-------------------------------

t1 = 'VLF12 (dB=10log(V**2/Hz))'
t2 = 'VDC12 (dB=10log(V**2/Hz))'

fig, ax = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-45,-30],yr=[1000,8000],xr=[0,900], yscale='linear',title=t1,ax=ax[0])
ps.plot_spectrogram(tspecDC12,fspecDC12,np.abs(powercDC12),vr=[-45,-35],yr=[1000,8000],xr=[0,900], yscale='linear',title=t2,ax=ax[1])








#----------------------------------------------------------------------------
#Modified non-aliased VLF data so that it would appear as aliased in the VDC channel 
#----------------------------------------------------------------------------

#Take freqs of VLF data > 5000 Hz (limit of VDC chananels) and modify them

goo = list(map(tuple, np.where((fspec > 5000) & (fspec < 10000))))
ftst = np.squeeze(fspec[goo])
ptmp = np.squeeze(powerc[goo,:])
ftst2 = [5000 - (i - 5000) for i in ftst]


#Compare inverted VLF spectrum (>5000 Hz) to VDC spectrum
t1 = 'VLF12 (dB=10log(V**2/Hz))'
t2 = 'VLF12 inverted (dB=10log(V**2/Hz))'
t3 = 'VDC12 (dB=10log(V**2/Hz))'

fig, ax = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[1000,8000],xr=[100,200], yscale='linear',title=t1,ax=ax[0])
ps.plot_spectrogram(tspec,ftst2,np.abs(ptmp),vr=[-60,-20],yr=[1000,8000],xr=[100,200], yscale='linear',invert=True,title=t2,ax=ax[1])
ps.plot_spectrogram(tspecDC12,fspecDC12,np.abs(powerDC12),vr=[-60,-20],yr=[1000,8000],xr=[100,200], yscale='linear',title=t3,ax=ax[2])

fig, ax = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[1000,8000],xr=[680,800], yscale='linear',title=t1,ax=ax[0])
ps.plot_spectrogram(tspec,ftst2,np.abs(ptmp),vr=[-60,-20],yr=[1000,8000],xr=[680,800], yscale='linear',invert=True,title=t2,ax=ax[1])
ps.plot_spectrogram(tspecDC12,fspecDC12,np.abs(powerDC12),vr=[-45,-35],yr=[1000,8000],xr=[680,800], yscale='linear',title=t3,ax=ax[2])




#----------------------------------------------------------------------------------------
#First let's compare the phase channels b/t EDC and VLF to see where they are most different 
#----------------------------------------------------------------------------------------

prad1, Hmag1, f1 = gainphase.end_load_gainphase("Endurance_Analog 1_VLF12D_6-30000-100.txt")
prad2, Hmag2, f2 = gainphase.end_load_gainphase("Endurance_Analog 1_V12D_10-10000-100.txt")

Hmag1n = [i/17.5 for i in Hmag1]
Hmag2n = [i/1.75 for i in Hmag2]


plt.plot(f1,Hmag1n, f2,Hmag2n)
plt.xlim(0,12000)




#--------------------------------------
#TEST WAVES
#--------------------------------------

#--------------------------------------
#Wave 1: 140-160 sec; ~5000-7000 Hz (VLF12 channel)
#In this VDC12 channel this is seen b/t 3000-5000 Hz, t=120-180s
#According to the relative gain curves of the two channels, the amplitude of this bandpassed VDC12 
#waveform should be about 0.1 times that of the VLF12 channel
#---RESULT - 
#--------------------------------------


wfDCbp = filt.butter_bandpass_filter(wfDC12, 3000, 4900, fsDC, order= 10)
wfbp = filt.butter_bandpass_filter(wf, 5000, 7000, fs, order= 10)


#Test bandpassing
fspecgoo1, tspecgoo1, powergoo1 = signal.spectrogram(wfDCbp, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude',scaling='density')
fspecgoo2, tspecgoo2, powergoo2 = signal.spectrogram(wfbp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude',scaling='density')

t1 = 'VLF12 (dB=10log(V**2/Hz))'
t2 = 'VDC12 (dB=10log(V**2/Hz))'

fig, ax = plt.subplots(2)

ps.plot_spectrogram(tspecgoo2,fspecgoo2,powergoo2,vr=[-80,-20],yr=[1000,8000],xr=[100,200], yscale='linear', title=t1, ax=ax[0])
ps.plot_spectrogram(tspecgoo1,fspecgoo1,powergoo1,vr=[-80,-20],yr=[1000,8000],xr=[100,200], yscale='linear', title=t2, ax=ax[1])


fig,axs = plt.subplots(2)
axs[0].plot(tdat,wfbp)
axs[1].plot(tdatDC,wfDCbp)
axs[0].set_xlim(160,180)
axs[1].set_xlim(160,180)
axs[0].set_ylim(-0.15,0.15)
axs[1].set_ylim(-0.04,0.04)












"""
#HF channels 
Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt
Endurance_Analog 1_HF12_1000-20000000-100.txt
Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt
Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt


#VLF digitial files
Endurance_Analog 1_VLF12D_6-30000-100.txt
Endurance_Analog 1_VLF13D_6-30000-100.txt
Endurance_Analog 1_VLF41D_6-30000-100.txt
Endurance_Analog 1_VLF24D_6-30000-100.txt
Endurance_Analog 1_VLF42D_6-30000-100.txt
Endurance_Analog 1_VLF32D_6-30000-100.txt
Endurance_Analog 1_VLF34D_6-30000-100.txt


Endurance_Analog 1_VLF12A_6-100000-100.txt

Endurance_Analog 1_V42D_10-10000-100.txt
Endurance_Analog 1_V41D_10-10000-100.txt
Endurance_Analog 1_V34D_10-10000-100.txt
Endurance_Analog 1_V34A_10-10000-100.txt
Endurance_Analog 1_V32D_10-10000-100.txt
Endurance_Analog 1_V24D_10-10000-100.txt
Endurance_Analog 1_V13D_10-10000-100.txt
Endurance_Analog 1_V12D_10-10000-100.txt
Endurance_Analog 1_V12A_10-10000-100.txt
Endurance_Analog 1_V4SD_10-10000-100.txt
Endurance_Analog 1_V4SA_10-10000-100.txt
Endurance_Analog 1_V3SD_10-10000-100.txt
Endurance_Analog 1_V3SA_10-10000-100.txt
Endurance_Analog 1_V2SD_10-10000-100.txt
Endurance_Analog 1_V2SA_10-10000-100.txt
Endurance_Analog 1_V1SD_10-10000-100.txt
Endurance_Analog 1_V1SA_10-10000-100.txt

"""


