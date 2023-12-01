
"""
Compare e- PES data (1 eV to 1 keV) to Bernstein waves

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader



#Get timeline of data to separate out science collection times 
tl, gsS, gsE = end_data_loader.load_timeline()

#Get PES e- data
#erange = [100,400]
#dfin0, dfout0 = end_data_loader.load_PES([106,230])
dfin0, dfout0 = end_data_loader.load_PES([400,420])


slp = end_data_loader.load_slp()

ig = end_data_loader.load_ig()

fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,400],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(ig.Alt,np.log10(ig.Dens_neutral))
axs[1].set_xlim(100,400)


fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,800],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(dfin0.sum(axis=1))



dfin1, dfout1 = end_data_loader.load_PES([10,200])
dfin2, dfout2 = end_data_loader.load_PES([200,400])
dfin3, dfout3 = end_data_loader.load_PES([400,600])
dfin4, dfout4 = end_data_loader.load_PES([600,800])
dfin5, dfout5 = end_data_loader.load_PES([800,1000])


#Energy bins
"""
0.1370,        8.4128,        9.1497,        9.9431,       11.5373,       13.3872,       15.5305,       
18.0196,       20.4224,       23.6898,       28.1379,       32.6274,       37.8539,       43.9113,       
50.9463,       59.1018,       68.5621,       79.5376,       92.2691,      107.0311,      123.9898,      
143.8542,      166.9112,      193.6813,      224.6471,      260.6623,      302.3216,      350.7751,     
406.9147,      472.0391,      547.5972,      635.2592,      736.9189,      854.8783,      991.6627,
"""


ephem = end_data_loader.load_ephemeris()

fig, axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-20],xr=[100,900],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(ephem['Flight Time'], ephem[' Altitde (km)'])
axs[1].set_xlim(100,900)
axs[2].plot(slp['ToF [s]'], slp['SLP Ni [/m3]'])
axs[3].plot(slp['ToF [s]'], slp['SLP Ti [K]'])
axs[4].plot(slp['ToF [s]'], slp['SLP Te [K]'])
axs[5].plot(slp['ToF [s]'], slp['SLP Vsc [V]'])

#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

v12DC = EFL('V12D')
wf12DC, tdatDC = v12DC.load_data_gainphase_corrected()

fsDC = v12DC.chnspecs['fs']

v12 = EFL('VLF12D')

wf12, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']


#----------------------------------------------------------------------
#Get spectral data for finding waves (use VLF12)
#----------------------------------------------------------------------

fspec, tspec, powerc12 = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powerc12DC = signal.spectrogram(wf12DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')






fig, axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,800],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(dfin1.sum(axis=1))
axs[2].plot(dfin2.sum(axis=1))
axs[3].plot(dfin3.sum(axis=1))
axs[4].plot(dfin4.sum(axis=1))
axs[5].plot(dfin5.sum(axis=1))

fig, axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,800],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(dfout1.sum(axis=1))
axs[2].plot(dfout2.sum(axis=1))
axs[3].plot(dfout3.sum(axis=1))
axs[4].plot(dfout4.sum(axis=1))
axs[5].plot(dfout5.sum(axis=1))

fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,800],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(dfout1.sum(axis=1))
axs[1].plot(dfout2.sum(axis=1))
axs[1].plot(dfout3.sum(axis=1))
axs[1].plot(dfout4.sum(axis=1))
axs[1].plot(dfout5.sum(axis=1))
axs[1].set_yscale('log')

fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-60,-30],xr=[100,800],yr=[1,12000],yscale='linear',ax=axs[0])
axs[1].plot(np.abs(dfout1.sum(axis=1)))
axs[1].plot(np.abs(dfout2.sum(axis=1)))
axs[1].plot(np.abs(dfout3.sum(axis=1)))
axs[1].plot(np.abs(dfout4.sum(axis=1)))
axs[1].plot(np.abs(dfout5.sum(axis=1)))
axs[1].set_yscale('log')




