"""
Endurance analysis of coupling b/t DC spikes and Bernstein power


"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np

from scipy.interpolate import interp1d
#from scipy.fft import rfft, irfft
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL





#---------------------------------------------
#Load VLF gain/phase corrected data
#---------------------------------------------

vDC = EFL('V12D')
wfDC, tdatDC = vDC.load_data_gainphase_corrected()
fsDC = vDC.chnspecs['fs']

v12 = EFL('VLF12D')
wf12, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']




#----------------------------------------------------------------------
#Get spectral data for finding waves (use VLF12)
#----------------------------------------------------------------------

fspec, tspec, powerc = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powercDC = signal.spectrogram(wfDC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')




#--------------------------------------
#Bernstein vertical spikes
#--------------------------------------
fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-30],yr=[4000,8000],xr=[110,130], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,0],yr=[0,100],xr=[110,130], yscale='linear',ax=axs[1])

wfDCbp = filt.butter_highpass_filter(wfDC,1,fsDC,order=10)
wf12bp = filt.butter_highpass_filter(wf12,4000,fs,order=10)

tr = [117, 120] 
tr = [117.575, 117.65] 
fig, axs = plt.subplots(2)
axs[0].plot(tdatDC,wfDCbp)
axs[1].plot(tdat,wf12bp)
for i in range(2): 
    axs[i].set_xlim(tr)

axs[0].set_ylim(-1,2.5)
axs[1].set_ylim(-0.4,0.4)





