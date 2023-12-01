"""
See if there's any association b/t Bernstein power and near-DC power 


Todo:

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
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader




#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

vDC = EFL('V12D')
wfDC, tdatDC = vDC.load_data_gainphase_corrected()
fsDC = vDC.chnspecs['fs']

v12 = EFL('VLF12D')
wf, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']

v1s = EFL('V1SD')
v2s = EFL('V2SD')
wf1s, tdats = v1s.load_data_gainphase_corrected()
wf2s, tdats = v2s.load_data_gainphase_corrected()


#----------------------------------------------------------------------
#Mission timeline data
#----------------------------------------------------------------------

tl, gsS, gsE = end_data_loader.load_timeline()

"""
125 195
205 275
285 355
365 435
445 515
525 595
625 695
705 775
785 855
865 900.6
"""

fspec, tspec, powerc = signal.spectrogram(wf, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,10000],xr=[100,850], yscale='log')


#-------------------------------------------------------------------------------------
#Density structures
#--can be somewhat difficult to identify in skin data
#-------------------------------------------------------------------------------------




#wDC = filt.butter_bandpass_filter(wfDC, 1, 80, fsDC, order= 8)
wDC = filt.butter_bandpass_filter(wfDC, 10,80, fsDC, order= 8)
w = filt.butter_bandpass_filter(wf, 5000, 8000, fs, order= 8)
w1s = filt.butter_bandpass_filter(wf1s, 5, 80, fs, order= 8)
w2s = filt.butter_bandpass_filter(wf2s, 5, 80, fs, order= 8)


fspec, tspec, powerc = signal.spectrogram(wf, fs, nperseg=256,noverlap=128,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powercDC = signal.spectrogram(wfDC, fsDC, nperseg=256,noverlap=128,window='hann',return_onesided=True,mode='complex')


figs,axss = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=[700,900], yscale='linear',ax=axss[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[1,100],xr=[700,900], yscale='linear',ax=axss[1])



#tr = [705,775]
tr = [110,115]
fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=tr, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[0,100],xr=tr, yscale='linear',ax=axs[1])
axs[3].plot(tdat,wf)
axs[2].plot(tdatDC,wfDC)
axs[3].set_ylim(-1.3,1.3)
axs[2].set_ylim(-15.5,15.5)
for i in range(2): axs[i+2].set_xlim(tr)


tr = [113,115]
fig,axs = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=tr, yscale='linear',ax=axs[0])
axs[1].plot(tdat,wf)
axs[2].plot(tdat,w)
axs[1].set_ylim(-1.3,1.3)
axs[1].set_xlim(tr)
axs[2].set_ylim(-0.1,0.1)
axs[2].set_xlim(tr)




plt.plot(tdat,w, tdatDC,wDC)
plt.xlim(705,775)
plt.ylim(-0.25,0.25)





