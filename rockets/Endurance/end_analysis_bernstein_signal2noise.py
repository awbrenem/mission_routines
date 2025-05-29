"""
Test the idea that V1 sees "extra" Bernstein power in the top band b/c it's slightly more sensitive than the other
probes and the signal is nearly at the noise level. 



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

#vDC = EFL('V12D')
#wfDC, tdatDC = vDC.load_data_gainphase_corrected()
#fsDC = vDC.chnspecs['fs']

v12 = EFL('VLF12D')
wfc, tdat = v12.load_data_gainphase_corrected()
wf, tdat = v12.load_data()
fs = v12.chnspecs['fs']

v1s = EFL('V1SD')
v2s = EFL('V2SD')
wf1sc, tdats = v1s.load_data_gainphase_corrected()
wf2sc, tdats = v2s.load_data_gainphase_corrected()

wf1s, tdats = v1s.load_data()
wf2s, tdats = v2s.load_data()



#Calibrated
v1 = EFL('V1SD')
v1s, t, = v1.load_data_gainphase_corrected()
v2 = EFL('V2SD')
v2s, t, = v2.load_data_gainphase_corrected()
v3 = EFL('V3SD')
v3s, t, = v3.load_data_gainphase_corrected()
v4 = EFL('V4SD')
v4s, t, = v4.load_data_gainphase_corrected()

#Raw
v1 = EFL('V1SD')
v1s, t, = v1.load_data()
v2 = EFL('V2SD')
v2s, t, = v2.load_data()
v3 = EFL('V3SD')
v3s, t, = v3.load_data()
v4 = EFL('V4SD')
v4s, t, = v4.load_data()

fs = v1.chnspecs['fs']





fspec, tspec, powerc = signal.spectrogram(wf, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-30,-5],yr=[4000,10000],xr=[100,200], yscale='log')




fspec, tspec, powerc1 = signal.spectrogram(v1s, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc2 = signal.spectrogram(v2s, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc3 = signal.spectrogram(v3s, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc4 = signal.spectrogram(v4s, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')


ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=[1e-7,1e-5],yr=[4000,10000],xr=[100,200], yscale='linear')



fig, axs = plt.subplots(4)



ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=[-13,5],yr=[4000,10000],xr=[100,200], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc2),vr=[-13,5],yr=[4000,10000],xr=[100,200], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc3),vr=[-13,5],yr=[4000,10000],xr=[100,200], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc4),vr=[-10,5],yr=[4000,10000],xr=[100,200], yscale='linear',ax=axs[3])



#-------------------------------------------------------------------------------------
#SLP fixed bias data (end of mission only)
#-------------------------------------------------------------------------------------

slpfb = end_data_loader.load_slp_fixedbias()


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


figs,axss = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=[700,900], yscale='linear',ax=axss[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[1,100],xr=[700,900], yscale='linear',ax=axss[1])
axss[2].plot(slpfb['ToF [s]'], slpfb['Normalized Density [\m3]'])
axss[2].set_xlim(700,900)

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





