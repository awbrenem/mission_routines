"""
Test the Endurance transfer function on DC and VLF data

"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import end_load_data as end
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
import plot_spectrogram as ps






fig,axs = plt.subplots(3)
axs[0].plot(tdat,wavedat)
axs[1].plot(tdat,wf)
axs[2].plot(tdat,wf_corr)

#for i in range(3): axs[i].set_ylim([-0.5,0.5])
#for i in range(3): axs[i].set_xlim([393.1+0.0625,393.1+0.065])


plt.plot(tdat,wf, tdat, wf_corr)
plt.xlim([393.1+0.0625,393.1+0.065])

##Phase shift test using normalized gain versions.
#goomax = np.max(wf)
#wfnorm = wf/goomax
#goomax = np.max(wf_corr)
#wf_corrnorm = wf_corr/goomax

#plt.plot(tdat,wfnorm, tdat,wf_corrnorm)
#plt.xlim([393.1+0.0625,393.1+0.065])
#plt.ylim(-0.02,0.02)



#Test VLF data now that transfer function has been applied by comparing to EDC12 channel. 
#There is a ~1000 Hz wave from 165.5 - 165.54. The EDC12 and VLF12 gains match nicely, but the phases are different in the uncorrected data. 
#See if they line up here. 

import filter_wave_frequency as filt

#----------------------------------------------------------------------------------------
#First let's compare the phase channels b/t EDC and VLF to see where they are most different 
#----------------------------------------------------------------------------------------

prad1, Hmag1, f1 = end_load_gainphase("Endurance_Analog 1_VLF12D_6-30000-100.txt")
prad2, Hmag2, f2 = end_load_gainphase("Endurance_Analog 1_V12D_10-10000-100.txt")

#unwrap phase angles
prad1u = np.unwrap(prad1)
prad2u = np.unwrap(prad2)

interp = interp1d(f2,prad2u,kind='cubic', bounds_error=False)
prad2ui = interp(f1)


#Plot that clearly shows the phase differences b/t EDC and VLF channels (third panel)
prad_diff = prad1u - prad2ui
fig,axs = plt.subplots(3)
axs[0].plot(f1,prad1, f2,prad2)
axs[1].plot(f1,prad1u, f2,prad2u, f1,prad2ui)
axs[2].plot(f1,prad_diff)
for i in range(3): axs[i].set_xlim(0,10000)



#----------------------------------------------------------------------------------------
#Load EDC data and plot EDC and VLF spectrograms to find common wave events to analyze
#----------------------------------------------------------------------------------------

vlf12c = wf_corr

#--------------------------------------
#TEST WAVES:
#Wave 1: 470-510 sec; ~4200 Hz (***NEED TO GAIN/PHASE CORRECT THE EDC DATA FOR THIS COMPARISON)
#--------------------------------------

edc = end.efield_dc()
edc12 = edc['dv12_mvm']
tedc = edc['times']
fs_edc = edc['samplerate']
plt.plot(tedc,edc12)

fsedc, tsedc, poweredc = signal.spectrogram(edc12, fs_edc, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tsedc,fsedc,np.abs(poweredc),vr=[-80,-20],yr=[10,6000],xr=[450,520], yscale='linear')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-80,-20],yr=[10,6000],xr=[450,520], yscale='linear')


edc12bp = filt.butter_bandpass_filter(edc12, 2000, 4500, fs_edc, order= 10)
vlf12cbp = filt.butter_bandpass_filter(vlf12c, 2000, 4500, fs, order= 10)

fig,axs = plt.subplots(2)
axs[0].plot(tedc,edc12bp)
axs[1].plot(tdat,vlf12cbp)
for i in range(2): axs[i].set_xlim(510,510.02)
for i in range(2): axs[i].set_ylim(-0.05,0.05)


"""
goo = list(map(tuple, np.where((tedc > 509) & (tedc < 511))))
edc12bpz = edc12bp[goo]
tedcz = tedc[goo]

goo = list(map(tuple, np.where((tdat > 509) & (tdat < 511))))
vlf12bpz = vlf12bp[goo]
tdatz = tdat[goo]

plt.plot(tedcz,edc12bpz, tdatz,vlf12bpz)
"""


#ft, tt, pt = signal.spectrogram(edc12bp, fs_edc, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
#ps.plot_spectrogram(tt,ft,np.abs(pt),vr=[-80,-20],yr=[10,6000],xr=[450,520], yscale='linear')
#ft, tt, pt = signal.spectrogram(vlf12bp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
#ps.plot_spectrogram(tt,ft,np.abs(pt),vr=[-80,-20],yr=[10,6000],xr=[450,520], yscale='linear')




highpass = filt.butter_bandpass_filter(waveform, 1, 10, fs, order= 10)



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


