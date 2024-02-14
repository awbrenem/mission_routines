"""
Analyze polarization properties of Bernstein waves using the spectra. 
Specifically, the quantity (E12 - E34) / (E12 + E34) is very useful for seeing small differences


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
import interferometry_routines as interf





#--------------------------------------------------------------
#Get timeline of data to separate out science collection times 
#--------------------------------------------------------------

tl, gsS, gsE, bsS, bsE = end_data_loader.load_timeline()
ephem1, ephem2 = end_data_loader.load_ephemeris()



#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

Bo = EFL('mag')
bx,by,bz,tdatb = Bo.load_data()
bmag = np.sqrt(bx**2 + by**2 + bz**2)
fcH = 28*bmag / 1836


v12 = EFL('VLF12D')
v13 = EFL('VLF13D')
v41 = EFL('VLF41D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')

fs = v12.chnspecs['fs']
wf12, tdat = v12.load_data_gainphase_corrected()
wf13, tdat = v13.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf41, tgoo = v41.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
wf42 = -wf24
wf14 = -wf41



#----------------------------------------------------------------------
#Get spectral data for finding waves
#----------------------------------------------------------------------

fspec, tspec, powerc12 = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc13 = signal.spectrogram(wf13, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc34 = signal.spectrogram(wf34, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc42 = signal.spectrogram(wf42, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc14 = signal.spectrogram(wf14, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc32 = signal.spectrogram(wf32, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')



#----------------------------------------------------------------------------------------
#Calculate the fractional difference b/t different spectra to identify artificial waves
#----------------------------------------------------------------------------------------

ptmp_diff = np.abs(powerc12) - np.abs(powerc34)
ptmp_sum = np.abs(powerc12) + np.abs(powerc34)
ptmp_fracdiff = ptmp_diff/ptmp_sum

#version with low values removed
ptmp_fracdiff2 = ptmp_fracdiff.copy()
valmin = 0.05
ptmp_fracdiff2[np.abs(ptmp_fracdiff2) < valmin] = float('nan')



#yr = [10,15000]
#xr = [100,200]
yr = [0,10000]
xr = [100,900]

ys = 'linear'
cmap = 'RdYlGn'
title = '(E12-E34)/(E12+E34)\nGreen is more power on E12\nRed is more power on E34'

ps.plot_spectrogram(tspec,fspec,ptmp_fracdiff2,
                    vr=[-1,1],yr=yr,xr=xr, yscale=ys,xlabel='time(s)',
                    ylabel='f(Hz)',zscale='linear',title=title + '\nvals > |' + str(valmin) + '| only',
                    plot_kwargs={'cmap':'RdYlGn'})



#Compare the rising <1 kHz signal on (mostly) E34 to the falling Berstein power

xr = [100,220]
fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,ptmp_fracdiff2,vr=[-0.5,0.5],ax=axs[0],
                    yr=[2000,-1000],xr=xr, yscale=ys,xlabel='time(s)',
                    ylabel='f(Hz)',zscale='linear',plot_kwargs={'cmap':'bwr'})
ps.plot_spectrogram(tspec, fspec, np.abs(powerc12),vr=[-40,-30],ax=axs[1],
                    yr=[5000,8000],xr=xr, yscale=ys,xlabel='time(s)',
                    ylabel='f(Hz)',zscale='log')
axs[0].plot(tdatb,fcH)


#Compare to attitude changes

xr = [100,900]
fig, axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,ptmp_fracdiff2,vr=[-1,1],ax=axs[0],
                    yr=yr,xr=xr, yscale=ys,xlabel='time(s)',
                    ylabel='f(Hz)',zscale='linear',plot_kwargs={'cmap':cmap},title=title + '\nvals > |' + str(valmin) + '| only')
axs[1].plot(ephem2['Time'],ephem2['Yaw'])
axs[2].plot(ephem2['Time'],ephem2['Roll'])
axs[3].plot(ephem2['Time'],ephem2['Pitch'])
axs[4].plot(ephem2['Time'],ephem2['RollRate'])
axs[5].plot(ephem2['Time'],ephem2['AoA_T'])

for i in range(len(axs)):
    axs[i].set_xlim(xr)

axs[1].set_ylabel('yaw')
axs[2].set_ylabel('roll')
axs[3].set_ylabel('pitch')
axs[4].set_ylabel('roll rate')
axs[5].set_ylabel('AoA_T')
