"""
Test the Endurance transfer function by comparing two different channels 
(e.g. VLF data with DC data for identified waves)


#--------------------------------------
Possible waves to test (VLF12 and V12DC)
~630 Hz wave from 480-490 s

#--------------------------------------
#Test ~1000 Hz wave from 165.5 - 165.54. The V12DC and VLF12 gains match nicely, but the phases are different in the uncorrected data. 
#See if they line up here. 

#--------------------------------------
#TEST WAVES:
#Wave 1: 470-510 sec; ~4200 Hz 
#--------------------------------------

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import end_load_data as end
import end_load_gainphase as gainphase

from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
from end_fields_loader import Endurance_Fields_Loader as EFL
from math import remainder


#---------------------------------------------
#Load corrected data
#---------------------------------------------

c1 = 'V12D'
c2 = 'V34D'


v1c = EFL(c1)
wf1c, tdat1c = v1c.load_data_gainphase_corrected()
fs1 = v1c.chnspecs['fs']

v2c = EFL(c2)
wf2c, tdat2c = v2c.load_data_gainphase_corrected()
fs2 = v2c.chnspecs['fs']



#---------------------------------------------
#Load uncorrected data
#---------------------------------------------

v1 = EFL(c1)
v2 = EFL(c2)
wf1, tdat1 = v1.load_data()
wf2, tdat2 = v2.load_data()




#----------------------------------------------------------------------------------------
#Compare phase b/t two channels at frequency of interest
#----------------------------------------------------------------------------------------

def phase_comparison(vline_hz=0):


    v1c.plot_gainphase()
    v2c.plot_gainphase()


    #--Different channels have different x-axis (freq) values. Need to interpolate 
    #--to common base. To do this, need to unwrap phase angles to remove sharp discontinuities
    prad1u = np.unwrap(v1c.phase)
    prad2u = np.unwrap(v2c.phase)


    freq_interp = range(0,50000,1)
    interp1 = interp1d(v1c.freq_gainphase,prad1u,kind='cubic', bounds_error=False)
    prad1u = interp1(freq_interp)
    interp2 = interp1d(v2c.freq_gainphase,prad2u,kind='cubic', bounds_error=False)
    prad2u = interp2(freq_interp)

    prad_diff = prad1u - prad2u


    #--rewrap
    prad_diff = (prad_diff + np.pi) % (2 * np.pi) - np.pi



    #--Plot that clearly shows the phase differences b/t EDC and VLF channels (third panel)
    fig,axs = plt.subplots(2)
    axs[0].plot(v1c.freq_gainphase,v1c.phase, v2c.freq_gainphase,v2c.phase)
    axs[1].plot(freq_interp,prad_diff)
    axs[0].set_ylabel('phase\n(rad)\nBlue='+c1+' \nOrange='+c2)
    axs[1].set_ylabel('delta-phase (rad)')
    for i in range(2): axs[i].set_xlim(1,50000)
    for i in range(2): axs[i].set_xscale('log')

    for i in range(2):
        axs[i].axvline(vline_hz,linestyle='--', linewidth=0.8)
        axs[i].axhline(np.pi, linestyle='--', linewidth=0.8)
        axs[i].axhline(np.pi/2, linestyle='--', linewidth=0.8)
        axs[i].axhline(np.pi/4, linestyle='--', linewidth=0.4)
        axs[i].axhline(0, linestyle='--', linewidth=0.6)
        axs[i].axhline(-np.pi, linestyle='--', linewidth=0.8)
        axs[i].axhline(-np.pi/2, linestyle='--', linewidth=0.8)
        axs[i].axhline(-np.pi/4, linestyle='--', linewidth=0.4)


    print("here")

#--------------------------------------------------------------------
#Comparison plot of uncalibrated vs calibrated data for both channels
#--------------------------------------------------------------------

def wave_plot(fl, fh, x0=0, x1=800, y0=-10, y1=10, **kwargs):



    wf1cbp = filt.butter_bandpass_filter(wf1c, fl, fh, fs1, order= 10)
    wf1bp = filt.butter_bandpass_filter(wf1, fl, fh, fs1, order= 10)
    wf2cbp = filt.butter_bandpass_filter(wf2c, fl, fh, fs2, order= 10)
    wf2bp = filt.butter_bandpass_filter(wf2, fl, fh, fs2, order= 10)


    #--limit to times of interest for faster plotting
    goo1 = list(np.where((tdat1c >= x0) & (tdat1c <= x1)))[0]
    goo2 = list(np.where((tdat2c >= x0) & (tdat2c <= x1)))[0]



    fig,axs = plt.subplots(2,2,figsize=(8,8))
    axs[0,0].plot(tdat2c[goo2],wf2cbp[goo2], tdat2[goo2],wf2bp[goo2])
    axs[0,0].title.set_text(v2.chn + ' (corrected vs non-corrected)')
    axs[0,1].plot(tdat1c[goo1],wf1cbp[goo1], tdat1[goo1],wf1bp[goo1])
    axs[0,1].title.set_text(v1.chn + ' (corrected vs non-corrected)')
    axs[1,0].plot(tdat2[goo2],wf2bp[goo2], tdat1[goo1],wf1bp[goo1])
    axs[1,0].title.set_text(v2.chn + ' vs ' + v1.chn + ' (non-corrected)')
    axs[1,1].plot(tdat2c[goo2],wf2cbp[goo2], tdat1c[goo1],wf1cbp[goo1])
    axs[1,1].title.set_text(v2.chn + ' vs ' + v1.chn + ' (corrected)')
    for i in range(2):
        for j in range(2):
            axs[i,j].set_xlim(x0,x1)
            axs[i,j].set_ylim(y0,y1)


    plt.show()


#----------------------------------------------
#Non-bandpassed plots showing calibrated non-calibrated data
#Make sure there's no artificial sine waves present that can be introduced by 
#bandpassing or improper application of transfer function, etc. 
#----------------------------------------------

#fig,axs = plt.subplots(2)
#axs[0].plot(tdatDC,wfDC)
#axs[1].plot(tdatDC,wfDCc)
#for i in range(2): axs[i].set_xlim(500,500.4)
#axs[0].set_ylim(-5,5)
#axs[1].set_ylim(-5,5)



#-------------------------------
#Find wave events that can be used to test calibration
#-------------------------------

fspec1c, tspec1c, powerc1c = signal.spectrogram(wf1c, fs1, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
fspec1, tspec1, powerc1 = signal.spectrogram(wf1, fs1, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
fspec2c, tspec2c, powerc2c = signal.spectrogram(wf2c, fs2, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
fspec2, tspec2, powerc2 = signal.spectrogram(wf2, fs2, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')






#--------------------------------------
#TEST WAVES
#--------------------------------------


#--------------------------------------
#Wave 1: 470-510 sec; ~4200 Hz
#--------------------------------------

#wave 1
#ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=[-50,-35],yr=[0,14000],xr=[100,800], yscale='linear')
#fig,axs = plt.subplots(2)
#ps.plot_spectrogram(tspec1c,fspec1c,np.abs(powerc1c),vr=[-50,-25],yr=[3000,5000],xr=[400,600], yscale='linear',ax=axs[0])
#ps.plot_spectrogram(tspec2c,fspec2c,np.abs(powerc2c),vr=[-50,-25],yr=[3000,5000],xr=[400,600], yscale='linear',ax=axs[1])


phase_comparison(4200)
wave_plot(3000, 4500, x0=480, x1=480.02, y0=-0.05, y1=0.05)



#--------------------------------------
#Wave 2: 794-795 sec; ~100-300 Hz
#--------------------------------------

#ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=[-60,-20],yr=[50,200],xr=[794,795], yscale='linear')
#ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=[-60,-30],yr=[50,200],xr=[794,795], yscale='linear')

phase_comparison(300)
wave_plot(100, 300, x0=793.7, x1=793.9, y0=-0.02, y1=0.02)




#--------------------------------------
#Wave 3: 165.46 - 165.56 sec; ~80 Hz
#--------------------------------------

#fig, axs = plt.subplots(2)
#ps.plot_spectrogram(tspec1c,fspec1c,np.abs(powerc1c),vr=[-60,-20],yr=[50,2000],xr=[150,160], yscale='linear', ax=axs[0])
#ps.plot_spectrogram(tspec2c,fspec2c,np.abs(powerc2c),vr=[-60,-30],yr=[50,2000],xr=[150,160], yscale='linear', ax=axs[1])

phase_comparison(80)
wave_plot(60, 3000, x0=165.4, x1=165.6, y0=-0.8, y1=0.8)
wave_plot(60, 3000, x0=165.2, x1=165.8, y0=-0.8, y1=0.8)


print("end")










