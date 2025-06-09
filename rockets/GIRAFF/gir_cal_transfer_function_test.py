"""
Test the GIRAFF transfer function by comparing two different channels 
(e.g. VLF data with DC data for identified waves)


#--------------------------------------
Possible waves to test (VLF12 and V12DC)
...
#--------------------------------------

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from math import remainder
import fft_spectrum_piecewise as fftspec


#---------------------------------------------
#Load corrected data
#---------------------------------------------

#Load E-fields data
pld = '381'
c1 = 'VLF34D'
c2 = 'V34D'
v1 = GFL(pld,c1)
v2 = GFL(pld,c2)


#Load uncorrected data
wf1, tdat1 = v1.load_data()
wf2, tdat2 = v2.load_data()
#Load gain/phase corrected data
wf1c, tdat1c = v1.load_data_gainphase_corrected()
wf2c, tdat2c = v2.load_data_gainphase_corrected()





#----------------------------------------------------------------------------------------
#Compare phase b/t two channels at frequency of interest
#----------------------------------------------------------------------------------------

def phase_comparison(vline_hz=0):

    v1.plot_gainphase()
    v2.plot_gainphase()


    #--Different channels have different x-axis (freq) values. Need to interpolate 
    #--to common base. To do this, need to unwrap phase angles to remove sharp discontinuities
    prad1u = np.unwrap(v1.phase)
    prad2u = np.unwrap(v2.phase)


    freq_interp = range(0,50000,1)
    interp1 = interp1d(v1.freq_gainphase,prad1u,kind='cubic', bounds_error=False)
    prad1u = interp1(freq_interp)
    interp2 = interp1d(v2.freq_gainphase,prad2u,kind='cubic', bounds_error=False)
    prad2u = interp2(freq_interp)

    prad_diff = prad1u - prad2u


    #--rewrap
    prad_diff = (prad_diff + np.pi) % (2 * np.pi) - np.pi



    #--Plot that clearly shows the phase differences b/t EDC and VLF channels (third panel)
    fig,axs = plt.subplots(2)
    axs[0].plot(v1.freq_gainphase,v1.phase, v2.freq_gainphase,v2.phase)
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

def wave_plot(fl, fh, x0=0, x1=500, y0=-10, y1=10, **kwargs):

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
    axs[1,0].title.set_text('non-corrected: ' + v2.chn + ' vs ' + v1.chn)
    axs[1,1].plot(tdat2c[goo2],wf2cbp[goo2], tdat1c[goo1],wf1cbp[goo1])
    axs[1,1].title.set_text('corrected: ' + v2.chn + ' vs ' + v1.chn)
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



fspec1, tspec1, powerc1, fs1 = fftspec.fft_spectrum_piecewise(tdat1, wf1, fs_thres=0.1, nfft=16384, noverlap=8)
fspec1c, tspec1c, powerc1c, fs1c = fftspec.fft_spectrum_piecewise(tdat1c, wf1c, fs_thres=0.1, nfft=16384, noverlap=8)
fspec2, tspec2, powerc2, fs2 = fftspec.fft_spectrum_piecewise(tdat2, wf2, fs_thres=0.3, nfft=16384, noverlap=8)
fspec2c, tspec2c, powerc2c, fs2c = fftspec.fft_spectrum_piecewise(tdat2c, wf2c, fs_thres=0.3, nfft=16384, noverlap=8)

fs1 = fs1[0]
fs1c = fs1c[0]
fs2 = fs2[0]
fs2c = fs2c[0]



#--------------------------------------
#TEST WAVES (get these from gir_cal_identify_waves_for_transfer_function.py)
#--------------------------------------


#--------------------------------------
#Wave 1; t=172, f=1100 Hz
#Should see very little phase difference b/t VLF12D and V12D
#--------------------------------------


xr = [170,175]
#yr = [30,50000]
yr = [800,1400]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)



phase_comparison(1100)
#wave_plot(3000, 4500, x0=480, x1=480.02, y0=-0.05, y1=0.05)
wave_plot(800, 1400, x0=170, x1=170.03, y0=-0.02, y1=0.02)



#--------------------------------------
#Wave 2: 270 sec; ~15000 Hz
#High freq test - NOT WORKING WELL
#--------------------------------------

xr = [265,275]
yr = [10000,20000]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)

phase_comparison(15000)
wave_plot(14800, 15000, x0=270.001, x1=270.003, y0=-0.06, y1=0.06)




#------------------------
#Wave 3; t=282, f=100 Hz
#------------------------

xr = [280,285]
yr = [10,400]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)


phase_comparison(100)
wave_plot(100, 120, x0=282.0, x1=282.3, y0=-0.15, y1=0.15)
wave_plot(100, 120, x0=283.0, x1=283.3, y0=-0.15, y1=0.15)


print("end")


#------------------------
#Wave 4; t=463.5, f=100 Hz (broadband DC)
#------------------------

xr = [463,465]
yr = [30,500]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)


phase_comparison(60)
wave_plot(60, 70, x0=463.6, x1=464, y0=-0.1, y1=0.1)
phase_comparison(100)
wave_plot(100, 410, x0=463.65, x1=463.75, y0=-1.4, y1=1.4)
wave_plot(100, 410, x0=463.72, x1=463.78, y0=-0.2, y1=0.2)
phase_comparison(1000)
wave_plot(1000, 1010, x0=463.7, x1=463.71, y0=-0.1, y1=0.1)
#Following has no phase modification of the DC channel, but significant mod of the VLF channel
wave_plot(13500, 13510, x0=463.7, x1=463.703, y0=-0.01, y1=0.01)


#------------------------
#Wave 5; t=138, f=100 Hz (broadband DC)
#------------------------

xr = [135,140]
yr = [30,500]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)


#------------------------
#Wave 6; t=107, f=100 Hz (boom deploy)
#------------------------

xr = [106,108]
yr = [30,5000]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)

phase_comparison(400)
wave_plot(400, 402, x0=106.5, x1=106.58, y0=-0.1, y1=0.1)
