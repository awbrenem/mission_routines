#Identify waves to test the gain/phase transfer functions.

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
#from scipy import signal
import numpy as np 
#import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
#import plasma_params_get_flhr_freq as dflh
#import pyIGRF
#import pickle
#from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
from scipy.interpolate import interp1d





#Load E-fields data
pld = '381'
c1 = 'VLF12D'
c2 = 'V12D'
v1 = GFL(pld,c1)
wf1, tdat1 = v1.load_data()
v2 = GFL(pld,c2)
wf2, tdat2 = v2.load_data()


v1.plot_gainphase()
v2.plot_gainphase()




#----------------------------------------------
#Spectra of channels to find waves seen on both
#----------------------------------------------

fspec1, tspec1, powerc1, fs = fftspec.fft_spectrum_piecewise(tdat1, wf1, fs_thres=0.1, nfft=16384, noverlap=8)
fspec2, tspec2, powerc2, fs = fftspec.fft_spectrum_piecewise(tdat2, wf2, fs_thres=0.3, nfft=16384, noverlap=8)



#----------------------------------
#General plot for identifying waves
#----------------------------------

xr = [100,550]
yr = [30,50000]
#yr = [800,1400]
ys = 'log'
vr1 = [-90,-50]
vr2 = [-90,-50]
#minzval = -80
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)



#------------------------
#Wave 1; t=172, f=1100 Hz
#------------------------

xr = [170,175]
#yr = [30,50000]
yr = [800,1400]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
#minzval = -80
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)




#------------------------
#Wave 2; t=270, f=15000 Hz
#------------------------

xr = [265,275]
#yr = [30,50000]
yr = [10000,20000]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
#minzval = -80
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)



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



#------------------------
#Wave 4; t=463.5, f=100 Hz (broadband DC)
#------------------------

xr = [463,465]
yr = [30,500]
#yr = [800,1400]
ys = 'linear'
vr1 = [-90,-50]
vr2 = [-90,-50]
#minzval = -80
fig, axs = plt.subplots(2, figsize=(9,7))
ps.plot_spectrogram(tspec1,fspec1,np.abs(powerc1),vr=vr1,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c1+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[0])
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec2,fspec2,np.abs(powerc2),vr=vr2,yscale=ys,yr=yr,xr=xr,ylabel='power spectrum '+pld+' '+c2+'\nfreq(Hz)\ndB of (mV/m)^2/Hz',ax=axs[1])
fig.tight_layout(pad=0)


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

