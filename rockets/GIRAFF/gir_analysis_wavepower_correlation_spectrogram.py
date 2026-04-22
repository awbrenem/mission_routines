"""
Calculate spectra of the coherence and time lags of GIRAFF signals with respect to a reference signal 
NOTE: don't use DC coupled data b/c it's dominated by spin tone

NOTE: To see correlation with spin tone, use gir_analysis_spinphase_correlation_spectrogram.py
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import numpy as np 
import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec



pld = '381'


#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()



#determine nfft value 
nfft=2*4096
#nfft=512


fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=4)
power12 = np.abs(powerc)



#----------------------------------------------
#Define reference signal to correlate wave data against
#----------------------------------------------


freq_desired = 700 #Select frequency index to use as reference signal
goo = np.argmin(np.abs(fspec - freq_desired))

fref = f"{fspec[goo]:.0f}"
reference_signal = power12[goo,:]
#plt.plot(tspec, reference_signal)
#plt.xlim(200,260)
#plt.ylim(-0.01,0.01)

title='Correlation w/ 36.'+pld+' power at ' + fref + ' Hz\ngir_analysis_spinphase_correlation_spectrogram.py'


#------------------------------------------------------------
#Compute the coherence b/t different bands (e.g. power at all freqs vs power at 88 Hz)
#------------------------------------------------------------




bandw = 400  #freq chunk size of correlation window (Hz)
tchunk = 2 #time chunk size of correlation window (sec)
#NOTE: tchunk must be greater than tchunk_min = nfft * (tdat[1]-tdat[0])


f0 = 10 #min/max freqs to include in correlation analysis
f1 = 50000
corrmin = 0.4


lagmax = 1/bandw #max lag to consider in correlation (sec)


#calculate correlation/phase spectrograms
tvals, fvals, corrmax, tlagmax2 = ca.correlate_against_reference_signal_spectrogram(tspec, fspec, power12, reference_signal, bandw, f0, f1, tchunk, lagmax, corrmin)


#Shift tvals so that the times are centered on the middle of each chunk instead of the beginning
tvals = tvals - (tchunk/2)


yr = [50,50000]
vr = [-90,-50]
xr = [0,500]
fig,axs = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,power12,vr=vr,yscale='log',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",colorbar=1,ax=axs[0],title=title)
ps.plot_spectrogram(tvals,fvals,corrmax,vr=[0,1],yscale='log',zscale='linear',yr=yr,xr=xr,ylabel='Correlation\nfreq(Hz)',colorbar=1,ax=axs[1])
ps.plot_spectrogram(tvals,fvals,tlagmax2,vr=[-lagmax,lagmax],yscale='log',zscale='linear',yr=yr,xr=xr,ylabel='Lag (sec)\nfreq(Hz)',colorbar=1,ax=axs[2])
#add a horizontal line to the above power spectrum indicating the reference frequency
for i in range(3):
    axs[i].axhline(fspec[goo], color='magenta', linestyle='--', linewidth=0.9)




print('h')