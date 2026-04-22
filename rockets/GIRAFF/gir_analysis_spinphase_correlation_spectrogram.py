"""
Calculate spectra of the coherence and time lags of GIRAFF signals with the rocket spin period
or with other reference signals ("Option 2" and reference_signal2; e.g. 88 Hz power)

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
import pickle
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
import filter_wave_frequency as filt

pld = '381'
spinperiod = 1.5



#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()



#determine nfft value 
nfft=2*4096
#print('number of chunks/spin = ', spinperiod / tchunk)

fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=4)
power12 = np.abs(powerc)



#----------------------------------------------
#Define reference signal to correlate wave data against
#Use half-spin tone directly from DC magnetic field data
#----------------------------------------------


#Load magnetic field and interpolate to times of spectra
magv = GFL(pld,'mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]

#Highpass filter Bo so that I can see the spin period more easily on plots
reference_signal = filt.butter_highpass_filter(Bo, 0.05, 1/(Bot[1]-Bot[0]))
reference_signal = np.interp(tspec, Bot, reference_signal)
#Half spin tone
reference_signal = np.abs(reference_signal)



#------------------------------------------------------------
#Compute the coherence b/t different bands (e.g. power at all freqs vs power at 88 Hz)
#------------------------------------------------------------


bandw = 100  #freq chunk size of correlation window (Hz)
tchunk = 6   #time chunk size of correlation window (sec)
#NOTE: tchunk must be greater than tchunk_min = nfft * (tdat[1]-tdat[0])


#min/max freqs to include in correlation analysis
f0 = 10
f1 = 50000
corrmin = 0.4

lagmax = spinperiod   #max lag to consider in correlation (sec)

tvals, fvals, corrmax, tlagmax2 = ca.correlate_against_reference_signal_spectrogram(tspec, fspec, power12, reference_signal, bandw, f0, f1, tchunk, lagmax, corrmin)


#Shift tvals so that the times are centered on the middle of each chunk instead of the beginning
tvals = tvals - (tchunk/2)

title='Correlation w/ 36.'+pld+' half spin tone\ngir_analysis_spinphase_correlation_spectrogram.py'


yr = [50,50000]
vr = [-90,-50]
xr = [0,500]
fig,axs = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,power12,vr=vr,yscale='log',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",colorbar=1,ax=axs[0],title=title)
ps.plot_spectrogram(tvals,fvals,corrmax,vr=[0,1],yscale='log',zscale='linear',yr=yr,xr=xr,ylabel='Correlation\nfreq(Hz)',colorbar=1,ax=axs[1])
ps.plot_spectrogram(tvals,fvals,tlagmax2,vr=[-lagmax,lagmax],yscale='log',zscale='linear',yr=yr,xr=xr,ylabel='Lag (sec)\nfreq(Hz)',colorbar=1,ax=axs[2])



print('h')