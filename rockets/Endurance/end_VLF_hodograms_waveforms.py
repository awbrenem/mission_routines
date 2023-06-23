"""
Test the Endurance transfer function on DC and VLF data


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
#from scipy.fft import rfft, irfft
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle
import correlation_analysis





#---------------------------------------------
#Load VLF gain/phase corrected data
#Endurance_Analog 1_VLF12D_6-30000-100_gainphase_corrected.pkl 
#Endurance_Analog 1_VLF24D_6-30000-100_gainphase_corrected.pkl
#Endurance_Analog 1_VLF32D_6-30000-100_gainphase_corrected.pkl
#Endurance_Analog 1_VLF34D_6-30000-100_gainphase_corrected.pkl
#---------------------------------------------

pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'

fnsav = 'Endurance_Analog 1_VLF12D_6-30000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
wfc12 = wf_corr_load['wf']

tdat = wf_corr_load['tvals']
fs = np.mean([1/(tdat[i+1]-tdat[i]) for i in range(len(tdat)-1)])


fnsav = 'Endurance_Analog 1_VLF34D_6-30000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
wfc34 = wf_corr_load['wf']

fnsav = 'Endurance_Analog 1_VLF24D_6-30000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
wfc24 = wf_corr_load['wf']

fnsav = 'Endurance_Analog 1_VLF32D_6-30000-100_gainphase_corrected'
wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))
wfc32 = wf_corr_load['wf']




#Get spectral data for finding waves (use VLF12)
fspec, tspec, powerc = signal.spectrogram(wfc12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')


#----------------------------------------------------------------------
#Compare all VLF waveforms for desired timerange
#----------------------------------------------------------------------

def plot_wf(tr):

    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wfc12bp[goot]), list(wfc34bp[goot]), list(wfc24bp[goot]), list(wfc32bp[goot])])


    fig,axs = plt.subplots(6,figsize=(8,8))
    axs[0].plot(tdat[goot],wfc12bp[goot],color='b')
    axs[0].plot(tdat[goot],wfc34bp[goot],color='g')
    axs[0].title.set_text('VLF (12 vs 34)')
    axs[1].plot(tdat[goot],wfc12bp[goot],color='b')
    axs[1].plot(tdat[goot],wfc24bp[goot],color='r')
    axs[1].title.set_text('VLF (12 vs 24)')
    axs[2].plot(tdat[goot],wfc12bp[goot],color='b')
    axs[2].plot(tdat[goot],wfc32bp[goot],color='c')
    axs[2].title.set_text('VLF (12 vs 32)')
    axs[3].plot(tdat[goot],wfc34bp[goot],color='g')
    axs[3].plot(tdat[goot],wfc24bp[goot],color='r')
    axs[3].title.set_text('VLF (34 vs 24)')
    axs[4].plot(tdat[goot],wfc34bp[goot],color='g')
    axs[4].plot(tdat[goot],wfc32bp[goot],color='c')
    axs[4].title.set_text('VLF (34 vs 32)')
    axs[5].plot(tdat[goot],wfc24bp[goot],color='r')
    axs[5].plot(tdat[goot],wfc32bp[goot],color='c')
    axs[5].title.set_text('VLF (24 vs 32)')
    for i in range(6):
        axs[i].set_ylim(-maxv, maxv)


#----------------------------------------------------------------------
#Compare all hodograms for desired timerange
#----------------------------------------------------------------------

def plot_hod(tr):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wfc12bp[goot]), list(wfc34bp[goot]), list(wfc24bp[goot]), list(wfc32bp[goot])])

    fig,axs = plt.subplots(3,2,figsize=(8,8))
    axs[0,0].plot(wfc12bp[goot],wfc34bp[goot])    
    axs[0,1].plot(wfc12bp[goot],wfc24bp[goot])    
    axs[1,0].plot(wfc12bp[goot],wfc32bp[goot])    
    axs[1,1].plot(wfc34bp[goot],wfc24bp[goot])    
    axs[2,0].plot(wfc34bp[goot],wfc32bp[goot])    
    axs[2,1].plot(wfc24bp[goot],wfc32bp[goot])    
    axs[0,0].title.set_text('12 vs 34')
    axs[0,1].title.set_text('12 vs 24')
    axs[1,0].title.set_text('12 vs 32')
    axs[1,1].title.set_text('34 vs 24')
    axs[2,0].title.set_text('34 vs 32')
    axs[2,1].title.set_text('24 vs 32')
    for i in range(3):
        for j in range(2):
            axs[i,j].set(aspect='equal')
            axs[i,j].set_xlim(-maxv, maxv)
            axs[i,j].set_ylim(-maxv, maxv)



#Cross-spectral density
def plot_csd(tr, xlm=[0,12000],ylm=[-180,180],nps=512,wfs=['12','34']):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    window = np.hanning(len(goot[0]))
    #wf1z = wfc12[goot]*window
    wf1z = wfc24[goot]*window
    wf2z = wfc34[goot]*window
    csd, angle, f = correlation_analysis.cross_spectral_density(wf1z,wf2z,fs,nperseg=nps)

    goodv = np.where(csd > 0.6)


    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(f, csd)
    ax1.plot(f[goodv],csd[goodv],'.')
    ax1.set_xlim(xlm)
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(f, angle)
    ax2.plot(f[goodv],angle[goodv],'.')
    ax2.set_xlim(xlm)
    ax2.set_ylim(ylm)


#--------------------------------------
#Wave 1: Berstein (~120-130 sec at 5000-8000 Hz)
#--------------------------------------

ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='linear')


wfc12bp = filt.butter_bandpass_filter(wfc12, 4000, 6000, fs, order= 10)
wfc34bp = filt.butter_bandpass_filter(wfc34, 4000, 6000, fs, order= 10)
wfc24bp = filt.butter_bandpass_filter(wfc24, 4000, 6000, fs, order= 10)
wfc32bp = filt.butter_bandpass_filter(wfc32, 4000, 6000, fs, order= 10)

fspecz, tspecz, powercz = signal.spectrogram(wfc12bp, fs, nperseg=4*2048,noverlap=2*2048,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspecz,fspecz,np.abs(powercz),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='linear')



#Snapshot 1
tr = [122.62,122.72]  #zoomed out
tr = [122.67718,122.68]  #zoomed in
plot_wf(tr)
plot_csd([tr[0]-3, tr[1]+3])
plot_hod(tr) #Fairly circular



#Snapshot 2
tr = [123.5,123.54]
tr = [123.522,123.526] #Zoomed in
plot_wf(tr)
plot_csd([tr[0]-0.2, tr[1]+0.2])
plot_hod(tr)




plot_wf([123.00275,123.00375])
plot_wf([115.00275,128.00375])
plot_hod([123.00275,123.00375])






#NOTES: angle for Bernstein waves (wf12, wf34) is > 0 and increases towards 90. 
#Assuming wave vector propagates perp to Bo, this implies wavelengths < 5m * sin(90) = 5 m [Kintner+98]. 
#The DC power is at phase = 0, which means long wavelength (or large structures)


#Dynamic cross-spectral density
goot = np.where((tdat >= 100) & (tdat <= 200))
window = np.hanning(len(goot[0]))
wf1z = wfc12[goot]*window
wf2z = wfc34[goot]*window
timechunk = 1  #sec
Pxy, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf1z,wf2z,tdat[goot],fs,timechunk,nperseg=1024,plot=True)





"""
    from matplotlib import mlab

    nperseg = 1024

    # First create power spectral densities for normalization
    (ps1, f) = mlab.psd(wf1z, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    (ps2, f) = mlab.psd(wf2z, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    #plt.plot(f, ps1)
    #plt.plot(f, ps2)

    
    # Then calculate cross spectral density
    (csd, f) = mlab.csd(wf1z, wf2z,NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=False)
    csd_norm = np.absolute(csd)**2 / (ps1 * ps2)

    
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    # Normalize cross spectral absolute values by auto power spectral density
    ax1.plot(f, csd_norm)
    ax2 = fig.add_subplot(1, 2, 2)
    angle = np.angle(csd, deg=True)
    angle[angle<-90] += 360
    ax2.plot(f, angle)

"""

print("here")











#--------------------------------------
#Wave 3: 165.46 - 165.56 sec; ~100-1000 Hz
#--------------------------------------

wfc12bp = filt.butter_bandpass_filter(wfc12, 60, 300, fs, order= 10)
wfc34bp = filt.butter_bandpass_filter(wfc34, 60, 300, fs, order= 10)
wfc24bp = filt.butter_bandpass_filter(wfc24, 60, 300, fs, order= 10)
wfc32bp = filt.butter_bandpass_filter(wfc32, 60, 300, fs, order= 10)

plot_wf([165.20,165.8])
plot_hod([165.20,165.8])

