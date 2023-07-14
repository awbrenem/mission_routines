"""
Interferometry analysis on Endurance


Todo:
--Test different Bernstein bands and compare
--Test cross-correlation b/t one band and another on single VLF channel. 
--Test aliased Bernstein waves with Skins. This has advantages for doing interferometry
--Plot difference in PSD for long and short booms (see Pfaff+96). Any systematic difference
    can indicate that the shorter booms are not responding the same.
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
#import end_load_data as end
#import end_load_gainphase as gainphase

from scipy.interpolate import interp1d
#from scipy.fft import rfft, irfft
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL





#---------------------------------------------
#Load VLF gain/phase corrected data
#---------------------------------------------

vDC = EFL('V12D')
wfDC, tdatDC = vDC.load_data_gainphase_corrected()
fsDC = vDC.chnspecs['fs']

v12 = EFL('VLF12D')
v34 = EFL('VLF34D')
#v13 = EFL('VLF13D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')
#v41 = EFL('VLF41D')

wf12, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']
wf34, tgoo = v34.load_data_gainphase_corrected()
#wf13, tgoo = v13.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
#wf41, tgoo = v41.load_data_gainphase_corrected()





"""
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

"""


#----------------------------------------------------------------------
#Get spectral data for finding waves (use VLF12)
#----------------------------------------------------------------------

fspec, tspec, powerc = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powercDC = signal.spectrogram(wfDC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')


#----------------------------------------------------------------------
#Compare all VLF waveforms for desired timerange
#----------------------------------------------------------------------

def plot_wf(tr):

    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wf12bp[goot]), list(wf34bp[goot]), list(wf24bp[goot]), list(wf32bp[goot])])


    fig,axs = plt.subplots(6,figsize=(8,8))
    axs[0].plot(tdat[goot],wf12bp[goot],color='b')
    axs[0].plot(tdat[goot],wf34bp[goot],color='g')
    axs[0].title.set_text('VLF (12 vs 34)')
    axs[1].plot(tdat[goot],wf12bp[goot],color='b')
    axs[1].plot(tdat[goot],wf24bp[goot],color='r')
    axs[1].title.set_text('VLF (12 vs 24)')
    axs[2].plot(tdat[goot],wf12bp[goot],color='b')
    axs[2].plot(tdat[goot],wf32bp[goot],color='c')
    axs[2].title.set_text('VLF (12 vs 32)')
    axs[3].plot(tdat[goot],wf34bp[goot],color='g')
    axs[3].plot(tdat[goot],wf24bp[goot],color='r')
    axs[3].title.set_text('VLF (34 vs 24)')
    axs[4].plot(tdat[goot],wf34bp[goot],color='g')
    axs[4].plot(tdat[goot],wf32bp[goot],color='c')
    axs[4].title.set_text('VLF (34 vs 32)')
    axs[5].plot(tdat[goot],wf24bp[goot],color='r')
    axs[5].plot(tdat[goot],wf32bp[goot],color='c')
    axs[5].title.set_text('VLF (24 vs 32)')
    for i in range(6):
        axs[i].set_ylim(-maxv, maxv)


#----------------------------------------------------------------------
#Compare all hodograms for desired timerange
#----------------------------------------------------------------------

def plot_hod(tr):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wf12bp[goot]), list(wf34bp[goot]), list(wf24bp[goot]), list(wf32bp[goot])])

    fig,axs = plt.subplots(3,2,figsize=(8,8))
    axs[0,0].plot(wf12bp[goot],wf34bp[goot])    
    axs[0,1].plot(wf12bp[goot],wf24bp[goot])    
    axs[1,0].plot(wf12bp[goot],wf32bp[goot])    
    axs[1,1].plot(wf34bp[goot],wf24bp[goot])    
    axs[2,0].plot(wf34bp[goot],wf32bp[goot])    
    axs[2,1].plot(wf24bp[goot],wf32bp[goot])    
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
def plot_csd(tr, xlm=[0,12000],ylm=[-180,180],nps=512,mincoh=0.5):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    window = np.hanning(len(goot[0]))
    wf12z = wf12[goot]*window
    wf24z = wf24[goot]*window
    wf34z = wf34[goot]*window
    wf32z = wf32[goot]*window
    csd1234, angle1234, f = correlation_analysis.cross_spectral_density(wf12z,wf34z,fs,nperseg=nps)
    csd1224, angle1224, f = correlation_analysis.cross_spectral_density(wf12z,wf24z,fs,nperseg=nps)
    csd1232, angle1232, f = correlation_analysis.cross_spectral_density(wf12z,wf32z,fs,nperseg=nps)
    csd3424, angle3424, f = correlation_analysis.cross_spectral_density(wf34z,wf24z,fs,nperseg=nps)
    csd3432, angle3432, f = correlation_analysis.cross_spectral_density(wf34z,wf32z,fs,nperseg=nps)
    csd2432, angle2432, f = correlation_analysis.cross_spectral_density(wf24z,wf32z,fs,nperseg=nps)

    goodv1234 = np.where(csd1234 >= mincoh)
    goodv1224 = np.where(csd1224 >= mincoh)
    goodv1232 = np.where(csd1232 >= mincoh)
    goodv3424 = np.where(csd3424 >= mincoh)
    goodv3432 = np.where(csd3432 >= mincoh)
    goodv2432 = np.where(csd2432 >= mincoh)


    fig, axs = plt.subplots(4,3,figsize=(12,8))
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)

    axs[0,0].plot(f, csd1234)
    axs[0,0].plot(f[goodv1234], csd1234[goodv1234],'.')
    axs[0,0].set_ylabel('coherency\nVLF12-VLF34')
    axs[0,0].set_xlabel('freq(Hz)')
    axs[1,0].plot(f, angle1234, color='r')
    axs[1,0].plot(f[goodv1234], angle1234[goodv1234],'.')
    axs[1,0].set_ylim(ylm)
    axs[1,0].set_ylabel('Phase\nVLF12-VLF34')

    axs[0,1].plot(f, csd1224)
    axs[0,1].plot(f[goodv1224], csd1224[goodv1224],'.')
    axs[0,1].set_ylabel('coherency\nVLF12-VLF24')
    axs[0,1].set_xlabel('freq(Hz)')
    axs[1,1].plot(f, angle1224, color='r')
    axs[1,1].plot(f[goodv1224], angle1224[goodv1224],'.')
    axs[1,1].set_ylim(ylm)
    axs[1,1].set_ylabel('Phase\nVLF12-VLF24')

    axs[0,2].plot(f, csd1232)
    axs[0,2].plot(f[goodv1232], csd1232[goodv1232],'.')
    axs[0,2].set_ylabel('coherency\nVLF12-VLF32')
    axs[0,2].set_xlabel('freq(Hz)')
    axs[1,2].plot(f, angle1232, color='r')
    axs[1,2].plot(f[goodv1232], angle1232[goodv1232],'.')
    axs[1,2].set_ylim(ylm)
    axs[1,2].set_ylabel('Phase\nVLF12-VLF32')

    axs[2,0].plot(f, csd3424)
    axs[2,0].plot(f[goodv3424], csd3424[goodv3424],'.')
    axs[2,0].set_ylabel('coherency\nVLF34-VLF24')
    axs[2,0].set_xlabel('freq(Hz)')
    axs[3,0].plot(f, angle3424, color='r')
    axs[3,0].plot(f[goodv3424], angle3424[goodv3424],'.')
    axs[3,0].set_ylim(ylm)
    axs[3,0].set_ylabel('Phase\nVLF34-VLF24')

    axs[2,1].plot(f, csd3432)
    axs[2,1].plot(f[goodv3432], csd3432[goodv3432],'.')
    axs[2,1].set_ylabel('coherency\nVLF34-VLF32')
    axs[2,1].set_xlabel('freq(Hz)')
    axs[3,1].plot(f, angle3432, color='r')
    axs[3,1].plot(f[goodv3432], angle3432[goodv3432],'.')
    axs[3,1].set_ylim(ylm)
    axs[3,1].set_ylabel('Phase\nVLF34-VLF32')

    axs[2,2].plot(f, csd2432)
    axs[2,2].plot(f[goodv2432], csd2432[goodv2432],'.')
    axs[2,2].set_ylabel('coherency\nVLF24-VLF32')
    axs[2,2].set_xlabel('freq(Hz)')
    axs[3,2].plot(f, angle2432, color='r')
    axs[3,2].plot(f[goodv2432], angle2432[goodv2432],'.')
    axs[3,2].set_ylim(ylm)
    axs[3,2].set_ylabel('Phase\nVLF24-VLF32')


    for i in range(4):
        for j in range(3):
            axs[i,j].set_xlim(xlm)
            #axs[i,j].set_ylim(ylm)



#--------------------------------------
#Bernstein vertical spikes
#--------------------------------------
fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-30],yr=[4000,8000],xr=[110,130], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,0],yr=[0,100],xr=[110,130], yscale='linear',ax=axs[1])

wfDCbp = filt.butter_highpass_filter(wfDC,1,fsDC,order=10)
wf12bp = filt.butter_highpass_filter(wf12,4000,fs,order=10)

tr = [117, 120] 
fig, axs = plt.subplots(2)
axs[0].plot(tdatDC,wfDCbp)
axs[1].plot(tdat,wf12bp)
for i in range(2): 
    axs[i].set_xlim(tr)

axs[0].set_ylim(-1,2.5)
axs[1].set_ylim(-0.4,0.4)


#--------------------------------------
#Wave 1: Berstein (~120-130 sec at 5000-8000 Hz)
#--------------------------------------

ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='linear')


fmin = 4000
fmax = 6000
wf12bp = filt.butter_bandpass_filter(wf12, fmin, fmax, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, fmin, fmax, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, fmin, fmax, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, fmin, fmax, fs, order= 10)

fspecz, tspecz, powercz = signal.spectrogram(wf12bp, fs, nperseg=4*2048,noverlap=2*2048,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspecz,fspecz,np.abs(powercz),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='log')



#Snapshot 1
tr = [122.62,122.72]  #zoomed out
#tr = [122.67718,122.68]  #zoomed in
plot_wf(tr)

#tr = [122.62,122.72]  #zoomed out
tr = [124,130]  #zoomed out
S34, f = correlation_analysis.psd(wf34bp, tdat, fs, tr)
S32, f = correlation_analysis.psd(wf32bp, tdat, fs, tr)

plt.semilogy(f,S34, f, S32)
plt.ylim(1e-12, 1e-4)
#plt.xscale('log')
plt.xscale('linear')
plt.xlim(10,12000)
plt.xlim(4000,8000)


#import scipy.signal
#(f, S) = scipy.signal.periodogram(wf, fs, scaling='density')


plot_csd([tr[0]-3, tr[1]+3],nps=256)
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
wf1z = wf12[goot]*window
wf2z = wf34[goot]*window
timechunk = 1  #sec
Pxy, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf1z,wf2z,tdat[goot],fs,timechunk,nperseg=1024,plot=True)






print("here")











#--------------------------------------
#Wave 3: 165.46 - 165.56 sec; ~100-1000 Hz
#--------------------------------------

wf12bp = filt.butter_bandpass_filter(wf12, 60, 300, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, 60, 300, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, 60, 300, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, 60, 300, fs, order= 10)

plot_wf([165.20,165.8])
plot_hod([165.20,165.8])

