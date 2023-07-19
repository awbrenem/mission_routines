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
from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader



#--------------------------------------------------------------
#Get timeline of data to separate out science collection times 
#--------------------------------------------------------------

tl, gsS, gsE = end_data_loader.load_timeline()


#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

v12DC = EFL('V12D')
v34DC = EFL('V34D')
wf12DC, tdat12DC = v12DC.load_data_gainphase_corrected()
wf34DC, tdat34DC = v34DC.load_data_gainphase_corrected()
fsDC = v12DC.chnspecs['fs']

v12 = EFL('VLF12D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')

wf12, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()

v1 = EFL('V1SD')
v2 = EFL('V2SD')
v3 = EFL('V3SD')
v4 = EFL('V4SD')

wf1, tdats = v1.load_data_gainphase_corrected()
wf2, tdats = v2.load_data_gainphase_corrected()
wf3, tdats = v3.load_data_gainphase_corrected()
wf4, tdats = v4.load_data_gainphase_corrected()
fss = v1.chnspecs['fs']

#----------------------------------------------------------------------
#Get spectral data for finding waves (use VLF12)
#----------------------------------------------------------------------

fspec, tspec, powerc = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs, tspecs, powercs = signal.spectrogram(wf1, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs2, tspecs2, powercs2 = signal.spectrogram(wf1-wf2, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec12DC, tspec12DC, powerc12DC = signal.spectrogram(wf12DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')


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



#Cross-spectral density for skins 
def plot_csd_skins(tr, xlm=[0,1000],ylm=[-180,180],nps=512,mincoh=0.5,ylmpsd=[0,0]):
    goot = np.where((tdats >= tr[0]) & (tdats <= tr[1]))
    window = np.hanning(len(goot[0]))
    wf1z = wf1[goot]*window
    wf2z = wf2[goot]*window
    wf3z = wf3[goot]*window
    wf4z = wf4[goot]*window
    
    csd12, angle12, f = correlation_analysis.cross_spectral_density(wf1z,wf2z,fss,nperseg=nps)
    #csd34, angle34, f = correlation_analysis.cross_spectral_density(wf3z,wf4z,fss,nperseg=nps)
 
    goodv12 = np.where(csd12 >= mincoh)
    #goodv34 = np.where(csd34 >= mincoh)


    psdS1, psdf1 = correlation_analysis.psd(wf1z, tdats[goot], fss, tr)
    psdS2, psdf2 = correlation_analysis.psd(wf2z, tdats[goot], fss, tr)


    fig, axs = plt.subplots(3,figsize=(12,8))
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)

    axs[0].plot(f, csd12)
    axs[0].plot(f[goodv12], csd12[goodv12],'.')
    axs[0].set_ylabel('coherency\nV1,V2')
    axs[0].set_ylim(0,1)
    axs[0].set_xlabel('freq(Hz)')
    axs[1].plot(f, angle12, color='r')
    axs[1].plot(f[goodv12], angle12[goodv12],'.')
    axs[1].set_ylim(ylm)
    axs[1].set_ylabel('Phase\nV1-V2')
    axs[2].plot(psdf1, psdS1, psdf2, psdS2)
    axs[2].set_ylabel('PSD')
    axs[2].set_yscale('log')

    for i in range(3):
        axs[i].set_xlim(xlm)
    
    if ylmpsd != [0,0]:
        axs[2].set_ylim(ylmpsd)




#Cross-spectral density for VLF
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





"""
#--wave 1 possibility (small amp)
plt.plot(tdats[good],w1,tdats[good],w2)
plt.xlim(397,398)
plt.ylim(-0.001,0.001)
plot_csd_skins([395,398])


#--wave 2 possibility (larger amp than wave 1. Roughly 25 Hz blip)
plt.plot(tdats[good],w1,tdats[good],w2)
#plt.xlim(651.5,654)
#plt.xlim(650.2,653.5)
plt.xlim(649.2,654)
plt.ylim(-0.005,0.005)
plot_csd_skins([649.7,654],mincoh=0.6,xlm=[0,200],nps=1024,ylmpsd=[1e-13,1e-8])



#--prior to density structures
plt.plot(tdats[good],w1,tdats[good],w2)
plt.xlim(390,392)
plt.ylim(-0.001,0.001)


#artificial signal 
plt.plot(tdats[good],w1,tdats[good],w2)
plt.xlim(399,400)
plt.ylim(-0.01,0.01)




fig,axs = plt.subplots(2)
axs[0].plot(tdats[good],w1)
axs[1].plot(tdats[good],w2)
"""

#----------------------------------------------------------------------------
#Plot coherence and phase 
#----------------------------------------------------------------------------


tchunk = 2  #sec
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.5,nperseg=2048)

fig,axs = plt.subplots(2)
ps.plot_spectrogram(tchunks,freqs,coh,vr=[0.9,1], zscale='linear',xr=[120,220],yr=[5500,7000],yscale='linear',ax=axs[0])
ps.plot_spectrogram(tchunks,freqs,np.abs(phase),vr=[0,140], zscale='linear',xr=[120,220],yr=[5500,7000],yscale='linear',ax=axs[1])


#Slice the phase vs freq for particular time(s)

tz = 150
goo = np.where(tchunks >= tz)

navg = 8  #number of consecutive bins to average phase over


pgoo = phase[:,goo[0][0]:goo[0][navg]]
pavg = [0] * np.shape(pgoo)[0]
for i in range(np.shape(pgoo)[0]):
    pavg[i] = np.average(pgoo[i,:])    


plt.plot(freqs,phase[:,goo[0][0:navg]])
plt.plot(freqs,pavg,'.',color='black')
plt.xlim(5000,7500)
plt.ylim(0,180)

#6040; 70.6
#6524; 125.8 

df = 484
dp = 55



fig,axs = plt.subplots(2)
ps.plot_spectrogram(tchunks,freqs,coh,vr=[0.7,1], zscale='linear',xr=[100,800],yr=[0,120],yscale='linear',ax=axs[0])
ps.plot_spectrogram(tchunks,freqs,np.abs(phase),vr=[0,180], zscale='linear',xr=[100,800],yr=[0,120],yscale='linear',ax=axs[1])



#--------------------------------------
#Wave TEST - long wavelength LH 
#--------------------------------------

#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-10],yr=[6000,12000],xr=[430,432], yscale='linear')
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[6000,12000],xr=[800,810], yscale='linear')
fmin = 6000
fmax = 12000
wf12bp = filt.butter_bandpass_filter(wf12, fmin, fmax, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, fmin, fmax, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, fmin, fmax, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, fmin, fmax, fs, order= 10)


tr = [430.8222,430.825] 
#tr = [802.003,802.006] 
plot_wf(tr)



#--------------------------------------
40#Wave 1: Berstein (~120-130 sec at 5000-8000 Hz)
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

plt.plot(tdat,wf12bp,tdat,-wf24bp)
plt.xlim(122.67775,122.6800)
plt.ylim(-0.05,0.05)



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



#-------------------------------------------------------------------------------------
#Density structures
#--can be somewhat difficult to identify in skin data
#-------------------------------------------------------------------------------------

ps.plot_spectrogram(tspec12DC,fspec12DC,np.abs(powerc12DC),vr=[-60,10],yr=[10,1000],xr=[100,850], yscale='log')


fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[10,100],xr=[630,660], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecs2,fspecs2,np.abs(powercs2),vr=[-60,-40],yr=[10,100],xr=[630,660], yscale='linear',ax=axs[1])
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-10],yr=[10,1000],xr=[300,340], yscale='log')
#ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[10,1000],xr=[300,340], yscale='log')

#tr = [630,660]
tr = [100,300]
good = np.where((tdats >= tr[0]) & (tdats <= tr[1]))
goodDC = np.where((tdat12DC >= tr[0]) & (tdat12DC <= tr[1]))


w1 = filt.butter_bandpass_filter(wf1[good], 2, 80, fss, order= 10)
w2 = filt.butter_bandpass_filter(wf2[good], 2, 80, fss, order= 10)
w3 = filt.butter_bandpass_filter(wf3[good], 2, 80, fss, order= 10)
w4 = filt.butter_bandpass_filter(wf4[good], 2, 80, fss, order= 10)
w12DC = filt.butter_bandpass_filter(wf12DC[goodDC], 2, 80, fsDC, order= 10)
w34DC = filt.butter_bandpass_filter(wf34DC[goodDC], 2, 80, fsDC, order= 10)

"""
trz = [651,655]
#trz = [652.5,653] #25 Hz wiggles
fig, axs = plt.subplots(4)
ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[10,200],xr=trz, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecs2,fspecs2,np.abs(powercs2),vr=[-60,-40],yr=[10,200],xr=trz, yscale='linear',ax=axs[1])
axs[2].plot(tdats[good],w1)
axs[2].set_xlim(trz)
axs[2].set_ylim(-0.003,0.003)
axs[3].plot(tdats[good],w1-w2)
axs[3].set_xlim(trz)
axs[3].set_ylim(-0.003,0.003)
"""


fig,axs = plt.subplots(2)
axs[0].plot(tdat12DC[goodDC],w12DC)
axs[1].plot(tdat34DC[goodDC],w34DC)
for i in range(6): axs[i].set_xlim(150,151)
for i in range(2): axs[i].set_ylim(-0.001,0.001)
for i in range(2): axs[i+2].set_ylim(-0.0006,0.0006)
axs[4].set_ylim(-0.4,0.4)
axs[5].set_ylim(-0.4,0.4)



fig,axs = plt.subplots(6)
axs[0].plot(tdats[good],w1,tdats[good],w2)
axs[1].plot(tdats[good],w3,tdats[good],w4)
axs[2].plot(tdats[good],w1-w2)
axs[3].plot(tdats[good],w3-w4)
axs[4].plot(tdat12DC[goodDC],w12DC)
axs[5].plot(tdat34DC[goodDC],w34DC)
for i in range(6): axs[i].set_xlim(150,151)
for i in range(2): axs[i].set_ylim(-0.001,0.001)
for i in range(2): axs[i+2].set_ylim(-0.0006,0.0006)
axs[4].set_ylim(-0.4,0.4)
axs[5].set_ylim(-0.4,0.4)

#Density structure events 

#150-170 - lots of few Hz structuring
#166.4 (12)

#212.5 (12-34)  - large
#240-260 (12-34) - lots of 1/3 Hz fluctuations
#260.8
#268.6 (12-34)
#286.4 (34)

#304.75 (12)
#308.1 (12-34)
#312.75 (34)
#338.75 (12)
#376.7 (12-34)
#395.5 (12)
#396.8 (12-34)
#397.3 (12)

#400-600 (intermittent ~quarter to half second turbulence - e.g. 430-434)
#647.2-647.6 (12)
#648.4 (34)
#650.2 (12-34)
#652.7 (12-34)
#692.3 (12-34)
#705.5 (12-34)
#706.5 (12-34)
#722.8 (12,34)




#Timing analysis on density structure (**only good one I've found so far**)
tshift = 0.003
plt.plot(tdats[good],w1,tdats[good]-tshift,-w2)
plt.xlim(652.4,653.0)
plt.ylim(-0.002,0.002)


plt.plot(tdats[good],w3,tdats[good],-w4)
plt.xlim(648.0,648.5)
plt.ylim(-0.001,0.001)





"""
#-----------------------------------
#Artifical signal with harmonics that slowly rise in tone over entire mission
#-----------------------------------

ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[80,200],xr=[380,420], yscale='linear')

tr = [388,395]
good = np.where((tdats >= tr[0]) & (tdats <= tr[1]))


w1 = filt.butter_highpass_filter(wf1[good], 1, fss, order= 10)
w2 = filt.butter_highpass_filter(wf2[good], 1, fss, order= 10)
w3 = filt.butter_highpass_filter(wf3[good], 1, fss, order= 10)

w4 = filt.butter_highpass_filter(wf4[good], 1, fss, order= 10)



plt.plot(tdats[good],w1,tdats[good],w2,tdats[good],w3,tdats[good],w4)
plt.xlim(390.4,390.5)

fig,axs = plt.subplots(2)
axs[0].plot(tdats[good],w1)
axs[1].plot(tdats[good],w2)


plot_csd_skins([390,392])

"""
