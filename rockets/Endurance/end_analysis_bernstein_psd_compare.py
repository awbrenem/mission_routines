"""
Compare power spectral densities for the Bernstein waves for different VLF channels

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import filter_wave_frequency as filt


v12 = EFL('VLF12D')
v13 = EFL('VLF13D')
v41 = EFL('VLF41D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')

fs = v12.chnspecs['fs']


wf12, tdat = v12.load_data_gainphase_corrected()
wf13, tgoo = v13.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf41, tgoo = v41.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
wf42 = -wf24
wf14 = -wf41




#tr = [140,900]
#tr = [140,145]
#tr = [220,260]
#tr = [750,750.2]
tr = [720,725]
#tr = [725,760]
#tr = [380,385]
goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))

wf12z = wf12[goot]
wf13z = wf13[goot]
wf34z = wf34[goot]
wf42z = wf42[goot]
wf14z = wf14[goot]
wf32z = wf32[goot]


#Windowed (hanning) FFT
nfft = 1024
psd12, psdf = correlation_analysis.psd(wf12z, tdat[goot], fs, tr, nft=nfft)
psd13, psdf = correlation_analysis.psd(wf13z, tdat[goot], fs, tr, nft=nfft)
psd34, psdf = correlation_analysis.psd(wf34z, tdat[goot], fs, tr, nft=nfft)
psd42, psdf = correlation_analysis.psd(wf42z, tdat[goot], fs, tr, nft=nfft)
psd14, psdf = correlation_analysis.psd(wf14z, tdat[goot], fs, tr, nft=nfft)
psd32, psdf = correlation_analysis.psd(wf32z, tdat[goot], fs, tr, nft=nfft)


#All pairs compared
fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
axs[0].plot(psdf,psd12)  #blue
axs[0].set_yscale('linear')
axs[0].set_ylim(-90,-40)
axs[0].set_xlim(5000,7000)
axs[0].set_xscale('linear')
axs[0].plot(psdf,psd34) #orange
axs[0].plot(psdf,psd14) #green
axs[0].plot(psdf,psd13) #red
axs[0].plot(psdf,psd42) #purple
axs[0].plot(psdf,psd32) #brown



#version without V1
fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
axs[0].plot(psdf,psd34)  #blue
axs[0].set_yscale('linear')
axs[0].set_ylim(-80,-40)
axs[0].set_xlim(5000,8000)
axs[0].set_xscale('linear')
axs[0].plot(psdf,psd42) #orange
axs[0].plot(psdf,psd32) #green



#Comparison of V1 and V2
fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

#axs[0].plot(psdf,psd12,psdf,psd34)
axs[0].plot(psdf,psd12)  #blue
axs[0].set_yscale('linear')
axs[0].set_ylim(1e-8,4e-6)
axs[0].set_xlim(5000,8000)
axs[0].set_xscale('linear')
axs[0].plot(psdf,psd34) #orange













vr = [-45,-20]
ys = 'log'
yr = [1,15000]
xr = [350,450]
ftmp, ttmp, ptmp = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')

ftmp, ttmp, ptmp = signal.spectrogram(wf14, fs, nperseg=256,noverlap=256/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')


goo = np.where(ttmp >= tr[0])

plt.plot(ftmp,np.abs(ptmp[:,goo[0][0]]))



axs[2].plot(psdf1, psdS1, psdf2, psdS2)
axs[2].set_ylabel('PSD')
axs[2].set_yscale('log')

for i in range(3):
    axs[i].set_xlim(xlm)

if ylmpsd != [0,0]:
    axs[2].set_ylim(ylmpsd)

