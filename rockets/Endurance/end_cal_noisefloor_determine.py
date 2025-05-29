"""
Determine noise floor for each channel. 
RESULT: the diagonals have a higher noise floor due to shorter antenna length.

NOTE: this is not actually the instrument noise floor but the plasma noise. 
I need to get instrument noise floor from Paulo. 

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



nfft = 2048

#tr = [140,160]
#tr = [750,750.2]
#tr = [720,720.6]
#tr = [720,725]
#tr = [725,760]
#tr = [380,385]
#tr = [880,890]
#tr = [460,490]
#tr = [200,400]
#tr = [130,850]
tr = [500,505]
#tr = [41.25,41.252]

v12 = EFL('VLF12D')
v13 = EFL('VLF13D')
v41 = EFL('VLF41D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')

fs = v12.chnspecs['fs']



#--calibrated data
wf12, tdat = v12.load_data_gainphase_corrected()
wf13, tgoo = v13.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf41, tgoo = v41.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
wf42 = -wf24
wf14 = -wf41


goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
#window = np.hanning(len(goot[0]))



wf12z = wf12[goot]#*window
wf13z = wf13[goot]#*window
wf34z = wf34[goot]#*window
wf42z = wf42[goot]#*window
wf14z = wf14[goot]#*window
wf32z = wf32[goot]#*window


plt.plot(tdat[goot],wf12z)
plt.ylim(-3,3)
plt.xlim(300,400)
plt.plot(tdat[goot],wf34z)
plt.plot(tdat[goot],wf13z)


psd12, psdf = correlation_analysis.psd(wf12z, tdat[goot], fs, tr, nft=nfft)
psd13, psdf = correlation_analysis.psd(wf13z, tdat[goot], fs, tr, nft=nfft)
psd34, psdf = correlation_analysis.psd(wf34z, tdat[goot], fs, tr, nft=nfft)
psd42, psdf = correlation_analysis.psd(wf42z, tdat[goot], fs, tr, nft=nfft)
psd14, psdf = correlation_analysis.psd(wf14z, tdat[goot], fs, tr, nft=nfft)
psd32, psdf = correlation_analysis.psd(wf32z, tdat[goot], fs, tr, nft=nfft)

amp = 0.3 mV/m
amp2 = 0.00025 mV/m / sqrt(Hz)
df = 14.65  #bin size

amp22 = 0.00025 * df**2

#dt = 501.00077 - 501.000921
#6600 Hz
#amp = 0.1  mV/m
#psd amp = 2 uV/m 


#1.5e-7 mV/m / sqrt(Hz)

df = 14.6 

plt.plot(tdat[goot],wf12[goot])
plt.xlim(501,501.01)

#dt = 41.250177 - 41.250310
#amp = 13 uV/m
#@7500 Hz


plt.plot(psdf,psd12)
plt.ylim(-0.001,0.001)

v = plt.magnitude_spectrum(wf12z, Fs=fs, scale='dB')
Pxx, freqs = plt.psd(wf12z, Fs=fs, scale_by_freq=False)

mvm = 10**(-60/20)
uvm = 1000*mvm




#--uncalibrated data
wf12o, tdat = v12.load_data()
wf13o, tgoo = v13.load_data()
wf34o, tgoo = v34.load_data()
wf24o, tgoo = v24.load_data()
wf41o, tgoo = v41.load_data()
wf32o, tgoo = v32.load_data()
wf42o = -wf24o
wf14o = -wf41o

wf12oz = wf12o[goot]#*window
wf13oz = wf13o[goot]#*window
wf34oz = wf34o[goot]#*window
wf42oz = wf42o[goot]#*window
wf14oz = wf14o[goot]#*window
wf32oz = wf32o[goot]#*window


psd12o, psdf = correlation_analysis.psd(wf12oz, tdat[goot], fs, tr, nft=nfft)
psd13o, psdf = correlation_analysis.psd(wf13oz, tdat[goot], fs, tr, nft=nfft)
psd34o, psdf = correlation_analysis.psd(wf34oz, tdat[goot], fs, tr, nft=nfft)
psd42o, psdf = correlation_analysis.psd(wf42oz, tdat[goot], fs, tr, nft=nfft)
psd14o, psdf = correlation_analysis.psd(wf14oz, tdat[goot], fs, tr, nft=nfft)
psd32o, psdf = correlation_analysis.psd(wf32oz, tdat[goot], fs, tr, nft=nfft)


vr = [-50,-20]
yr = [4000,4500]
xr = [475,525]
fspec, tspec, powspec = signal.spectrogram(wf12, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
fspec, tspec, powspec = signal.spectrogram(wf34, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
fspec, tspec, powspec = signal.spectrogram(wf14, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
fspec, tspec, powspec = signal.spectrogram(wf32, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")


fspec, tspec, powspec = signal.spectrogram(wf42, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
fspec, tspec, powspec = signal.spectrogram(wf13, fs, nperseg=nfft,noverlap=nfft/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powspec),xr=xr,vr=vr,yr=yr,yscale='linear',xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")



#--plot calibrated data
fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[0].plot(psdf,np.sqrt(psd12))
axs[0].set_yscale('log')
axs[0].set_xlim(0,8000)
axs[0].set_xscale('linear')
axs[0].set_ylim(1e-9,1e-6)

axs[0].plot(psdf,np.sqrt(psd34))
axs[0].plot(psdf,np.sqrt(psd14))
axs[0].plot(psdf,np.sqrt(psd13))
axs[0].plot(psdf,np.sqrt(psd42))
axs[0].plot(psdf,np.sqrt(psd32))


#--plot uncalibrated data
axs[1].plot(psdf,np.sqrt(psd12o))
axs[1].set_yscale('log')
axs[1].set_xlim(10,15000)
axs[1].set_xscale('log')
axs[1].set_ylim(1e-9,1e-4)

axs[1].plot(psdf,np.sqrt(psd34o))
axs[1].plot(psdf,np.sqrt(psd14o))
axs[1].plot(psdf,np.sqrt(psd13o))
axs[1].plot(psdf,np.sqrt(psd42o))
axs[1].plot(psdf,np.sqrt(psd32o))


print('h')



#Determine how to adjust the amplitude offset so that the channels match in gain
fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[0].plot(psdf,np.sqrt(psd12))
axs[0].set_yscale('linear')
axs[0].set_xlim(10,100)
axs[0].set_xscale('log')
axs[0].set_ylim(1e-5,1e-4)

axs[0].plot(psdf,np.sqrt(psd34))
axs[0].plot(psdf,np.sqrt(psd14))
axs[0].plot(psdf,np.sqrt(psd13))
axs[0].plot(psdf,np.sqrt(psd42))
axs[0].plot(psdf,np.sqrt(psd32))






