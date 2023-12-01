"""
VLF channels involving V1 may have an issue with extra noise during the Bernstein waves. 

The gain calibrations for the VLF channels are all within about 1-3% of each other 
at both DC and ~6 kHz (the gain curves are ideally flat from ~30 to 7 kHz)
But, the power in VLF channels involving V1 is up to 1.6x higher. 
There's also spurious noise in at the upper end of the Bernstein range. 

Test to see if this is in all the VLF channels involving V1. 


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
#wf12, tdat = v12.load_data()
#wf13, tgoo = v13.load_data()
#wf34, tgoo = v34.load_data()
#wf24, tgoo = v24.load_data()
#wf41, tgoo = v41.load_data()
#wf32, tgoo = v32.load_data()
wf42 = -wf24
wf14 = -wf41


fs = 30000



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 1024
fspec, tspec, powerc12 = signal.spectrogram(wf12, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc13 = signal.spectrogram(wf13, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc14 = signal.spectrogram(wf14, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc34 = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc42 = signal.spectrogram(wf42, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc32 = signal.spectrogram(wf32, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')


vr = [-40,-20]
xr = [140,180]
yr = [5000,8000]

#fig,axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF12')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc13),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF13')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc14),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF14')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF34')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc42),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF42')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc32),vr=vr,xr=xr,yr=yr,yscale='linear',title='VLF32')



tr = [140,160]
goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))

wf12z = wf12[goot]
wf13z = wf13[goot]
wf14z = wf14[goot]
wf34z = wf34[goot]
wf42z = wf42[goot]
wf32z = wf32[goot]

nfft = 1024

psd12, psdf = correlation_analysis.psd(wf12z, tdat[goot], fs, tr, nft=nfft)
psd13, psdf = correlation_analysis.psd(wf13z, tdat[goot], fs, tr, nft=nfft)
psd14, psdf = correlation_analysis.psd(wf14z, tdat[goot], fs, tr, nft=nfft)
psd34, psdf = correlation_analysis.psd(wf34z, tdat[goot], fs, tr, nft=nfft)
psd42, psdf = correlation_analysis.psd(wf42z, tdat[goot], fs, tr, nft=nfft)
psd32, psdf = correlation_analysis.psd(wf32z, tdat[goot], fs, tr, nft=nfft)


fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
axs[0].plot(psdf,np.sqrt(psd12))  #blue
axs[0].set_yscale('log')
axs[0].set_ylim(1e-10,4e-6)
axs[0].set_xlim(5000,8000)
axs[0].set_xscale('linear')

axs[0].plot(psdf,np.sqrt(psd34)) #orange
axs[0].plot(psdf,np.sqrt(psd14)) #green
axs[0].plot(psdf,np.sqrt(psd13)) #red
axs[0].plot(psdf,np.sqrt(psd42)) #purple
axs[0].plot(psdf,np.sqrt(psd32)) #brown



