"""
Code to plot the coherence and phase spectra for selected boom pairs

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt


#Load data for two channels of interest
c1s = 'VLF24D'
c2s = 'VLF13D'

c1 = EFL(c1s)
fs = c1.chnspecs['fs']
wf1, tdat = c1.load_data_gainphase_corrected()
c2 = EFL(c2s)
wf2, tdat = c2.load_data_gainphase_corrected()

#------see if this needs to be flipped
wf2 = -1 * wf2
cs2 = 'VLF31D'
#-------------------------------------

#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 4096
fspec, tspec, powerc1 = signal.spectrogram(wf1, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc2 = signal.spectrogram(wf2, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')

cohmin = 0.01  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result
#cohmin = 0.3 



##NOTE: + sense of phase defined as pointing towards center of potential of "powerc1"
##Nval = 3
#Nval = 3
#gx,cohx,phasex = correlation_analysis.interferometric_coherence_2D(powerc1,powerc2,Nval,coh_min=cohmin)
#phasex = np.degrees(phasex)


tchunk = 0.5  #delta-time (sec) for each time chunk to divide up the spectra into (and average over)
nchunks = int(np.ceil((wf1.size/fs)/tchunk)) #number of chunks in ENTIRE timerange
cohx2, phasex2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wf1,wf2,tdat,fs,tchunk,coh_min=cohmin,nperseg=512)


ptmp_diff = np.abs(np.abs(powerc1) - np.abs(powerc2))
ptmp_sum = np.abs(powerc1) + np.abs(powerc2)
ptmp_fracdiff = ptmp_diff/ptmp_sum



#Plot values from Method 2
#yr = [0,15000]
yr = [4000,8000]
yscale='linear'
xr = [100,900]
#xr = [100,300]
vr = [-45,-25]

fig,axs = plt.subplots(5)  #,gridspec_kw={'height_ratios':[1,1,1,1,1,1,1,1,1]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=vr,xr=xr,yr=yr,yscale=yscale,ax=axs[0],xlabel='time(sec)',ylabel=c1s + "\nfreq(Hz)")
axs[0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc2),vr=vr,xr=xr,yr=yr,yscale=yscale,ax=axs[1],xlabel='time(sec)',ylabel=c2s + "\nfreq(Hz)")
axs[1].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspec,fspec,ptmp_fracdiff,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=yscale,ax=axs[2],xlabel='time(sec)',ylabel="percent diff\nfreq(Hz)")
axs[2].get_xaxis().set_visible(False)
ps.plot_spectrogram(tchunks2,freqs2,phasex2,vr=[-180,180],cmap='twilight_shifted',zscale='linear',xr=xr,yr=yr,yscale=yscale,ax=axs[3],xlabel='time(sec)',ylabel='phase\n('+c1s+'-'+c2s+')\nfreq(Hz)')
axs[3].get_xaxis().set_visible(False)
ps.plot_spectrogram(tchunks2,freqs2,cohx2**2,vr=[0,1], zscale='linear',xr=xr,yr=yr,yscale=yscale,ax=axs[4],xlabel='time(sec)',ylabel='coh**2\n(coh >'+str(cohmin)+')\nfreq(Hz)')

#ps.plot_spectrogram(tspec,fspec,cohx**2,vr=[0,1], zscale='linear',xr=xr,yr=yr,yscale=yscale,ax=axs[3],xlabel='time(sec)',ylabel='coh**2\nfreq(Hz)')
#ps.plot_spectrogram(tspec,fspec,phasex,vr=[-180,180],cmap='turbo',zscale='linear',xr=xr,yr=yr,yscale=yscale,ax=axs[4],xlabel='time(sec)',ylabel='phase (c1-c2)\nfreq(Hz)')


good = np.where(tchunks2 > 820)

fig,axs = plt.subplots(2)
axs[0].plot(freqs2,phasex2[:,good[0][0]])
axs[0].plot(freqs2,phasex2[:,good[0][0]],'.')
#axs[0].set_xlim(5000,8000)
axs[0].set_xlim(0,8000)
axs[1].plot(freqs2,cohx2[:,good[0][0]])
axs[1].plot(freqs2,cohx2[:,good[0][0]],'.')
#axs[1].set_xlim(5000,8000)
axs[1].set_xlim(0,8000)




yr = [4500,8000]
xr = [700,900]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=vr,xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel=c1s + "\nfreq(Hz)")
