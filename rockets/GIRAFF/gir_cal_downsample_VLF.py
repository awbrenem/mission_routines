
#Downsample the GIRAFF VLF data to be at a lower cadence for use with waves < 50 kHz (nearly all the interesting waves).
#Not doing this means the waveform files are HUGE and difficult to use. 


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
#sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from scipy import signal
import numpy as np 
import pickle



#Load E-fields data
pld = '380'
chn = 'VLF34D'
v12 = GFL(pld,chn)
wf12, tdat = v12.load_data()



#original and renamed files
if pld == '381':
    fnorig = '36' + pld + '_GIRAFF_TM2_LFDSP_IT_'+chn+'_mvm'
else:
    fnorig = '36' + pld + '_GIRAFF_TM2_LFDSP_'+chn+'_mvm'

fnnew = fnorig + '_downsampled_50kHz'



#Determine number of points to resample to. 
fs_max = 50000. #Max useable frequency 
nyq = 2*fs_max  #Nyquist
srt = nyq + nyq*0.1  #Sample rate - make 10% higher than Nyquist to avoid aliasing



#Resample in chunks so that we don't have to pass huge arrays (very slow).
chunksz = int(1e6)  #number of points for each chunk
nchunks = int(len(wf12) / chunksz)  
nsamp = int(np.round(srt * (tdat[int(chunksz)] - tdat[0]))) #samples that will be in each chunk after resampling 




#Downsample each chunk
wf = np.zeros((nsamp, nchunks))
tv = np.zeros((nsamp, nchunks))
for i in range(nchunks):
      y, x = signal.resample(wf12[i*chunksz:(i+1)*chunksz], nsamp, t=tdat[i*chunksz:(i+1)*chunksz], window='hann')
      wf[:,i] = y 
      tv[:,i] = x 




#Reshape the chunked arrays to an array of size [n]
wf = wf.reshape(nsamp*nchunks, order='F')
times = tv.reshape(nsamp*nchunks, order='F')



notes = {'VLF waveform data from ' + fnorig + ' downsampled to 50kHz from gir_cal_downsample_VLF.py'}


path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/efield_VLF/' + fnnew + '.pkl'
pickle.dump([times, wf, notes],open(path,'wb'))





"""
#Test the results
import plot_spectrogram as ps
import matplotlib.pyplot as plt

plt.plot(times)
plt.plot(times,wf)
#plt.plot(tdat,wf12)

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tv2, wf2, fs_thres=0.1, nfft=nfft, noverlap=2)
fspec = fspec[:,0]

vr=[-80,-50]
xr = [100,550]
yr = [0,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=[0,1400],xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='log',yr=[300,50000],xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")

"""


