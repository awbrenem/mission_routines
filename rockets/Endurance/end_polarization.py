"""
Determine polarization of various waves from Endurance.
"""



import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
import end_load_data
#import plasma_params_get_density_from_flhr_freq as dflh
import numpy as np 
import matplotlib.pyplot as plt
#import plasmapy
from astropy import units as u  
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
from scipy import signal



from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y




#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]

"""Enable auto module reloading"""
#%load_ext autoreload
#%autoreload 2



mag = end_load_data.mag_dc()
#Index(['tsec', 'Bx', 'By', 'Bz', 'Bmag'], dtype='object')
evlf = end_load_data.efield_vlf()
#Index(['tsec', 'amp'], dtype='object')
ehf = end_load_data.efield_hf()
#dict_keys(['afftpow12', 'afftpow34', 'afreq', 'atimesfft', 'fftsize', 'overlap', 'weight', 'samplerate', 'nfreq', 'nfftlines', 'in_file'])


#VLF data
tvals = evlf.tvlf
wf12 = evlf.dvlf12_mvm
wf34 = evlf.dvlf34_mvm
##HF data
#tvals = evlf.tvlf
#wf12 = evlf.dvlf12_mvm
#wf34 = evlf.dvlf12_mvm



#Spectral overview
fs = evlf.samplerate
freq12, tspec12, power12 = signal.spectrogram(wf12, fs, nperseg=16384,noverlap=16384/2) #, return_onesided=1)
freq34, tspec34, power34 = signal.spectrogram(wf34, fs, nperseg=16384,noverlap=16384/2) #, return_onesided=1)
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40], xr=[900,1000],yr=[6000,10000], yscale='linear')
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[0,10000],xr=[100,900], yscale='linear')
ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[10,1000],xr=[100,900], yscale='log')

#ps.plot_spectrogram(tspec34,freq34,power34,vr=[-100,-40], xr=[0,900],yr=[0,10000], yscale='linear')


#Isolate waveform for polarization

#tr = [400, 405]
tr = [146, 146.00015]
fr = [6100, 6300]
tr = [470,00]
#fr = [5000, 5900]
#fr = [6500, 6800]

pmask = np.full_like(power12, 1)

badt = np.squeeze(np.where(tspec12 < tr[0]))
if len(badt) != 0: pmask[:,badt] = 0
badt = np.squeeze(np.where(tspec12 > tr[1]))
if len(badt) != 0: pmask[:,badt] = 0
badf = np.squeeze(np.where((freq12 < fr[0])))
if len(badt) != 0: pmask[badf,:] = 0
badf = np.squeeze(np.where((freq12 > fr[1])))
if len(badt) != 0: pmask[badf,:] = 0

p12 = power12 * pmask

#ps.plot_spectrogram(tspec12,freq12,p12,vr=[-100,-40], xr=tr,yr=fr, yscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12,vr=[-100,-40], xr=tr,yr=fr, yscale='log')


pmask = np.full_like(power34, 1)

badt = np.squeeze(np.where(tspec34 < tr[0]))
pmask[:,badt] = 0
badt = np.squeeze(np.where(tspec34 > tr[1]))
pmask[:,badt] = 0
badf = np.squeeze(np.where((freq34 < fr[0])))
pmask[badf,:] = 0
badf = np.squeeze(np.where((freq34 > fr[1])))
pmask[badf,:] = 0

p34 = power34 * pmask

#ps.plot_spectrogram(tspec34,freq34,p34,vr=[-100,-40], xr=tr,yr=fr, yscale='linear')
ps.plot_spectrogram(tspec34,freq34,p34,vr=[-100,-40], xr=tr,yr=fr, yscale='log')




#Bandpass filter waveform to preferred freqs and times
filt12 = butter_bandpass_filter(wf12, fr[0],fr[1],fs,order=2)
filt34 = butter_bandpass_filter(wf34, fr[0],fr[1],fs,order=2)

#reduce to desired times only
goodt = np.squeeze(np.where((tvals > tr[0]) & (tvals < tr[1])))
filt12 = filt12[goodt]
filt34 = filt34[goodt]


#Plot hodograms 
fig, axs = plt.subplots()
axs.scatter(filt12,filt34)
#maxv = np.max(filt12) gt np.max(filt34)
maxv = max(np.max(filt12), np.max(filt34))
axs.set_xlim(-1*maxv, maxv)
axs.set_ylim(-1*maxv, maxv)
axs.set_aspect(1)




print("here")




