"""
Determine polarization of various waves from Endurance.

UPLEG: the magnetic field (Bo) is in the direction of V1 x V3 and rocket (main) velocity is 
in the V1 x V4 direction.

DOWNLEG: Bo in the same direction as rocket velocity (V1 x V4). 

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
#import end_load_data
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader

#import plasma_params_get_density_from_flhr_freq as dflh
import numpy as np 
import matplotlib.pyplot as plt
#import plasmapy
#from astropy import units as u  
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency
import plot_spectrogram as ps
import plot_hodogram_dynamic as hod
from scipy import signal




#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]

"""Enable auto module reloading"""
#%load_ext autoreload
#%autoreload 2


#evlf = end_load_data.efield_vlf()
#Index(['tsec', 'amp'], dtype='object')


v12 = EFL('VLF12D')
v34 = EFL('VLF34D')

fs = v12.chnspecs['fs']
wf12, tvals = v12.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()


#########
#Artificially shift wf12 one point to the right to test for anomalous instrument timing issue
#wf12 = np.roll(wf12,1)



#Spectral overview
nps = 2048
freq12, tspec12, power12 = signal.spectrogram(wf12, fs, nperseg=nps,noverlap=nps/2,window='hann') #, return_onesided=1)
freq34, tspec34, power34 = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann') #, return_onesided=1)
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40], xr=[900,1000],yr=[6000,10000], yscale='linear')
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[0,10000],xr=[100,900], yscale='linear')
ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[6000,14000],xr=[280,320], yscale='linear')



#Isolate waveform for polarization

#tr = [400, 405]
#tr = [146, 146.00015]
#fr = [5500, 6000]
#fr = [5900, 6100]
tr = [160, 180]
fr = [5000, 5900]

#tr = [301.5,301.65]
#fr = [6000,14000]

#-Left hand waves
#tr = [705, 715]
#fr = [5300,5600]


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
ps.plot_spectrogram(tspec12,freq12,p12,vr=[-80,-30], xr=tr,yr=fr, yscale='log')


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
ps.plot_spectrogram(tspec34,freq34,p34,vr=[-80,-30], xr=tr,yr=fr, yscale='log')




#Bandpass filter waveform to preferred freqs and times
filt12 = filter_wave_frequency.butter_bandpass_filter(wf12, fr[0],fr[1],fs,order=2)
filt34 = filter_wave_frequency.butter_bandpass_filter(wf34, fr[0],fr[1],fs,order=2)






#all the points to be plotted
trz= [165.00, 165.005]
#trz= [301.56, 301.58]
goodt = np.squeeze(np.where((tvals > trz[0]) & (tvals < trz[1])))
filt12z = filt12[goodt]
filt34z = filt34[goodt]

ptitle = 't=' + str(trz[0]) + '-' + str(trz[1]) + ' sec'
num = 1

plot_kwargs = {'xlim':[-0.1,0.1],
               'xlabel':'VLF12',
               'ylabel':'VLF34',
               'title':ptitle}


plt.plot(filt12z,filt34z)
plt.plot(np.roll(filt12z,-1),filt34z)


hod.plot_hodogram_dynamic(filt12z, filt34z, npts=num, gap=1, plot_kwargs=plot_kwargs)#,pauseT=0.3)

















print("here")





