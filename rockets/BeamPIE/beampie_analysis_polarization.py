"""
Determine polarization of triggered waves from the BeamPIE sounding rocket

Analyze swooshes from 280-317 sec (5-10 kHz)

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')


#import plasma_params_get_density_from_flhr_freq as dflh
import numpy as np 
import matplotlib.pyplot as plt
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency
import plot_spectrogram as ps
import plot_hodogram_dynamic as hod
from scipy import signal

#from os.path import dirname, join as pjoin
#import scipy.io as sio
from scipy.io import readsav



#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]



path = '/Users/abrenema/Desktop/Research/Rocket_missions/BeamPIE/data/bfield_AC/'
fn = '52009_TM2_LFDSP_S1-DeltaBxyz_counts.sav'
sav_data = readsav(path+fn)

bx = sav_data['dmagxcounts'] #on 3-4 axis near V3
by = sav_data['dmagycounts'] #on 1-2 axis near V1
bz = sav_data['dmagzcounts'] #on spin axis near V11
btimes = sav_data['tmag']

#;Subtract baseline offset
bx = bx - 16384.
by = by - 16384.
bz = bz - 16384.

wf12 = by 
wf34 = bx 
wfz = bz
times = btimes


#Orient Bw data these so that the antenna coord system points along Bo.
#See: --boom_diagram.jpg
#This is a RH coord. system with Bo along the z-direction.
wf12 *= -1 #x-hat
wfz *=  1  #z-hat
wf34 *= 1  #y-hat


#wf12 = sav_data12['data_vlf12d']
#wf34 = sav_data34['data_vlf34a']
#times12 = sav_data12['time_vlf12d']
#times34 = sav_data34['time_vlf34a']

#Remove values at negative times
wf12 = wf12[times >= 0]
wf34 = wf34[times >= 0]
times = times[times >= 0]

fs = 64000.

#tdelta = times - np.roll(times,1)


#Spectral overview
nps = 256
freq12, tspec12, power12 = signal.spectrogram(wf12, fs, nperseg=nps,noverlap=nps/2,window='hann') #, return_onesided=1)
freq34, tspec34, power34 = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann') #, return_onesided=1)
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40], xr=[900,1000],yr=[6000,10000], yscale='linear')
#ps.plot_spectrogram(tspec12,freq12,power12,vr=[-80,-40],yr=[0,10000],xr=[100,900], yscale='linear')

p12log = np.log10(power12)


#SWOOP times
tswoop = [261.,263.75,263.95,265.72,267.6,269.3,271.1,272.93,274.77,276.55,278.35,280.15,281.95,283.7,285.5,287.35,289.1,290.92,292.79,294.52,296.3,298.0,299.88,301.6,303.4,305.2,307.0]
tdelta = tswoop - np.roll(tswoop,1)

i=27
fig,axs = plt.subplots(2)
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[20,32000],xr=[tswoop[i]-0.17,tswoop[i]+0.17],vr=[-4,-1], yscale='linear', zscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,10000],xr=[tswoop[i]-0.17,tswoop[i]+0.17],vr=[-4,-1], yscale='linear', zscale='linear',ax=axs[1])


ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[263.6,264.2],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[265.5,266],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[267.,268.],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[269.,269.5],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[270.75,271.5],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[272.5,273.2],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[274.6,274.9],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[276.4,276.75],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[278.2,278.6],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[280.,280.3],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[281.5,282.2],vr=[-4,-1], yscale='linear', zscale='linear')
ps.plot_spectrogram(tspec12,freq12,p12log,yr=[2000,32000],xr=[250.0,270.],vr=[-4,-2], yscale='linear', zscale='linear')





#Isolate waveform for polarization

#tr = [629.5, 631]
#fr = [3250, 3450]
#tr = [265.65, 265.7]
tr = [265.75, 265.85]
fr = [4000, 8000]



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
ps.plot_spectrogram(tspec12,freq12,p12,vr=[-40,0], xr=tr,yr=fr, yscale='log')


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
ps.plot_spectrogram(tspec34,freq34,p34,vr=[-40,0], xr=tr,yr=fr, yscale='log')




#Bandpass filter waveform to preferred freqs and times
filt12 = filter_wave_frequency.butter_bandpass_filter(wf12, fr[0],fr[1],fs,order=2)
filt34 = filter_wave_frequency.butter_bandpass_filter(wf34, fr[0],fr[1],fs,order=2)






#all the points to be plotted
#trz= [629.60, 629.605]
#trz = [265.75, 265.8]
#trz = [265.65, 265.68]
trz = [265.75, 265.85]

goodt = np.squeeze(np.where((times > trz[0]) & (times < trz[1])))
filt12z = filt12[goodt]
filt34z = filt34[goodt]

ptitle = 't=' + str(trz[0]) + '-' + str(trz[1]) + ' sec\nf='+str(fr[0]) + '-' + str(fr[1]) + 'Hz'
num = 1


plt.plot(filt12z,filt34z)
plt.ylim(-8,8)
plt.xlim(-8,8)
#plt.plot(np.roll(filt12z,-1),filt34z)

plot_kwargs = {'xlim':[-8,8],
               'xlabel':'VLF12',
               'ylabel':'VLF34',
               'title':ptitle}


hod.plot_hodogram_dynamic(filt12z, filt34z, npts=4, gap=1, plot_kwargs=plot_kwargs,pauseT=0.03)

















print("here")






