"""
Determine polarization of Bernstein waves from KiNETX sounding rocket

UPLEG: the magnetic field (Bo) is in the direction of V1 x V3 and rocket (main) velocity is 
in the V1 x V4 direction.

DOWNLEG: Bo in the same direction as rocket velocity (V1 x V4). 

"""

import sys 
#sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
#import pyspedas
#from pytplot import tplot


#import plasma_params_get_density_from_flhr_freq as dflh
import numpy as np 
import matplotlib.pyplot as plt
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency
import plot_spectrogram as ps
import plot_hodogram_dynamic as hod
from scipy import signal
import correlation_analysis as ca
import fft_spectrum_piecewise as fsp


#from os.path import dirname, join as pjoin
#import scipy.io as sio
from scipy.io import readsav



#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]


#------------------------------------------------------------------------------
#VLF data of the Barium releases
#------------------------------------------------------------------------------

path = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/efield_VLF/'
fn = 'VLF12D_mvm_onVLF34ATimeTags_40kHz_Flight52007.sav'
sav_data12 = readsav(path+fn)

#fn = 'VLF34A_counts_TMwords13579_40kHz_Flight52007.sav'
fn = 'VLF34A_counts_mvm_TMwords13579_40kHz_Flight52007.sav'
sav_data34 = readsav(path+fn)


wf12 = sav_data12['data_vlf12d']
wf34 = sav_data34['data_vlf34a']
times12 = sav_data12['time_vlf12d']
times34 = sav_data34['time_vlf34a']




#-------------------------------------------
#polarization data from IDL Chaston routine
#-------------------------------------------

degpol = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/Bernstein_degpol.sav')
elip = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/Bernstein_elip.sav')
hlip = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/Bernstein_hlip.sav')
pow = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/Bernstein_pow.sav')


times_chaston = degpol['pol2']['x'][0]
freqs_chaston = degpol['pol2']['v'][0]

degpolspec = degpol['pol2']['y'][0]
elipspec = elip['elip2']['y'][0]
hlipspec = hlip['hlip2']['y'][0]
powspec = pow['pow']['y'][0]



#------------------------------------------------------------------------------
#DC-coupled data of the Barium releases
#------------------------------------------------------------------------------

"""
path = '/Users/abrenema/Desktop/Research/Rocket_missions/KiNETX/data/'
fn = '52007_TM2_S1-DCEfield_V12D_V34D_V56D_V15D_V62D_Volts_mvm.sav'
sav_data = readsav(path+fn)

times = sav_data['times']
wf12 = sav_data['v12d_mvm']
wf34 = sav_data['v34d_mvm']
wfz = sav_data['v15d_mvm']
"""
#------------------------------------------------------------------------------




#Remove values at negative times
wf12 = wf12[times12 >= 0]
wf34 = wf34[times34 >= 0]
times12 = times12[times12 >= 0]
times34 = times34[times34 >= 0]

times = times12
#sr = sav_data12['samplerate']

tdelta = times - np.roll(times,1)

fs = 40000.

#Spectral overview
nps = 256
freq12, tspec12, power12 = signal.spectrogram(wf12, fs, nperseg=nps,noverlap=nps/2,window='hann', return_onesided=False)
freq34, tspec34, power34 = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann', return_onesided=False)


#turn the two-sided power array into a complex array
indextmp = int(nps/2)






#Isolate waveform for polarization

#tr = [629.5, 631]
#fr = [3250, 3450]
tr = [630, 631]
fr = [3400, 3600]



#freqs12_fin, tcenter12_fin, spec12_fin, fs12_fin = fsp.fft_spectrum_piecewise(times, wf12, nfft=512, noverlap=2, fs_thres=0.4)
#freqs34_fin, tcenter34_fin, spec34_fin, fs34_fin = fsp.fft_spectrum_piecewise(times, wf34, nfft=512, noverlap=2, fs_thres=0.4)


#Reduce waveforms to only times needed
t0 = 600
t1 = 700
goodt = np.where((times > t0 ) & (times < t1))
wf12z = wf12[goodt]
wf34z = wf34[goodt]
timesz = times[goodt]



freqs_fin, tcenter_fin, csd_fin,coh_fin,phase_fin,spec_fin1,spec_fin2, fs_fin = ca.csd_spectrum_piecewise(times,wf12,wf34,fs_thres=2)
phase_fin = phase_fin * (180/3.14)

#xr = [627.5,632.5]
xr = [620,640]
#xr = [580,610]
yr = [1000,10000]
fig,axs = plt.subplots(3)
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],np.abs(spec_fin1),vr=[-80,-20],yr=yr,xr=xr, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],coh_fin,vr=[0,1],yr=yr,xr=xr,yscale='linear',zscale='linear',ax=axs[1])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],phase_fin,vr=[-180,180],yr=yr,xr=xr,yscale='linear',zscale='linear',ax=axs[2])






#Limit the plots by the degree of polarization
degpolspec2 = np.copy(degpolspec)
elipspec2 = np.copy(elipspec)
hlipspec2 = np.copy(hlipspec)
powspec2 = np.copy(powspec)

#polarization limit
pollim = 0.8
for i in range(len(times_chaston)):
    badd = np.where(degpolspec[:,i] < pollim)
    degpolspec2[badd,i] = np.nan
    elipspec2[badd,i] = np.nan
    hlipspec2[badd,i] = np.nan
    powspec2[badd,i] = np.nan


#Limit the plots by wave power 
degpolspec3 = np.copy(degpolspec)
elipspec3 = np.copy(elipspec)
hlipspec3 = np.copy(hlipspec)
powspec3 = np.copy(powspec)

#polarization limit
#powlim = 1e-11  #first barium release (isolates Bernstein waves)
powlim = 5e-10  #second barium release (isolates Bernstein waves)
for i in range(len(times_chaston)):
    badd = np.where(powspec[:,i] < powlim)
    degpolspec3[badd,i] = np.nan
    elipspec3[badd,i] = np.nan
    hlipspec3[badd,i] = np.nan
    powspec3[badd,i] = np.nan





#tslice = 597.19  #4600 Hz peak
#tslice = 598.653  #4600 Hz peak
#tslice = 630.6  #4600 Hz peak
tslice = 630.66  #4600 Hz peak

#xr = [592,610]
#xr = [593,600]
xr = [627,633]

#-----------------------------------
#Plot values filtered by deg pol
#-----------------------------------


fig,axs = plt.subplots(8, figsize=(14,9))
axs[0].set_title('Bernstein from Barium release (kinetx_analysis_polarization.py)\ndegPol_min='+str(pollim))
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],np.abs(spec_fin1),vr=[-60,-10],yr=[2000,8000],xr=xr, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],np.abs(spec_fin2),vr=[-50,-15],yr=[2000,8000],xr=xr, yscale='linear',ax=axs[1])
ps.plot_spectrogram(times_chaston,freqs_chaston,np.abs(powspec3),xr=xr,yr=[2000,8000],vr=[-140,-60],yscale='linear',zscale='log',ax=axs[2])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],phase_fin,vr=[-180,180],yr=[2000,8000],xr=xr,yscale='linear',zscale='linear',ax=axs[3])
ps.plot_spectrogram(times_chaston,freqs_chaston,degpolspec2,xr=xr,yr=[2000,8000],vr=[pollim,1],yscale='linear',zscale='linear',ax=axs[4])
ps.plot_spectrogram(times_chaston,freqs_chaston,elipspec2,xr=xr,yr=[2000,8000],vr=[-1,1],yscale='linear',zscale='linear',ax=axs[5])
ps.plot_spectrogram(times_chaston,freqs_chaston,hlipspec2,xr=xr,yr=[2000,8000],vr=[0,1],yscale='linear',zscale='linear',ax=axs[6])
powavg, powarr, tarr = ps.slice_spectrogram(tslice, tcenter_fin, np.abs(spec_fin1), nsec=0.2)
axs[7].plot(freqs_fin[256:,0], powavg[256:])
axs[7].set_xscale('linear')
axs[7].set_yscale('log')
axs[7].set_xlim(500,10000)
axs[7].set_ylim(1e-6,1e-1)
axs[7].set_ylabel('power slice at\n'+ str(tslice))
for i in range(6):
    axs[i].vlines(tslice,0,10000)
axs[0].set_ylabel('sonogram\nVLF12')
axs[1].set_ylabel('sonogram\nVLF34')
axs[3].set_ylabel('phase')
axs[4].set_ylabel('deg pol')
axs[5].set_ylabel('ellipticity')
axs[6].set_ylabel('helicity')
axs[7].set_xlabel('freq (Hz)')





#-----------------------------------
#Plot values filtered by wave power
#-----------------------------------

fig,axs = plt.subplots(7, figsize=(14,9))
axs[0].set_title('Bernstein from Barium release (kinetx_analysis_polarization.py)\nPower_min='+str(powlim))
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],np.abs(spec_fin1),vr=[-60,-10],yr=[2000,8000],xr=xr, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],np.abs(spec_fin2),vr=[-50,-15],yr=[2000,8000],xr=xr, yscale='linear',ax=axs[1])
ps.plot_spectrogram(times_chaston,freqs_chaston,np.abs(powspec3),xr=xr,yr=[2000,8000],vr=[-90,-70],yscale='linear',zscale='log',ax=axs[2])
ps.plot_spectrogram(tcenter_fin,freqs_fin[:,0],phase_fin,vr=[-180,180],yr=[2000,8000],xr=xr,yscale='linear',zscale='linear',ax=axs[3])
ps.plot_spectrogram(times_chaston,freqs_chaston,degpolspec3,xr=xr,yr=[2000,8000],vr=[pollim,1],yscale='linear',zscale='linear',ax=axs[4])
ps.plot_spectrogram(times_chaston,freqs_chaston,elipspec3,xr=xr,yr=[2000,8000],vr=[-1,1],yscale='linear',zscale='linear',ax=axs[5])
ps.plot_spectrogram(times_chaston,freqs_chaston,hlipspec3,xr=xr,yr=[2000,8000],vr=[0,1],yscale='linear',zscale='linear',ax=axs[6])
#powavg, powarr, tarr = ps.slice_spectrogram(tslice, tcenter_fin, np.abs(spec_fin1), nsec=0.2)
#axs[7].plot(freqs_fin[256:,0], powavg[256:])
#axs[7].set_xscale('linear')
#axs[7].set_yscale('log')
#axs[7].set_xlim(500,10000)
#axs[7].set_ylim(1e-6,1e-1)
#axs[7].set_ylabel('power slice at\n'+ str(tslice))
#for i in range(6):
#    axs[i].vlines(tslice,0,10000)
axs[0].set_ylabel('sonogram\nVLF12')
axs[1].set_ylabel('sonogram\nVLF34')
axs[2].set_ylabel('sonogram\nChaston')
axs[3].set_ylabel('phase')
axs[4].set_ylabel('deg pol')
axs[5].set_ylabel('ellipticity')
axs[6].set_ylabel('helicity')
#axs[7].set_xlabel('freq (Hz)')



















#Bandpass filter waveform to preferred freqs and times
filt12 = filter_wave_frequency.butter_bandpass_filter(wf12, fr[0],fr[1],fs,order=2)
filt34 = filter_wave_frequency.butter_bandpass_filter(wf34, fr[0],fr[1],fs,order=2)






#all the points to be plotted
#trz= [629.60, 629.605]
trz= [630.0, 631.0]

goodt = np.squeeze(np.where((times > trz[0]) & (times < trz[1])))
filt12z = filt12[goodt]
filt34z = filt34[goodt]

ptitle = 't=' + str(trz[0]) + '-' + str(trz[1]) + ' sec\nf='+str(fr[0]) + '-' + str(fr[1]) + 'Hz'
num = 1


plt.plot(filt12z,filt34z)
plt.ylim(-5,5)
plt.xlim(-5,5)
#plt.plot(np.roll(filt12z,-1),filt34z)

plot_kwargs = {'xlim':[-5,5],
               'xlabel':'VLF12',
               'ylabel':'VLF34',
               'title':ptitle}


hod.plot_hodogram_dynamic(filt12z, filt34z, npts=2, gap=2, plot_kwargs=plot_kwargs,pauseT=0.02)

















print("here")






