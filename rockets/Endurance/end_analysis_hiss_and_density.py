"""
Compare VLF hiss with density cavities.

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
#import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import plot_hodogram_dynamic as hod

#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

#vDC = EFL('V12D')
#wfDC, tdatDC = vDC.load_data_gainphase_corrected()
#fsDC = vDC.chnspecs['fs']

v12 = EFL('VLF12D')
wf, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']

v34 = EFL('VLF34D')
wf34, tdat34 = v34.load_data_gainphase_corrected()






#v1s = EFL('V1SD')
#v2s = EFL('V2SD')
#wf1s, tdats = v1s.load_data_gainphase_corrected()
#wf2s, tdats = v2s.load_data_gainphase_corrected()



#----------------------------------------------------------------------
#Test for LHSS
#----------------------------------------------------------------------

fr = [4000,12000]
filt12 = filt.butter_bandpass_filter(wf,fr[0],fr[1],fs,order=2)
filt34 = filt.butter_bandpass_filter(wf34, fr[0],fr[1],fs,order=2)



t0 = 300
t1 = 400
good = np.where((tdat > t0) & (tdat < t1))[0]

wf12z = filt12[good]
wf34z = filt34[good]

fig,axs=plt.subplots(2)
axs[0].plot(tdat[good],wf12z)
axs[1].plot(tdat[good],wf34z)
for i in range(2):
    axs[i].set_xlim(330.923,331.026)
    axs[i].set_ylim(-1,1)


#hod.plot_hodogram_dynamic(wf12z, wf34z, npts=3, gap=1, pauseT=0.01)
t0z = 329.7470
t1z = 329.7655
goodz = np.where((tdat[good] > t0z) & (tdat[good] < t1z))[0]


i=0
di = 50
wf12zz = wf12z[goodz[i]:goodz[i+di]]
wf34zz = wf34z[goodz[i]:goodz[i+di]]
plt.plot(wf12zz,wf34zz)


#----------------------------------------------------------------------
#Mission timeline data
#----------------------------------------------------------------------

tl, gsS, gsE, bsS, bsE = end_data_loader.load_timeline()

fspec, tspec, powerc = signal.spectrogram(wf, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,10000],xr=[100,850], yscale='log')


#-------------------------------------------------------------------------------------
#SLP fixed bias data (end of mission only)
#-------------------------------------------------------------------------------------

slpfb = end_data_loader.load_slp_fixedbias()
timesdens = slpfb['ToF [s]']
dens = slpfb['Normalized Density [\m3]']


#5 sec cadence SLP data
slp5 = end_data_loader.load_slp()

#5 sec cadence SLP potential
slpp = end_data_loader.load_slp_potential()





#-------------------------------------------------------------------------------------
#Density structures
#--can be somewhat difficult to identify in skin data
#-------------------------------------------------------------------------------------


#wDC = filt.butter_bandpass_filter(wfDC, 1, 80, fsDC, order= 8)
#wDC = filt.butter_bandpass_filter(wfDC, 10,80, fsDC, order= 8)
w = filt.butter_bandpass_filter(wf, 5000, 12000, fs, order= 8)
#w1s = filt.butter_bandpass_filter(wf1s, 5, 80, fs, order= 8)
#w2s = filt.butter_bandpass_filter(wf2s, 5, 80, fs, order= 8)

fspec, tspec, powerc = signal.spectrogram(wf, fs, nperseg=256,noverlap=128,window='hann',return_onesided=True,mode='complex')
#fspecDC, tspecDC, powercDC = signal.spectrogram(wfDC, fsDC, nperseg=256,noverlap=128,window='hann',return_onesided=True,mode='complex')



tr = [200,300]
figs,axss = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,12000],xr=tr, yscale='linear',ax=axss[0],colorbar=0)
axss[1].plot(slp5['ToF [s]'], slp5['SLP Ni [/m3]'])
axss[1].set_xlim(tr)



#Test timetags for SLP data 
#Endurance_SLP_05_07_2023_witherrorbars.csv
#Endurance_GPS_velocity_position_altitude.csv

ephem, ephem2 = end_data_loader.load_ephemeris()


plt.plot(slp5['ToF [s]'],slp5['Alt [km]'],'.',ephem['Time'],ephem['Altitude'],'.')
#plt.xlim(100,110)
plt.plot(ephem['Time'],ephem['Altitude'])











figs,axss = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,12000],xr=[860,880], yscale='linear',ax=axss[0],colorbar=0)
#ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[1,100],xr=[700,900], yscale='linear',ax=axss[1])
axss[1].plot(timesdens, dens,'.')
axss[1].set_xlim(830,880)
axss[1].set_xlim(850,860)

#tr = [705,775]
tr = [110,115]
fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=tr, yscale='linear',ax=axs[0])
#ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[0,100],xr=tr, yscale='linear',ax=axs[1])
axs[3].plot(tdat,wf)
#axs[2].plot(tdatDC,wfDC)
axs[3].set_ylim(-1.3,1.3)
axs[2].set_ylim(-15.5,15.5)
for i in range(2): axs[i+2].set_xlim(tr)


tr = [113,115]
fig,axs = plt.subplots(3)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[4000,8000],xr=tr, yscale='linear',ax=axs[0])
axs[1].plot(tdat,wf)
axs[2].plot(tdat,w)
axs[1].set_ylim(-1.3,1.3)
axs[1].set_xlim(tr)
axs[2].set_ylim(-0.1,0.1)
axs[2].set_xlim(tr)




plt.plot(tdat,w, tdatDC,wDC)
plt.xlim(705,775)
plt.ylim(-0.25,0.25)





