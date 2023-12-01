"""
Analyze artificial (probably) signals on Endurance

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
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import interferometry_routines as interf





#--------------------------------------------------------------
#Get timeline of data to separate out science collection times 
#--------------------------------------------------------------

tl, gsS, gsE = end_data_loader.load_timeline()


#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

"""
v12DC = EFL('V12D')
v34DC = EFL('V34D')
v32DC = EFL('V32D')
v24DC = EFL('V24D')
wf12DC, tdatDC = v12DC.load_data_gainphase_corrected()
wf34DC, tdatDC = v34DC.load_data_gainphase_corrected()
wf32DC, tdatDC = v32DC.load_data_gainphase_corrected()
wf24DC, tdatDC = v24DC.load_data_gainphase_corrected()
fsDC = v12DC.chnspecs['fs']
"""



v12 = EFL('VLF12D')
v13 = EFL('VLF13D')
v41 = EFL('VLF41D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')

fs = v12.chnspecs['fs']
wf12, tdat = v12.load_data_gainphase_corrected()
wf13, tdat = v13.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf41, tgoo = v41.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
wf42 = -wf24
wf14 = -wf41


"""
v1 = EFL('V1SD')
v2 = EFL('V2SD')
v3 = EFL('V3SD')
v4 = EFL('V4SD')

wf1, tdats = v1.load_data_gainphase_corrected()
wf2, tdats = v2.load_data_gainphase_corrected()
wf3, tdats = v3.load_data_gainphase_corrected()
wf4, tdats = v4.load_data_gainphase_corrected()
fss = v1.chnspecs['fs']
"""




#---------------------------------------------------------------------------
#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
#nps = 1024
#fspecx, tspecx, powercAx = signal.spectrogram(wfAx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
#fspecx, tspecx, powercBx = signal.spectrogram(wfBx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
#fspecy, tspecy, powercAy = signal.spectrogram(wfAy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
#fspecy, tspecy, powercBy = signal.spectrogram(wfBy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')

#----------------------------------------------------------------------
#Get spectral data for finding waves
#----------------------------------------------------------------------

fspec, tspec, powerc12 = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc13 = signal.spectrogram(wf13, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc34 = signal.spectrogram(wf34, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc42 = signal.spectrogram(wf42, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc14 = signal.spectrogram(wf14, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc32 = signal.spectrogram(wf32, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')

#fspecDC, tspecDC, powerc12DC = signal.spectrogram(wf12DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecDC, tspecDC, powerc34DC = signal.spectrogram(wf34DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecDC, tspecDC, powerc32DC = signal.spectrogram(wf32DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecDC, tspecDC, powerc24DC = signal.spectrogram(wf24DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')

#fspecs, tspecs, powerc1s = signal.spectrogram(wf1, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecs, tspecs, powerc2s = signal.spectrogram(wf2, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecs, tspecs, powerc3s = signal.spectrogram(wf3, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
#fspecs, tspecs, powerc4s = signal.spectrogram(wf4, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')



#----------------------------------------------------------------------------------------
#Calculate the fractional difference b/t different spectra to identify artificial waves
#----------------------------------------------------------------------------------------

ptmp_diff = np.abs(np.abs(powerc12) - np.abs(powerc34))
ptmp_sum = np.abs(powerc12) + np.abs(powerc34)
ptmp_fracdiff = ptmp_diff/ptmp_sum

Nval = 6
gtmp,cohtmp,phasetmp = correlation_analysis.interferometric_coherence_2D(powerc12,powerc34,Nval)


yr = [4000,9000]
#xr = [0,900]
xr = [700,900]
vr = [-40,-25]
ys = 'linear'
fig, axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),ax=axs[0],vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power12\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),ax=axs[1],vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power34\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,np.abs(ptmp_fracdiff),ax=axs[2],vr=[0,1],yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='frac diff\nf(Hz)',zscale='linear')
ps.plot_spectrogram(tspec,fspec,cohtmp,ax=axs[3],vr=[0,1],yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='Coherence\nf(Hz)',zscale='linear')











cohmin = 0.5  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result




#NOTE: + sense of phase defined as pointing towards center of potential of "powercA"
#Nval = 3
#Nval = 10
#gx,cohx,phasex = correlation_analysis.interferometric_coherence_2D(powercAx,powercBx,Nval)
#gy,cohy,phasey = correlation_analysis.interferometric_coherence_2D(powercAy,powercBy,Nval)



goo = cohx < cohmin
cohx[goo] = float("nan")
phasex[goo] = float("nan")
phasex = np.degrees(phasex)

goo = cohy < cohmin
cohy[goo] = float("nan")
phasey[goo] = float("nan")
phasey = np.degrees(phasey)



#Reduce arrays to desired timerange
goo, powercAxz, tspecxz = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)
goo, powercBxz, tpowercBxz = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)
goo, cohxz, ttmp = ps.slice_spectrogram(tz,tspecx,cohx,nsec)
goo, phasexz, ttmp = ps.slice_spectrogram(tz,tspecx,phasex,nsec)

goo, powercAyz, tspecyz = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)
goo, powercByz, tpowercByz = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)
goo, cohyz, ttmp = ps.slice_spectrogram(tz,tspecy,cohy,nsec)
goo, phaseyz, ttmp = ps.slice_spectrogram(tz,tspecy,phasey,nsec)


#Endurance has ~3m booms, so the 
receiver_spacing = 2.1 #meters -- effective length of spaced receiver (=3*cos(45))






#Observation 1 - 5700 Hz signal almost entirely in x-hat' direction. 
#Test to see if E-field is only (mostly) contained in this component. 

v12 = EFL('VLF12D')
v34 = EFL('VLF34D')
wf12, tdatx = v12.load_data_gainphase_corrected()
wf34, tdatx = v34.load_data_gainphase_corrected()





#yr = [0,12000]
yr = [4000,14000]
#xr = [100,260]
#xr = [700,900]
xr = [700,775]
#vr = [-45,20]
vr = [-40,-25]
ftmp, ttmp, ptmp = signal.spectrogram(wfAx, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')
ftmp, ttmp, ptmp = signal.spectrogram(wfBx, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')

ftmp, ttmp, ptmp = signal.spectrogram(wfAy, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')
ftmp, ttmp, ptmp = signal.spectrogram(wfBy, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')

ftmp, ttmp, ptmp = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')
ftmp, ttmp, ptmp = signal.spectrogram(wf34, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')























#--------------------------------------------
#Oddity in V34 and V32 not seen in V12 or V24
#--------------------------------------------

fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc24),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc32),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[3])

fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powerc12DC),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powerc34DC),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powerc24DC),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powerc32DC),vr=[-40,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[3])

fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspecs,fspecs,np.abs(powerc1s),vr=[-60,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecs,fspecs,np.abs(powerc2s),vr=[-60,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspecs,fspecs,np.abs(powerc3s),vr=[-60,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspecs,fspecs,np.abs(powerc4s),vr=[-60,-25],yr=[0,2000],xr=[110,220], yscale='linear',ax=axs[3])


#tslice = 180
#nsec = 5
tslice = 450
nsec = 50
p12DCavg, p12DCvals = ps.slice_spectrogram(tslice,tspecDC,np.abs(powerc12DC),nsec)
p34DCavg, p34DCvals = ps.slice_spectrogram(tslice,tspecDC,np.abs(powerc34DC),nsec)
p24DCavg, p24DCvals = ps.slice_spectrogram(tslice,tspecDC,np.abs(powerc24DC),nsec)
p32DCavg, p32DCvals = ps.slice_spectrogram(tslice,tspecDC,np.abs(powerc32DC),nsec)
p12avg, p12vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc12),nsec)
p34avg, p34vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc34),nsec)
p24avg, p24vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc24),nsec)
p32avg, p32vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc32),nsec)


fig,axs = plt.subplots(4)
axs[0].plot(fspecDC,p12DCavg,color='black')
axs[0].plot(fspec,p12avg,color='blue')
axs[1].plot(fspecDC,p34DCavg,color='black')
axs[1].plot(fspec,p34avg,color='blue')
axs[2].plot(fspecDC,p24DCavg,color='black')
axs[2].plot(fspec,p24avg,color='blue')
axs[3].plot(fspecDC,p32DCavg,color='black')
axs[3].plot(fspec,p32avg,color='blue')
for i in range(4): axs[i].set_xlim(0,15000)
for i in range(4): axs[i].set_ylim(0,0.0005)




fig,axs = plt.subplots(1)
axs.plot(fspec,np.sqrt(p12avg),color='black')
axs.plot(fspec,np.sqrt(p34avg),color='blue')
axs.plot(fspec,np.sqrt(p24avg),color='red')
axs.plot(fspec,np.sqrt(p32avg),color='orange')
axs.set_xlim(0,5000)
#axs.set_xlim(0,200)
axs.set_ylim(0,0.04)
































"""
#-----------------------------------
#Artifical signal with harmonics that slowly rise in tone over entire mission
#-----------------------------------

ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[80,200],xr=[380,420], yscale='linear')

tr = [388,395]
good = np.where((tdats >= tr[0]) & (tdats <= tr[1]))


w1 = filt.butter_highpass_filter(wf1[good], 1, fss, order= 10)
w2 = filt.butter_highpass_filter(wf2[good], 1, fss, order= 10)
w3 = filt.butter_highpass_filter(wf3[good], 1, fss, order= 10)

w4 = filt.butter_highpass_filter(wf4[good], 1, fss, order= 10)



plt.plot(tdats[good],w1,tdats[good],w2,tdats[good],w3,tdats[good],w4)
plt.xlim(390.4,390.5)

fig,axs = plt.subplots(2)
axs[0].plot(tdats[good],w1)
axs[1].plot(tdats[good],w2)


plot_csd_skins([390,392])

"""
