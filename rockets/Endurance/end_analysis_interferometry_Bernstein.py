"""
Run interferometry analysis on Bernstein waves from Endurance
-2D power spec of f vs k
-1D k vs f 

(based on interferometry_routines_call.py)

NOTE: this technique only makes sense for parallel, spaced antennas. 
If the antennas are perpendicular, the ~90 deg phase shift in waves will be interpreted as a 
large k-value. 

NOTE: phase identification can be difficult if:

1) the wave is propagating nearly along one of the interferometry axes (e.g. if along the y-hat' direction then V32 and V14 will
measure very little signal, meaning the waveform can be noise-dominated)

2) if wavelength is on the order or less than boom length than results can't be trusted.


                       (y-hat)
                         V3
                      /  |  \             
           (y-hat') z    |     z (x-hat')         
                  /      |       \        
               V2--------O--------V1 (x-hat)    
                  \      |       /
                    z    |     z
                      \  |  /
                         V4
          

z-points represent the centers of potential of the interferometry diagonals          

Coordinate system (system of input test wave)
x-hat --> E12 = V1 - V2 direction (positive to right)
y-hat --> E34 = V3 - V4 direction (positive upwards)

This code outputs the spectrum of k-values in the kx' and ky' directions, where
x'-hat --> Ex' = V1V3z - V4V2z (45 deg inclined from x-hat)
y'-hat --> Ey' = V3V2z - V1V4z (45 deg inclined from y-hat)


#---interferometry pairs
#NOTE: positive sense of phase defined as pointing towards center of potential of "wfA"
#e.g. y-hat':  wfA = wf32; wfB = wf14
#So, for Endurance the interferometry pairs are 
#   vAstr = 'VLF32D' and vBstr = 'VLF41D'-->'VLF14D'
#or vAstr = 'VLF13D' and vBstr = 'VLF24D'-->'VLF42D'


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


#-------------------------------------------------------------
#Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using Endurance data
#-------------------------------------------------------------


#----------------------------------------------------------------
#Use long spaced receivers
#----------------------------------------------------------------
##x-hat' direction from phase analysis uses these two components
#vAstrx = 'VLF13D'
#vBstrx = 'VLF24D'  #-->42
#polarity_xA = 1
#polarity_xB = -1
##y-hat' direction from phase analysis uses these two components
#vAstry = 'VLF32D'
#vBstry = 'VLF41D'  #-->14
#polarity_yA = 1
#polarity_yB = -1
##Endurance has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
#receiver_spacing_xhat = 2.27 #m 
#receiver_spacing_yhat = 2.27 #m

#----------------------------------------------------------------
#Use one short spaced receiver and one long
#----------------------------------------------------------------
##x-hat' direction from phase analysis uses these two components
vAstrx = 'VLF13D'
vBstrx = 'VLF34D'  #-->43
polarity_xA = 1
polarity_xB = -1
##y-hat' direction from phase analysis uses these two components
vAstry = 'VLF32D'
vBstry = 'VLF41D' #-->14
polarity_yA = 1
polarity_yB = -1
##Endurance has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
receiver_spacing_xhat = 2.27/2 #m 
receiver_spacing_yhat = 2.27 #m

#----------------------------------------------------------------
#Avoid using V1 which seems to have spurious noise near 6 kHz
#----------------------------------------------------------------
##x-hat' direction from phase analysis uses these two components
#vAstrx = 'VLF34D'
#vBstrx = 'VLF24D'
#polarity_xA = 1
#polarity_xB = 1
##y-hat' direction from phase analysis uses these two components
#vAstry = 'VLF32D'
#vBstry = 'VLF34D'
#polarity_yA = 1
#polarity_yB = 1
##Endurance has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
#receiver_spacing_xhat = 2.27/2 #m 
#receiver_spacing_yhat = 2.27/2 #m





#Load Endurance waveform
vAx = EFL(vAstrx)
vBx = EFL(vBstrx)
wfAx, tdatx = vAx.load_data_gainphase_corrected()
wfBx, tdatx = vBx.load_data_gainphase_corrected()
vAy = EFL(vAstry)
vBy = EFL(vBstry)
wfAy, tdaty = vAy.load_data_gainphase_corrected()
wfBy, tdaty = vBy.load_data_gainphase_corrected()


#change polarity as needed
wfAx *= polarity_xA
wfBx *= polarity_xB
wfAy *= polarity_yA
wfBy *= polarity_yB


#for plot titles
xptitle = str(polarity_xA) + '*' + vAstrx + ' and\n' + str(polarity_xB) + '*' + vBstrx
yptitle = str(polarity_yA) + '*' + vAstry + ' and\n' + str(polarity_yB) + '*' + vBstry



#sample rate
fs = vAx.chnspecs['fs']



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 1024
fspecx, tspecx, powercAx = signal.spectrogram(wfAx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecx, tspecx, powercBx = signal.spectrogram(wfBx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecy, tspecy, powercAy = signal.spectrogram(wfAy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecy, tspecy, powercBy = signal.spectrogram(wfBy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')



cohmin = 0.3  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
vr = [-45,-20]
ys = 'linear'
kr = [-5,5]
#tz = 140 
#nsec = 5
#tz = 127
#nsec = 5
#tz = 115
#nsec = 2
#tz = 170 
#tz = 721.5
#nsec = 5
tz = 720
nsec = 5
#tz = 885
#nsec = 5
#tz = 381.4
#nsec = 3
yr = [0,16000]




#-------------------------------------------------------------
#Method 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using Endurance data
#-------------------------------------------------------------

#NOTE: + sense of phase defined as pointing towards center of potential of "powercA"
Nval = 3
#Nval = 10
gx,cohx,phasex = correlation_analysis.interferometric_coherence_2D(powercAx,powercBx,Nval)
gy,cohy,phasey = correlation_analysis.interferometric_coherence_2D(powercAy,powercBy,Nval)


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




#Get power spec of k vs freq
fkpowspecx, kvalsx, fvalsx, pmaxvalsx = interf.inter_fvsk(np.abs(powercAxz),tspecxz,fspecx, 
                                         np.radians(phasexz),tspecxz,fspecx,
                                         receiver_spacing_xhat,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-5,5])
fkpowspecy, kvalsy, fvalsy, pmaxvalsy = interf.inter_fvsk(np.abs(powercAyz),tspecyz,fspecy, 
                                         np.radians(phaseyz),tspecyz,fspecy,
                                         receiver_spacing_yhat,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-5,5])

#Turn k-values into wavelength
wl1 = np.zeros(len(fvalsx))
for i in range(len(fvalsx)):
    tmp = pmaxvalsx[i,:]
    kxgoo = kvalsx[np.where(tmp == 1)]
    tmp = pmaxvalsy[i,:]
    kygoo = kvalsy[np.where(tmp == 1)]
    kmaggoo = np.sqrt(kxgoo**2 + kygoo*2)
    if len(kmaggoo != 0): 
        wl1[i] = 2*np.pi/kmaggoo



#-------------------------------------------------------------
#Method 2: Do a freq (xaxis) vs k-value (yaxis) interferometry analysis using Endurance data
#-------------------------------------------------------------


tchunk = 0.1  #delta-time (sec) for each time chunk to divide up the spectra into
nchunks = int(np.ceil((wfAx.size/fs)/tchunk)) #number of chunks in ENTIRE timerange

cohx2, phasex2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAx,wfBx,tdatx,fs,tchunk,coh_min=cohmin,nperseg=512)
cohy2, phasey2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAy,wfBy,tdaty,fs,tchunk,coh_min=cohmin,nperseg=512)


#Reduce arrays to desired timerange
pavgx2, phasearrx2, tarr_phasex2 = ps.slice_spectrogram(tz,tchunks2,phasex2,nsec)
cavgx2, coharrx2, tarr_cohx2 = ps.slice_spectrogram(tz,tchunks2,cohx2,nsec)
powavgx2, powarrx2, tarr_powx2 = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)

pavgy2, phasearry2, tarr_phasey2 = ps.slice_spectrogram(tz,tchunks2,phasey2,nsec)
cavgy2, coharry2, tarr_cohy2 = ps.slice_spectrogram(tz,tchunks2,cohy2,nsec)
powavgy2, powarry2, tarr_powy2 = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)


"""
#Test out turning negative phase differences close to 180 to smaller values 
#phasearrx2[phasearrx2 <= -90] += 180
#phasearry2[phasearry2 <= -90] += 180
pavgx2 = np.asarray(pavgx2)
pavgy2 = np.asarray(pavgy2)

#pavgx2 = np.degrees(np.unwrap(np.radians(pavgx2)))

pavgx2[pavgx2 <= 0] += 180
pavgy2[pavgy2 <= 0] += 180
"""




#Turn the phase values into k-values 
kx2 = [np.radians(i) / receiver_spacing_xhat for i in pavgx2]
ky2 = [np.radians(i) / receiver_spacing_yhat for i in pavgy2]
#and then into wavelength perp to Bo
kmag = [np.sqrt(kx2[i]**2 + ky2[i]**2) for i in range(len(kx2))]
wl2 = [2*np.pi/i for i in kmag]

yr = [5000,9000]
vr=[-60,-20]
xrspec = [tz,tz+nsec]

#----------------------------------------------------------------------------------------
#Calculate the fractional difference b/t different spectra to identify artificial waves
#----------------------------------------------------------------------------------------

ptmp_diff = np.abs(np.abs(powarrx2) - np.abs(powarry2))
ptmp_sum = np.abs(powarrx2) + np.abs(powarry2)
ptmp_fracdiff = ptmp_diff/ptmp_sum


#Plot values from Method 2
fig,axs = plt.subplots(6,2, figsize=(8,8))  #,gridspec_kw={'height_ratios':[1,1,1,1,1,1,1,1,1]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.4,wspace=0.4)
#fig,axs = plt.subplots(5,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1],'width_ratios':[1,0.5]})
ps.plot_spectrogram(tspecxz,fspecx,np.abs(powarrx2),vr=vr,xr=xrspec,yr=yr,yscale='linear',ax=axs[0,0],xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
ps.plot_spectrogram(tarr_cohx2,freqs2,coharrx2,vr=[0.5,1], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[1,0],xlabel='time(sec)',ylabel='coherence\nfreq(Hz)')
ps.plot_spectrogram(tarr_phasex2,freqs2,phasearrx2,vr=[-140,140], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[2,0],xlabel='time(sec)',ylabel='phase(deg)\nfreq(Hz)')
ps.plot_spectrogram(tspecxz,fspecx,ptmp_fracdiff,vr=[0,1],xr=xrspec,yr=yr,zscale='linear',yscale='linear',ax=axs[3,0],xlabel='time(sec)',ylabel='fractional diff\npowspec1,\npowspec2\nfreq(Hz)')
for i in range(3): axs[i,0].axvline(tz, linestyle='--')
for i in range(3): axs[i,0].axvline(tz + nsec, linestyle='--')
axs[0,1].plot(fspecx,powarrx2)
axs[0,1].plot(fspecx,powavgx2,'.',color='black')
axs[0,1].set_xlim(yr)
axs[0,1].set_ylim(0,np.nanmax(powarrx2))
axs[0,1].set_ylabel('powermax')
axs[0,1].set_xlabel('freq(Hz)')
axs[1,1].plot(freqs2,coharrx2)
axs[1,1].plot(freqs2,cavgx2,'.',color='black')
axs[1,1].set_xlim(yr)
axs[1,1].set_ylim(0,1)
axs[1,1].set_ylabel('coherence')
axs[1,1].set_xlabel('freq(Hz)')
axs[2,1].plot(freqs2,phasearrx2)
axs[2,1].plot(freqs2,pavgx2,'.',color='black')
axs[2,1].set_xlim(yr)
axs[2,1].set_ylim(-180,180)
axs[2,1].set_ylabel("phase (x')")
axs[2,1].set_xlabel('freq(Hz)')
axs[3,1].plot(freqs2,phasearry2)
axs[3,1].plot(freqs2,pavgy2,'.',color='black')
axs[3,1].set_xlim(yr)
axs[3,1].set_ylim(-180,180)
axs[3,1].set_ylabel("phase (y')")
axs[3,1].set_xlabel('freq(Hz)')
axs[4,1].plot(freqs2,kx2,'.',freqs2,ky2,'.')
axs[4,1].set_xlim(yr)
axs[4,1].set_ylim(min(np.nanmin(kx2),np.nanmin(ky2)),max([np.nanmax(kx2),np.nanmax(ky2)]))
axs[4,1].set_ylabel("kx'(blue)\nky'(orange)")
axs[4,1].set_xlabel('freq(Hz)')
axs[5,1].plot(freqs2,wl2,'.',color='black')
axs[5,1].set_xlim(yr)
axs[5,1].set_ylim(0,50)
#axs[5,1].set_ylim(0,np.nanmax(wl2[wl2 != np.inf]))
axs[5,1].set_ylabel("wavelength(m)\nfrom |k|")
axs[5,1].set_xlabel('freq(Hz)')

fig.delaxes(axs[4,0])
fig.delaxes(axs[5,0])






#------------------------------------------------------------
#Plot final results
#------------------------------------------------------------

titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx
xr = [tz-3*nsec,tz+3*nsec]

krplot = [-2,2]

fig,axs = plt.subplots(5,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1,1],'width_ratios':[1,0.5]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.6,wspace=0.4)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercAx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[0,0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercBx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[1,0],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,cohx,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[2,0],xlabel='time(s)',ylabel='Coherence\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,np.abs(phasex),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[3,0],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
ps.plot_spectrogram(tspecxz,fspecx,ptmp_fracdiff,vr=[0,1],xr=xr,yr=yr,zscale='linear',yscale='linear',ax=axs[4,0],xlabel='time(sec)',ylabel='fractional diff\npowspec1,\npowspec2\nfreq(Hz)')

for i in range(5): axs[i,0].axvline(tz, linestyle='--')
for i in range(5): axs[i,0].axvline(tz+nsec, linestyle='--')

ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[0,1],minzval=0,maxzval=1,xlabel="kx'(rad/m)",ylabel=xptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[0,1],minzval=-120,maxzval=10,alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1,1],minzval=0,maxzval=1,xlabel="ky'(rad/m)",ylabel=yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1,1],minzval=-120,maxzval=10,alpha=0.5)
#Plot the limiting k-vector value where short wavelength effects start to occur. 
#i.e. location when wavelength equals about twice the interferometry receiver spacing
klimx = 2*np.pi / (2 * receiver_spacing_xhat)
klimy = 2*np.pi / (2 * receiver_spacing_yhat)
wlimx = (2 * receiver_spacing_xhat)
wlimy = (2 * receiver_spacing_yhat)
axs[0,1].axvline(klimx,linestyle="--")
axs[0,1].axvline(-klimx,linestyle="--")
axs[1,1].axvline(klimy,linestyle="--")
axs[1,1].axvline(-klimy,linestyle="--")

#Oplot the results from the initial (1d) analysis (Method 2)
axs[0,1].plot(kx2,freqs2,'*',color='black')
axs[1,1].plot(ky2,freqs2,'*',color='black')

axs[2,1].plot(wl1,fspecx,'.',color='black',markersize=2)
axs[2,1].plot(wl2,freqs2,'*',color='black')
axs[2,1].set_ylim(yr)
axs[2,1].set_xscale('linear')
axs[2,1].set_xlim(0,50)
axs[2,1].set_xlabel("wavelength(m)\n(from |k|)")
axs[2,1].set_ylabel("f(Hz)")
axs[2,1].axvline(wlimx,linestyle="--")
axs[2,1].axvline(wlimy,linestyle="--")

fig.delaxes(axs[3,1])
fig.delaxes(axs[4,1])










#-----------------------------
#Plot focused on the k spectra
#-----------------------------

fig,axs = plt.subplots(3)
fig.subplots_adjust(bottom=0.2,right=0.8,left=0.2,top=0.9,hspace=0.6,wspace=0.4)
ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[0],minzval=0,maxzval=1,xlabel="kx'(rad/m)",ylabel=xptitle+'\nf(Hz)',cmap='Greys',title=titlegoo)
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[0],minzval=-120,maxzval=10,alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1],minzval=0,maxzval=1,xlabel="ky'(rad/m)",ylabel=yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1],minzval=-120,maxzval=10,alpha=0.5)


axs[0].axvline(klimx, linestyle='--')
axs[0].axvline(-klimx, linestyle='--')
axs[1].axvline(klimy, linestyle='--')
axs[1].axvline(-klimy, linestyle='--')

#Oplot the results from the initial (1d) analysis
axs[0].plot(kx2,freqs2,'*',color='black')
axs[1].plot(ky2,freqs2,'*',color='black')

axs[2].plot(wl1,fspecx,'.',color='black',markersize=2)
axs[2].plot(wl2,freqs2,'*',color='black')
axs[2].set_ylim(yr)
axs[2].set_xscale('linear')
axs[2].set_xlim(0,50)
axs[2].set_xlabel("wavelength(m)\n(from |k|)")
axs[2].set_ylabel("f(Hz)")
axs[2].axvline(wlimx,linestyle="--")
axs[2].axvline(wlimy,linestyle="--")





print('h')
print('h')
print('h')
print('h')









#Observation 1 - 5700 Hz signal almost entirely in x-hat' direction. 
#Test to see if E-field is only (mostly) contained in this component. 

v12 = EFL('VLF12D')
v34 = EFL('VLF34D')
wf12, tdatx = v12.load_data_gainphase_corrected()
wf34, tdatx = v34.load_data_gainphase_corrected()




#get rid of DC offset
wfAxhp = filt.butter_highpass_filter(wfAx, 20, fs, order=5)
wfBxhp = filt.butter_highpass_filter(wfBx, 20, fs, order=5)
wfAyhp = filt.butter_highpass_filter(wfAy, 20, fs, order=5)
wfByhp = filt.butter_highpass_filter(wfBy, 20, fs, order=5)
wf12hp = filt.butter_highpass_filter(wf12, 20, fs, order=5)
wf34hp = filt.butter_highpass_filter(wf34, 20, fs, order=5)



wfAxbp = filt.butter_bandpass_filter(wfAxhp, 5000, 5900, fs, order=9)
wfBxbp = filt.butter_bandpass_filter(wfBxhp, 5000, 5900, fs, order=9)
wfAybp = filt.butter_bandpass_filter(wfAyhp, 5000, 5900, fs, order=9)
wfBybp = filt.butter_bandpass_filter(wfByhp, 5000, 5900, fs, order=9)
wf12bp = filt.butter_bandpass_filter(wf12hp, 5000, 5900, fs, order=9)
wf34bp = filt.butter_bandpass_filter(wf34hp, 5000, 5900, fs, order=9)

#...test filtering
ftmp, ttmp, ptmp = signal.spectrogram(wfAxbp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')
ftmp, ttmp, ptmp = signal.spectrogram(wfBxbp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')

ftmp, ttmp, ptmp = signal.spectrogram(wfAybp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')
ftmp, ttmp, ptmp = signal.spectrogram(wfBybp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(ttmp,ftmp,np.abs(ptmp),vr=vr,yr=yr,xr=xr, yscale=ys,xlabel='time(s)',ylabel='power\nf(Hz)')

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









#Test calibrations of VLF data 

#tst = wfAxbp[]

t0x = 162.
t1x = 167.
goo = np.where((tdatx >= t0x) & (tdatx <= t1x))

plt.plot(tdatx[goo[0]], wfAxbp[goo[0]])
plt.plot(tdatx[goo[0]], wfBxbp[goo[0]])

plt.plot(tdatx[goo[0]], wfBybp[goo[0]])
plt.plot(tdatx[goo[0]], wfAybp[goo[0]])

plt.plot(tdatx[goo[0]], wf12bp[goo[0]])
plt.plot(tdatx[goo[0]], wf34bp[goo[0]])


plt.plot(tdatx[goo[0]], wfAxhp[goo[0]])
plt.plot(tdatx[goo[0]], wfBxhp[goo[0]])
plt.plot(tdatx[goo[0]], wfByhp[goo[0]])
plt.plot(tdatx[goo[0]], wfAyhp[goo[0]])

plt.plot(tdatx[goo[0]], wf12hp[goo[0]])
plt.plot(tdatx[goo[0]], wf34hp[goo[0]])

