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
vAstrx = 'VLF13D'
vBstrx = 'VLF24D'  #-->42
polarity_xA = 1
polarity_xB = -1
##y-hat' direction from phase analysis uses these two components
vAstry = 'VLF32D'
vBstry = 'VLF41D'  #-->14
polarity_yA = 1
polarity_yB = -1
##Endurance has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
receiver_spacing_xhat = 2.27 #m 
receiver_spacing_yhat = 2.27 #m


#----------------------------------------------------------------
#Use one short spaced receiver and one long
#----------------------------------------------------------------
##x-hat' direction from phase analysis uses these two components
#vAstrx = 'VLF13D'
#vBstrx = 'VLF34D'  #-->43
#polarity_xA = 1
#polarity_xB = -1
##y-hat' direction from phase analysis uses these two components
#vAstry = 'VLF32D'
#vBstry = 'VLF41D' #-->14
#polarity_yA = 1
#polarity_yB = -1
##Endurance has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
#receiver_spacing_xhat = 2.27/2 #m 
#receiver_spacing_yhat = 2.27 #m

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



#Reduce waveforms to time of interest only for faster calculations
wfAxz = wfAx[np.where((tdatx > 875) & (tdatx < 900))]
wfAyz = wfAy[np.where((tdaty > 875) & (tdaty < 900))]
wfBxz = wfBx[np.where((tdatx > 875) & (tdatx < 900))]
wfByz = wfBy[np.where((tdaty > 875) & (tdaty < 900))]


#---------------------------------------------------------------------------




cohmin = 0.1  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
vr = [-45,-20]
ys = 'linear'
kr = [-5,5]
tz = 888.5
nsec = 1  #coherence spectra will be averaged over this time





#tchunk should be < nsec 

tchunk = 0.1  #delta-time (sec) for each FFT
nperseg = 1024 #determine freq resolution

#tchunk = 0.05  #delta-time (sec) for each time chunk to divide up the spectra into
nchunks = int(np.ceil((wfAx.size/fs)/tchunk)) #number of chunks in ENTIRE timerange
ndatapts = tchunk * fs #Make sure this number is >> nperseg
print(ndatapts/nperseg)


cohx2, phasex2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAxz,wfBxz,tdatx,fs,tchunk,coh_min=cohmin,nperseg=nperseg)
cohy2, phasey2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAyz,wfByz,tdaty,fs,tchunk,coh_min=cohmin,nperseg=nperseg)

phasearrx2 = phasex2 
coharrx2 = cohx2
phasearry2 = phasey2 
coharry2 = cohy2



t0z = 840
t1z = 845
wfAxz2 = wfAx[np.where((tdatx > t0z) & (tdatx < t1z))]
wfAyz2 = wfAy[np.where((tdaty > t0z) & (tdaty < t1z))]
wfBxz2 = wfBx[np.where((tdatx > t0z) & (tdatx < t1z))]
wfByz2 = wfBy[np.where((tdaty > t0z) & (tdaty < t1z))]



cohz, anglez, fz = correlation_analysis.cross_spectral_density(wfAxz2,wfBxz2,fs,nperseg=2048,plotshow=False)

fig, axs = plt.subplots(2)
axs[0].plot(fz, cohz)
axs[1].plot(fz, np.degrees(anglez))
axs[0].set_xscale('linear')
axs[1].set_xscale('linear')
axs[0].set_xlim(4000,8000)
axs[1].set_xlim(4000,8000)
axs[0].set_ylim(0,1.5)


"""
#Reduce arrays to desired timerange
pavgx2, phasearrx2, tarr_phasex2 = ps.slice_spectrogram(tz,tchunks2,phasex2,nsec)
cavgx2, coharrx2, tarr_cohx2 = ps.slice_spectrogram(tz,tchunks2,cohx2,nsec)
#powavgAx2, powarrAx2, tarr_powAx2 = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)
#powavgBx2, powarrBx2, tarr_powBx2 = ps.slice_spectrogram(tz,tspecx,np.abs(powercBx),nsec)

pavgy2, phasearry2, tarr_phasey2 = ps.slice_spectrogram(tz,tchunks2,phasey2,nsec)
cavgy2, coharry2, tarr_cohy2 = ps.slice_spectrogram(tz,tchunks2,cohy2,nsec)
#powavgAy2, powarrAy2, tarr_powAy2 = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)
#powavgBy2, powarrBy2, tarr_powBy2 = ps.slice_spectrogram(tz,tspecy,np.abs(powercBy),nsec)
"""



#Turn the phase values into k-values 
kx2 = [np.radians(i) / receiver_spacing_xhat for i in pavgx2]
ky2 = [np.radians(i) / receiver_spacing_yhat for i in pavgy2]
#and then into wavelength perp to Bo
kmag = [np.sqrt(kx2[i]**2 + ky2[i]**2) for i in range(len(kx2))]
wl2 = [2*np.pi/i for i in kmag]

#----------------------------------------------------------------------------------------
#Calculate the fractional difference b/t different spectra to identify artificial waves
#----------------------------------------------------------------------------------------

#ptmp_diff = np.abs(np.abs(powarrAx2) - np.abs(powarrBx2))
#ptmp_sum = np.abs(powarrAx2) + np.abs(powarrBx2)
#fracdiffx = ptmp_diff/ptmp_sum
#ptmp_diff = np.abs(np.abs(powarrAy2) - np.abs(powarrBy2))
#ptmp_sum = np.abs(powarrAy2) + np.abs(powarrBy2)
#fracdiffy = ptmp_diff/ptmp_sum


ptmp_diff = np.abs(np.abs(powercAx) - np.abs(powercBx))
ptmp_sum = np.abs(powercAx) + np.abs(powercBx)
fracdiffx = ptmp_diff/ptmp_sum
ptmp_diff = np.abs(np.abs(powercAy) - np.abs(powercBy))
ptmp_sum = np.abs(powercAy) + np.abs(powercBy)
fracdiffy = ptmp_diff/ptmp_sum

goo, fracdiffx2, goo2 = ps.slice_spectrogram(tz,tspecx,fracdiffx,nsec)
goo, fracdiffy2, goo2 = ps.slice_spectrogram(tz,tspecx,fracdiffy,nsec)



#---------------------------------------------------------------------------------------

yr = [10,8000]
vr=[-65,-20]
xrspec = [tz,tz+nsec]
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx


#Plot values from Method 2
fig,axs = plt.subplots(7,2, figsize=(11,8))#, sharex=True)  #,gridspec_kw={'height_ratios':[1,1,1,1,1,1,1,1,1]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tspecxz,fspecx,np.abs(powarrAx2),vr=vr,xr=xrspec,yr=yr,yscale='linear',ax=axs[0,0],ylabel="power\nspec (x')\nHz", title=titlegoo)
axs[0,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tarr_cohx2,freqs2,coharrx2**2,vr=[0,1],zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[1,0],ylabel="coh**2\n(x')\nHz")
axs[1,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tarr_phasex2,freqs2,np.degrees(phasearrx2),vr=[-180,180], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[2,0],ylabel="phase(deg)\n(x')\nHz", cmap='twilight_shifted')
axs[2,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tarr_cohy2,freqs2,coharry2**2,vr=[0,1],zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[3,0],ylabel="coh**2\n(y')\nHz")
axs[3,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tarr_phasey2,freqs2,np.degrees(phasearry2),vr=[-180,180], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[4,0],ylabel="phase(deg)\n(y')\nHz", cmap='twilight_shifted')
axs[4,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecxz,fspecx,fracdiffx2,vr=[0,1],xr=xrspec,yr=yr,zscale='linear',yscale='linear',ax=axs[5,0],ylabel="% diff (x')\nHz")
axs[5,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecxz,fspecx,fracdiffy2,vr=[0,1],xr=xrspec,yr=yr,zscale='linear',yscale='linear',ax=axs[6,0],xlabel='time(sec)',ylabel="% diff (y')\nHz")

for i in range(3): axs[i,0].axvline(tz, linestyle='--')
for i in range(3): axs[i,0].axvline(tz + nsec, linestyle='--')
axs[0,1].plot(fspecx,powarrAx2)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].plot(fspecx,powavgAx2,'.',color='black')
axs[0,1].set_xlim(yr)
axs[0,1].set_ylim(0,np.nanmax(powarrAx2))
axs[0,1].set_ylabel('powermax')

axs[1,1].plot(freqs2,coharrx2**2)
axs[1,1].plot(freqs2,cavgx2,'.',color='black')
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_xlim(yr)
axs[1,1].set_ylim(0,1)
axs[1,1].set_ylabel("coh**2 (x')")

axs[2,1].plot(freqs2,np.degrees(phasearrx2))
axs[2,1].plot(freqs2,pavgx2,'.',color='black')
axs[2,1].get_xaxis().set_visible(False)
axs[2,1].set_xlim(yr)
axs[2,1].set_ylim(-180,180)
axs[2,1].set_ylabel("phase (x')")

axs[3,1].plot(freqs2,coharry2**2)
axs[3,1].plot(freqs2,cavgy2,'.',color='black')
axs[3,1].get_xaxis().set_visible(False)
axs[3,1].set_xlim(yr)
axs[3,1].set_ylim(0,1)
axs[3,1].set_ylabel("coh**2 (y')")

axs[4,1].plot(freqs2,np.degrees(phasearry2))
axs[4,1].plot(freqs2,pavgy2,'.',color='black')
axs[4,1].get_xaxis().set_visible(False)
axs[4,1].set_xlim(yr)
axs[4,1].set_ylim(-180,180)
axs[4,1].set_ylabel("phase (y')")

axs[5,1].plot(freqs2,kx2,'.',freqs2,ky2,'.')
axs[5,1].get_xaxis().set_visible(False)
axs[5,1].set_xlim(yr)
axs[5,1].set_ylim(min(np.nanmin(kx2),np.nanmin(ky2)),max([np.nanmax(kx2),np.nanmax(ky2)]))
axs[5,1].set_ylabel("kx'(blue)\nky'(orange)")

axs[6,1].plot(freqs2,wl2,'.',color='black')
axs[6,1].set_xlim(yr)
axs[6,1].set_yscale('log')
axs[6,1].set_ylim(1,1000)
axs[6,1].set_ylabel("wavelength(m)\nfrom |k|")
axs[6,1].set_xlabel('Hz')




#------------------------------------------------------------
#Plot final results
#------------------------------------------------------------

titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx
xr = [tz-3*nsec,tz+3*nsec]

krplot = [-2,2]

fig,axs = plt.subplots(7,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1,1,1,1],'width_ratios':[1,0.5]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercAx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[0,0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
axs[0,0].get_xaxis().set_visible(False)

ps.plot_spectrogram(tspecx,fspecx,cohx**2,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[1,0],xlabel='time(s)',ylabel="Coh**2\n(x')\nf(Hz)")
axs[1,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecx,fspecx,np.degrees(phasex),vr=[-180,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[2,0],xlabel='time(s)',ylabel="Phase\n(x')\nf(Hz)",cmap='twilight_shifted')
axs[2,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecx,fspecx,cohy**2,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[3,0],xlabel='time(s)',ylabel="Coh**2\n(y')\nf(Hz)")
axs[3,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecx,fspecx,np.degrees(phasey),vr=[-180,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[4,0],xlabel='time(s)',ylabel="Phase\n(y')\nf(Hz)",cmap='twilight_shifted')
axs[4,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecx,fspecx,fracdiffx,vr=[0,1],xr=xr,yr=yr,zscale='linear',yscale='linear',ax=axs[5,0],xlabel='time(sec)',ylabel="% diff\n(x')\nf(Hz)")
axs[5,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecx,fspecx,fracdiffy,vr=[0,1],xr=xr,yr=yr,zscale='linear',yscale='linear',ax=axs[6,0],xlabel='time(sec)',ylabel="% diff\n(y')\nf(Hz)")

for i in range(7): axs[i,0].axvline(tz, linestyle='--')
for i in range(7): axs[i,0].axvline(tz+nsec, linestyle='--')

ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1,1],minzval=0,maxzval=1,xlabel="kx'(rad/m)",ylabel=xptitle+'\n(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1,1],minzval=-120,maxzval=10,alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[3,1],minzval=0,maxzval=1,xlabel="ky'(rad/m)",ylabel=yptitle+'\n(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[3,1],minzval=-120,maxzval=10,alpha=0.5)
#Plot the limiting k-vector value where short wavelength effects start to occur. 
#i.e. location when wavelength equals about twice the interferometry receiver spacing
klimx = 2*np.pi / (2 * receiver_spacing_xhat)
klimy = 2*np.pi / (2 * receiver_spacing_yhat)
wlimx = (2 * receiver_spacing_xhat)
wlimy = (2 * receiver_spacing_yhat)
axs[1,1].axvline(klimx,linestyle="--")
axs[1,1].axvline(-klimx,linestyle="--")
axs[3,1].axvline(klimy,linestyle="--")
axs[3,1].axvline(-klimy,linestyle="--")

#Oplot the results from the initial (1d) analysis (Method 2)
axs[1,1].plot(kx2,freqs2,'*',color='black')
axs[3,1].plot(ky2,freqs2,'*',color='black')
axs[5,1].plot(wl1,fspecx,'.',color='black',markersize=2)
axs[5,1].plot(wl2,freqs2,'*',color='black')
axs[5,1].set_ylim(yr)
axs[5,1].set_xscale('log')
axs[5,1].set_xlim(1,1000)
axs[5,1].set_xlabel("wavelength(m)\n(from |k|)")
axs[5,1].set_ylabel("(Hz)")
axs[5,1].axvline(wlimx,linestyle="--")
axs[5,1].axvline(wlimy,linestyle="--")

fig.delaxes(axs[0,1])
fig.delaxes(axs[2,1])
fig.delaxes(axs[4,1])
fig.delaxes(axs[6,1])










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



