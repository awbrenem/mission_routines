"""
Run interferometry analysis on DC density structures from GIRAFF
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

-----------------------------------------------
GIRAFF (NOTE: below have been updated from GIRAFF-sketch_wcallouts.png)
V1 = 0 deg
V2 = 180 deg
V3 = 90 deg
V4 = 270 deg
BELOW is view from aft looking forward.
                      
                         V4
                      /  |  \             
                    @    |     @         
                  /      |       \        
                 V2-------O--------V1  (x-hat)
                  \      |       /
           (y-hat') @    |     @  (x-hat')
                      \  |  /
                         V3
                      (y-hat)

                      
@-points represent the centers of potential of the interferometry diagonals          


Coordinate system (system of input test wave)
x-hat --> E12 = V1 - V2 direction (positive to right)
y-hat --> E34 = V3 - V4 direction (positive downwards)

This code outputs the spectrum of k-values in the kx' and ky' directions, where
x'-hat --> Ex' = V1V3z - V4V2z (45 deg inclined from x-hat)
y'-hat --> Ey' = V3V2z - V1V4z (45 deg inclined from y-hat)


#---interferometry pairs
#NOTE: positive sense of phase defined as pointing towards center of potential of "wfA"
#e.g. y-hat':  wfA = wf32; wfB = wf14
#So, for GIRAFF the interferometry pairs are 
#x-hat': vAstr = 'VLF13D' and vBstr = 'VLF42D' (from -1*'VLF24D' channel)
#y-hat': vAstr = 'VLF32D' and vBstr = 'VLF14D' (from -1*'VLF41D' channel)


"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import numpy as np 
import interferometry_routines as interf
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import fft_spectrum_piecewise as fftspec
import correlation_analysis as ca


#-------------------------------------------------------------
#Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using GIRAFF data
#-------------------------------------------------------------


#----------------------------------------------------------------
#Use long spaced receivers
#----------------------------------------------------------------
##x-hat' direction from phase analysis uses these two components
vAstrx = 'VLF13DF'
vBstrx = 'VLF24DF'  #-->42
polarity_xA = 1
polarity_xB = -1
##y-hat' direction from phase analysis uses these two components
vAstry = 'VLF32DF'
vBstry = 'VLF41DF'  #-->14
polarity_yA = 1
polarity_yB = -1
##GIRAFF has 3.25m booms, so the effective length of the diagonals is = 3.25*cos(45) = 2.3 m
receiver_spacing_xhat = 2.3 #m 
receiver_spacing_yhat = 2.3 #m


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
##GIRAFF has ~3m booms, so the effective length of the diagonals is = 3*cos(45) = 2.27
#receiver_spacing_xhat = 2.27/2 #m 
#receiver_spacing_yhat = 2.27 #m



pld = '381'

#Load GIRAFF waveform
vAx = GFL(pld,vAstrx)
vBx = GFL(pld,vBstrx)
wfAx, tdatx = vAx.load_data()
wfBx, tdatx = vBx.load_data()
vAy = GFL(pld,vAstry)
vBy = GFL(pld,vBstry)
wfAy, tdaty = vAy.load_data()
wfBy, tdaty = vBy.load_data()


#apply channel polarity
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
#Overview plot...
#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
#---------------------------------------------------------------------------

#******************************
#-->(SET)<--: size of following spectra AND of power spectra of k-value vs freq
nps = 512 
#******************************

fspec, tspec, powercAx, fsgoo = fftspec.fft_spectrum_piecewise(tdatx, wfAx, fs_thres=0.3, nfft=nps, noverlap=2)
fspec, tspec, powercBx, fsgoo = fftspec.fft_spectrum_piecewise(tdatx, wfBx, fs_thres=0.3, nfft=nps, noverlap=2)
fspec, tspec, powercAy, fsgoo = fftspec.fft_spectrum_piecewise(tdatx, wfAy, fs_thres=0.3, nfft=nps, noverlap=2)
fspec, tspec, powercBy, fsgoo = fftspec.fft_spectrum_piecewise(tdatx, wfBy, fs_thres=0.3, nfft=nps, noverlap=2)

ps.plot_spectrogram(tspec,fspec,np.abs(powercAx),vr=[-80,-20],yr=[100,50000],xr=[100,550], yscale='log')





#------------------------------------------------------------
#Reduce data to time range of interest 
# I'm using two timeranges. 
# (1) Arrays ending with "z": Phase/coherence calculations occur b/t [tz-dtslop, tz+nsec+dtslop], where dtslop is on the order of 10s sec. 
#     This is to ensure that the phase/coh calculations have enough data and is also a nice timerange for plots.
# (2) Arrays ending with "z2": I then average these (and extract slices) for the stricter timerange [tz, tz+nsec]
#------------------------------------------------------------


#******************************
#-->(SET)<--:
#Data to be considered lies b/t tz and tz + nsec

#----Hiss waves
#tz = 254
#nsec = 5
#tz =  403  #(403.8)  403.7 - 404
#nsec = 4
#----VLF waves associated with possible BBELF waves
tz = 268
nsec = 2
fplot_bracket = [10000,40000] #Used to set ylim for spectral plots (bracket wave range of interest)
cohmin = 0.6  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#******************************




#-------------------------------------------------------------
#Method 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using GIRAFF data
#-------------------------------------------------------------

#NOTE: + sense of phase defined as pointing towards center of potential of "powercA"

#******************************
#-->(SET)<--:
#Nval --> 2D blurring window. N must be large enough for the estimate to be robust, 
#         but small enough not to blur the transitions between different backscattering zones.
#         Try N~few 
Nval = 3
#Nval = 10
#******************************


#Reduce to only times needed plus some slop on either side
#...Later we will reduce versions of these to only the tz to tz + nsec
dtslop = 10
goodtslop = np.where((tdatx >= tz-dtslop) & (tdatx <= (tz+nsec+dtslop)))[0]
wfAxz = wfAx[goodtslop]
wfBxz = wfBx[goodtslop]
wfAyz = wfAy[goodtslop]
wfByz = wfBy[goodtslop]
tdatxz = tdatx[goodtslop]
tdatyz = tdaty[goodtslop]

fspecz, tspecz, csdxz, cohxz, phasexz, fsxz, powercAxz, powercBxz =  ca.csd_spectrum_piecewise(tdatxz, wfAxz, wfBxz, nfft=nps, noverlap=8, fs_thres=0.3)
fspecz, tspecz, csdyz, cohyz, phaseyz, fsyz, powercAyz, powercByz =  ca.csd_spectrum_piecewise(tdatyz, wfAyz, wfByz, nfft=nps, noverlap=8, fs_thres=0.3)

ps.plot_spectrogram(tspecz,fspecz,np.abs(powercAxz),vr=[-80,-20],yr=[300,50000],xr=[tz-dtslop,tz+nsec+dtslop], yscale='log')






#******************************
#-->(SET)<--:
#yrplot=[0,26000]
yrplot = fplot_bracket
xrplot=[tz-3,tz+3]
#******************************
fig, axs = plt.subplots(3)
ps.plot_spectrogram(tspecz,fspecz,np.abs(powercAxz),ax=axs[0],vr=[-90,-20],yr=yrplot,xr=xrplot)
ps.plot_spectrogram(tspecz,fspecz,cohxz,zscale='linear',vr=[0.8,1],ax=axs[1],yr=yrplot,xr=xrplot)
ps.plot_spectrogram(tspecz,fspecz,phasexz,zscale='linear',vr=[-3.14,3.14],ax=axs[2],yr=yrplot,xr=xrplot)



goo = cohxz < cohmin
cohxz[goo] = float("nan")
phasexz[goo] = float("nan")
goo = cohyz < cohmin
cohyz[goo] = float("nan")
phaseyz[goo] = float("nan")


#Reduce arrays to the timerange for strict analysis. 
powavgAx2, powercAxz2, tspecz2 = ps.slice_spectrogram(tz,tspecz,np.abs(powercAxz),nsec)
powavgBx2, powercBxz2, tspecz2 = ps.slice_spectrogram(tz,tspecz,np.abs(powercBxz),nsec)
cavgx2, cohxz2, ttmp = ps.slice_spectrogram(tz,tspecz,cohxz,nsec)
pavgx2, phasexz2, ttmp = ps.slice_spectrogram(tz,tspecz,phasexz,nsec)


powavgAy2, powercAyz2, tspecy2 = ps.slice_spectrogram(tz,tspecz,np.abs(powercAyz),nsec)
powavgBy2, powercByz2, tspecy2 = ps.slice_spectrogram(tz,tspecz,np.abs(powercByz),nsec)
cavgy2, cohyz2, ttmp = ps.slice_spectrogram(tz,tspecz,cohyz,nsec)
pavgy2, phaseyz2, ttmp = ps.slice_spectrogram(tz,tspecz,phaseyz,nsec)

fspecz2 = fspecz

#******************************
#-->(SET)<--:
klim_spec = [-1,1]
nkbins = 200
#******************************

#Get power spec of k vs freq
fkpowspecxz2, kvalsxz2, fvalsxz2, pmaxvalsxz2 = interf.inter_fvsk(np.abs(powercAxz2),tspecz2,fspecz2, 
                                         phasexz2,tspecz2,fspecz2,
                                         receiver_spacing_xhat,
                                         mean_max='max',
                                         nkbins=nkbins,
                                         klim=klim_spec)
fkpowspecyz2, kvalsyz2, fvalsyz2, pmaxvalsyz2 = interf.inter_fvsk(np.abs(powercAyz2),tspecz2,fspecz2, 
                                         phaseyz2,tspecz2,fspecz2,
                                         receiver_spacing_yhat,
                                         mean_max='max',
                                         nkbins=nkbins,
                                         klim=klim_spec)

#Turn k-values into wavelength
wl1 = np.zeros(len(fvalsxz2))
for i in range(len(fvalsxz2)):
    tmp = pmaxvalsxz2[i,:]
    kxgoo = kvalsxz2[np.where(tmp == 1)]
    tmp = pmaxvalsyz2[i,:]
    kygoo = kvalsyz2[np.where(tmp == 1)]
    kmaggoo = np.sqrt(kxgoo**2 + kygoo*2)
    if len(kmaggoo != 0): 
        wl1[i] = 2*np.pi/kmaggoo



#-------------------------------------------------------------
#Method 2: Do a freq (xaxis) vs k-value (yaxis) interferometry analysis
#-------------------------------------------------------------



#******************************
#-->(SET)<--:
#delta-time (sec) for each time chunk to average spectra over. This ultimately limits the time resolution.
#tchunk = 0.1
tchunk = 0.2
#choose based on desired freq resolution for coh and phase spectra
nperseg = 512
#******************************


#Turn the phase values into k-values 
kx2 = np.asarray(pavgx2) / receiver_spacing_xhat
ky2 = np.asarray(pavgy2) / receiver_spacing_yhat
#and then into wavelength perp to Bo
kmag = np.sqrt(kx2**2 + ky2**2) 
wl2 = 2*np.pi/kmag  #meters

plt.plot(fspecz, kmag)

#determine phase velocity vs freq 
#Vphase = 2*np.pi * freqs2 / np.asarray(kmag) / 1000  # km/s
Vphase = fspecz * wl2 / 1000
#plt.plot(freqs2, Vphase)

#----------------------------------------------------------------------------------------
#Calculate the fractional difference b/t different spectra to identify artificial waves
#----------------------------------------------------------------------------------------

ptmp_diff = np.abs(np.abs(powercAxz) - np.abs(powercBxz))
ptmp_sum = np.abs(powercAxz) + np.abs(powercBxz)
fracdiffx = ptmp_diff/ptmp_sum
ptmp_diff = np.abs(np.abs(powercAyz) - np.abs(powercByz))
ptmp_sum = np.abs(powercAyz) + np.abs(powercByz)
fracdiffy = ptmp_diff/ptmp_sum

goo, fracdiffx2, goo2 = ps.slice_spectrogram(tz,tspecz,fracdiffx,nsec)
goo, fracdiffy2, goo2 = ps.slice_spectrogram(tz,tspecz,fracdiffy,nsec)



#---------------------------------------------------------------------------------------

goodf = np.where((fspecz > fplot_bracket[0]) & (fspecz < fplot_bracket[1]))[0]


yr = [0,30000]
vr=[-90,-20]
xrspec = [tz,tz+nsec]
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx

plot_kwargs={'cmap':'turbo'}


#Plot values from Method 2
fig,axs = plt.subplots(9,2, figsize=(11,8))#, sharex=True)  #,gridspec_kw={'height_ratios':[1,1,1,1,1,1,1,1,1]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tspecz2,fspecz,np.abs(powercAxz2),vr=vr,xr=xrspec,yr=yr,yscale='linear',ax=axs[0,0],ylabel="power\nspec (x')\nHz", title=titlegoo)
axs[0,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,cohxz2**2,vr=[0,1],zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[1,0],ylabel="coh**2\n(x')\nHz")
axs[1,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,np.degrees(phasexz2),vr=[-180,180], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[2,0],ylabel="phase(deg)\n(x')\nHz",plot_kwargs={'cmap':'twilight_shifted'})
axs[2,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,cohyz2**2,vr=[0,1],zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[3,0],ylabel="coh**2\n(y')\nHz")
axs[3,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,np.degrees(phaseyz2),vr=[-180,180], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[4,0],ylabel="phase(deg)\n(y')\nHz",plot_kwargs={'cmap':'twilight_shifted'})
axs[4,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,fracdiffx2,vr=[0,1],xr=xrspec,yr=yr,zscale='linear',yscale='linear',ax=axs[5,0],ylabel="% diff (x')\nHz")
axs[5,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz2,fspecz,fracdiffy2,vr=[0,1],xr=xrspec,yr=yr,zscale='linear',yscale='linear',ax=axs[6,0],xlabel='time(sec)',ylabel="% diff (y')\nHz")

for i in range(3): axs[i,0].axvline(tz, linestyle='--')
for i in range(3): axs[i,0].axvline(tz + nsec, linestyle='--')
axs[0,1].plot(fspecz2,powercAxz2)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].plot(fspecz2,powavgAx2,'.',color='black')
axs[0,1].set_xlim(yr)
axs[0,1].set_ylim(0,np.nanmax(powercAxz2[goodf,:]))
axs[0,1].set_ylabel('powermax')

axs[1,1].plot(fspecz2,cohxz2**2)
axs[1,1].plot(fspecz2,cavgx2,'.',color='black')
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_xlim(yr)
axs[1,1].set_ylim(0,1)
axs[1,1].set_ylabel("coh**2\n(x')")

axs[2,1].plot(fspecz2,np.degrees(phasexz2))
axs[2,1].plot(fspecz2,np.degrees(pavgx2),'.',color='black')
axs[2,1].get_xaxis().set_visible(False)
axs[2,1].set_xlim(yr)
axs[2,1].set_ylim(-180,180)
axs[2,1].set_ylabel("phase\n(x')")

axs[3,1].plot(fspecz2,cohyz2**2)
axs[3,1].plot(fspecz2,cavgy2,'.',color='black')
axs[3,1].get_xaxis().set_visible(False)
axs[3,1].set_xlim(yr)
axs[3,1].set_ylim(0,1)
axs[3,1].set_ylabel("coh**2\n(y')")

axs[4,1].plot(fspecz2,np.degrees(phaseyz2))
axs[4,1].plot(fspecz2,np.degrees(pavgy2),'.',color='black')
axs[4,1].get_xaxis().set_visible(False)
axs[4,1].set_xlim(yr)
axs[4,1].set_ylim(-180,180)
axs[4,1].set_ylabel("phase\n(y')")

axs[5,1].plot(fspecz2,kx2,'.',fspecz,ky2,'.')
axs[5,1].get_xaxis().set_visible(False)
axs[5,1].set_xlim(yr)
axs[5,1].set_ylim(min(np.nanmin(kx2),np.nanmin(ky2)),max([np.nanmax(kx2),np.nanmax(ky2)]))
axs[5,1].set_ylabel("kx'(blue)\nky'(orange)")

axs[6,1].plot(fspecz2,wl2,'.',color='black')
axs[6,1].set_xlim(yr)
axs[6,1].set_yscale('linear')
axs[6,1].set_ylim(0,20)
axs[6,1].set_ylabel("wavelength\n(m)\nfrom |k|")
axs[6,1].set_xlabel('Hz')

axs[7,1].plot(fspecz2,Vphase,'.',color='black')
axs[7,1].set_xlim(yr)
axs[7,1].set_yscale('linear')
axs[7,1].set_ylim(0,2)
axs[7,1].set_ylabel("Vphase\n(km/s)")
axs[7,1].set_xlabel('Hz')

axs[8,1].plot(fspecz2,kmag/(2*3.14),'.',color='black')
axs[8,1].set_xlim(yr)
axs[8,1].set_yscale('linear')
axs[8,1].set_ylabel("|k|\n(1/m)")
axs[8,1].set_xlabel('Hz')

fig.delaxes(axs[7,0])
fig.delaxes(axs[8,0])


#------------------------------------------------------------
#Plot final results
#------------------------------------------------------------

titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx
xr = [tz-3*nsec,tz+3*nsec]

#******************************
#-->(SET)<--
krplot = [-1.2,1.2]
vr = [-85,-20]
ys = 'linear'
#******************************

fig,axs = plt.subplots(7,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1,1,1,1],'width_ratios':[1,0.5]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tspecz,fspecz,np.abs(powercAxz),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[0,0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
axs[0,0].get_xaxis().set_visible(False)

ps.plot_spectrogram(tspecz,fspecz,cohxz**2,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[1,0],xlabel='time(s)',ylabel="Coh**2\n(x')\nf(Hz)")
axs[1,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz,fspecz,np.degrees(phasexz),vr=[-180,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[2,0],xlabel='time(s)',ylabel="Phase\n(x')\nf(Hz)",plot_kwargs={'cmap':'twilight_shifted'})
axs[2,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz,fspecz,cohyz**2,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[3,0],xlabel='time(s)',ylabel="Coh**2\n(y')\nf(Hz)")
axs[3,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz,fspecz,np.degrees(phaseyz),vr=[-180,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[4,0],xlabel='time(s)',ylabel="Phase\n(y')\nf(Hz)",plot_kwargs={'cmap':'twilight_shifted'})
axs[4,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz,fspecz,fracdiffx,vr=[0,1],xr=xr,yr=yr,zscale='linear',yscale='linear',ax=axs[5,0],xlabel='time(sec)',ylabel="% diff\n(x')\nf(Hz)")
axs[5,0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tspecz,fspecz,fracdiffy,vr=[0,1],xr=xr,yr=yr,zscale='linear',yscale='linear',ax=axs[6,0],xlabel='time(sec)',ylabel="% diff\n(y')\nf(Hz)")

for i in range(7): axs[i,0].axvline(tz, linestyle='--')
for i in range(7): axs[i,0].axvline(tz+nsec, linestyle='--')

ps.plot_spectrogram(kvalsxz2/(2*3.14),fvalsxz2,pmaxvalsxz2,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1,1],minzval=0,maxzval=1,xlabel="kx'(1/m)",ylabel=xptitle+'\n(Hz)',plot_kwargs={'cmap':'Greys'})
ps.plot_spectrogram(kvalsxz2/(2*3.14),fvalsxz2,fkpowspecxz2,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1,1],minzval=-120,maxzval=10,plot_kwargs2={'origin':'lower','alpha':0.5,'interpolation':'nearest','aspect':'auto'})
ps.plot_spectrogram(kvalsyz2/(2*3.14),fvalsyz2,pmaxvalsyz2,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[3,1],minzval=0,maxzval=1,xlabel="ky'(1/m)",ylabel=yptitle+'\n(Hz)',plot_kwargs={'cmap':'Greys'})
ps.plot_spectrogram(kvalsyz2/(2*3.14),fvalsyz2,fkpowspecyz2,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[3,1],minzval=-120,maxzval=10,plot_kwargs2={'origin':'lower','alpha':0.5,'interpolation':'nearest','aspect':'auto'})
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
axs[1,1].plot(kx2,fspecz,'*',color='black')
axs[3,1].plot(ky2,fspecz,'*',color='black')
axs[5,1].plot(wl1,fspecz,'.',color='black',markersize=2)
axs[5,1].plot(wl2,fspecz,'*',color='black')
axs[5,1].set_ylim(yr)
axs[5,1].set_xscale('log')
axs[5,1].set_xlim(1,100)
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
ps.plot_spectrogram(kvalsxz2/(2*3.14),fvalsxz2,pmaxvalsxz2,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[0],minzval=0,maxzval=1,xlabel="kx'(1/m)",ylabel=xptitle+'\nf(Hz)',plot_kwargs={'cmap':'Greys'},title=titlegoo)
ps.plot_spectrogram(kvalsxz2/(2*3.14),fvalsxz2,fkpowspecxz2,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[0],minzval=-120,maxzval=10,plot_kwargs2={'origin':'lower','alpha':0.5,'interpolation':'nearest','aspect':'auto'})
ps.plot_spectrogram(kvalsyz2/(2*3.14),fvalsyz2,pmaxvalsyz2,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1],minzval=0,maxzval=1,xlabel="ky'(1/m)",ylabel=yptitle+'\nf(Hz)',plot_kwargs={'cmap':'Greys'})
ps.plot_spectrogram(kvalsyz2/(2*3.14),fvalsyz2,fkpowspecyz2,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1],minzval=-120,maxzval=10,plot_kwargs2={'origin':'lower','alpha':0.5,'interpolation':'nearest','aspect':'auto'})


axs[0].axvline(klimx, linestyle='--')
axs[0].axvline(-klimx, linestyle='--')
axs[1].axvline(klimy, linestyle='--')
axs[1].axvline(-klimy, linestyle='--')

#Oplot the results from the initial (1d) analysis
axs[0].plot(kx2,fspec,'*',color='black')
axs[1].plot(ky2,fspec,'*',color='black')

axs[2].plot(wl1,fspec,'.',color='black',markersize=2)
axs[2].plot(wl2,fspec,'*',color='black')
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


