#Make nice looking spectral plots for GIRAFF Bernstein waves on 380

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
#from scipy import signal
import numpy as np 
import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
#import plasma_params_get_flhr_freq as dflh
#import pyIGRF
import pickle
#import scipy.io as sio
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
#from scipy.interpolate import interp1d
import gir_load_data as gld






#--------------------------------------------------------------------------------------------------
#Read in polarization data of Bernstein waves from IDL save file
pathz = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/polarization_from_idl/'
if pld == '380':
    idldat1 = readsav(pathz + 'Giraff_380_pol_from_idl_fft=8192_100-300sec.sav')
    #idldat1 = readsav(pathz + 'Giraff_380_pol_from_idl_fft=8192_300-550sec.sav')
    #idldat1 = readsav(pathz + 'Giraff_380_pol_from_idl_fft=4096_400-500sec–EDC_not_gp_calibrated.sav')
elif pld == '381':
    #idldat1 = readsav(pathz + 'Giraff_381_pol_from_idl_fft=4096_100-300sec.sav')
    idldat1 = readsav(pathz + 'Giraff_381_pol_from_idl_fft=4096_300-550sec.sav')






data_pow = idldat1['data_pow']
data_elip = idldat1['data_elip']
data_pol = idldat1['data_pol']
data_hel = idldat1['data_hel']
times_elip = idldat1['times_elip']
freqs_elip = idldat1['freqs_elip']


#reduce ellipticity values to only those with a high deg polarization
data_elip2 = np.copy(data_elip)
for i in range(np.shape(data_pol)[0]):
      goodv = np.where(data_pol[i,:] < 0.7)[0]
      data_elip2[i,goodv] = np.nan

#--------------------------------------------------------------------------------------------------


#Get altitude data
ephem = gld.load_ephemeris(pld)
alts = ephem['altitude']
time_ephem = ephem['time']




#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))[0]



#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()


#Load SLP data
lp = gld.load_slp(pld)
lptimes = lp[0]
lpvals = lp[1]
#plt.plot(lptimes, lpvals)

#create sonogram of Langmuir Probe data
nfft=1024
fspecLP, tspecLP, powerLP, fsLP = fftspec.fft_spectrum_piecewise(lptimes, lpvals, fs_thres=0.1, nfft=nfft, noverlap=8)




#Reduce waveforms to times b/t 300-500 sec. 
goodt = np.where(tdat > 100)[0]
wf12 = wf12[goodt]
wf34 = wf34[goodt]
tdat = tdat[goodt]





#magnetic field 
#magv = GFL(pld,'mag')
#mag = magv.load_data()
#Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
#Bot = mag[3]

#Load IGRF data 
Bot,Bo = gld.load_igrf(pld)
#plt.plot(tigrf, Boigrf)



#---------------------------------
#Overall spectral plot (Fig 1)
#---------------------------------

#nfft=2048
nfft=1024
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=8)
fspec, tspec, powerc34, fs = fftspec.fft_spectrum_piecewise(tdat, wf34, fs_thres=0.1, nfft=nfft, noverlap=8)

nfft = 1024
freqs_fin, tcenter_fin, csd_fin, coh_fin, phase_fin, fs_fin, spec_fin1, spec_fin2 =  ca.csd_spectrum_piecewise(tdat[tdat > 100], wf34[tdat > 100], wf12[tdat > 100], nfft=nfft, noverlap=8, fs_thres=0.2)


phase_fin = phase_fin*(180/3.14)

pow = np.abs(spec_fin1)
minpow = 5e-10
for i in range(len(tcenter_fin)):
      goo = np.where(pow[:,i] < minpow)[0]
      phase_fin[goo,i] = np.nan




#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

#Gyrofrequencies (Hz)
fce = 28 * Bo   
#Interpolate high-cadence fce to low-cadence fpe
fce = np.interp(tspec, Bot, fce)
fcH = fce / 1836
fcO = fcH / 16



#-------------------------------------
#Filter polarization values (from Chaston's code)
#-------------------------------------

pmin = -140
polmin = 0.7

ptmp = 10.0 * np.log10(data_pow)
elip_tmp = data_elip2.copy()
hel_tmp = data_hel.copy()
elip_tmp[ptmp < pmin] = np.nan
elip_tmp[data_pol < polmin] = np.nan

hel_tmp[ptmp < pmin] = np.nan
hel_tmp[data_pol < polmin] = np.nan

#-------------------------------------
#Filter coherence/phase values
#-------------------------------------

pmin = -140
cohmin = 0.7

ptmp = 10.0 * np.log10(spec_fin1)
phase_tmp = phase_fin.copy()
coh_tmp = coh_fin.copy()
phase_tmp[ptmp < pmin] = np.nan
phase_tmp[coh_tmp < cohmin] = np.nan
coh_tmp[ptmp < pmin] = np.nan
coh_tmp[coh_tmp < cohmin] = np.nan






#--------------------------------------------
#Overall plots
#--------------------------------------------

#fspecLP, tspecLP, powerLP, fsLP


#xr = [100,300]
xr = [300,500]
vr = [-80,-60]
#vr_lp = [-250,-150] #380
vr_lp = [-230,-220] #381

yr = [500,16000]
fig, axs = plt.subplots(3,3, figsize=(16,9))
plt.title('gir_analysis_Bernstein_spectra.py')
plt.title('GIRAFF ' + pld + ' Bernstein (gir_analysis_Bernstein_spectra.py)')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[0,0])
axs[0,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_tmp,vr=[-120,120],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",ax=axs[1,0])
axs[1,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_tmp,vr=[cohmin,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,0])
ps.plot_spectrogram(times_elip,freqs_elip,hel_tmp,vr=[0,0.6],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="helicity (V12,V34)\nfreq(Hz)",ax=axs[0,1])
axs[0,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,elip_tmp,vr=[-1,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="ellipticity (V12,V34)\nfreq(Hz)",ax=axs[1,1],plot_kwargs={'cmap':'seismic'})
axs[1,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_pol,vr=[polmin,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="polarization (V12,V34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,1])

ps.plot_spectrogram(tspecLP, fspecLP, powerLP, vr=vr_lp, yscale='linear', yr=yr, xr=xr, ylabel="Langmuir Probe power\nfreq(Hz); dB", ax=axs[0,2])

fig.tight_layout(pad=2)
fig.subplots_adjust(left=0.15)

for i in range(3):
     for j in range(3):
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='black',s=12)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
            # Set font size to 8 pt for all text elements
            axs[i,j].tick_params(axis='both', labelsize=8)
            axs[i,j].xaxis.label.set_fontsize(8)
            axs[i,j].yaxis.label.set_fontsize(8)
            axs[i,j].plot(tspec,fcH,color='magenta')
            axs[i,j].plot(tspec,2*fcH,color='magenta')
            axs[i,j].plot(tspec,3*fcH,color='magenta')
            axs[i,j].plot(tspec,4*fcH,color='magenta')
            axs[i,j].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)

minvert = 0
for i in range(3):
      for j in range(3):
            ax2 = axs[i,j].twinx()
            ax2.plot(time_ephem, alts, color='blue', linewidth=2)
            ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
            ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
            ax2.yaxis.set_label_position("left")
            ax2.tick_params(axis='y', labelleft=True, labelright=False)
            ax2.spines['left'].set_position(('outward', 60))












#xr = [400,500]
#xr = [450,460]
xr = [400,420]
vr = [-80,-60]
fig, axs = plt.subplots(3,2, figsize=(16,9))
plt.title('GIRAFF ' + pld + ' Bernstein (gir_analysis_Bernstein_spectra.py)')
yr = [1000,8000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,0])
axs[0,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",ax=axs[1,0])
axs[1,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,0])
ps.plot_spectrogram(times_elip,freqs_elip,data_hel,vr=[0,0.6],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="helicity (V12,V34)\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,1])
axs[0,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,vr=[-1,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="ellipticity (V12,V34)\nfreq(Hz)",ax=axs[1,1])
axs[1,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_pol,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="polarization (V12,V34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,1])
fig.tight_layout(pad=0)

minvert = 0
for i in range(3):
      for j in range(2):
            if j == 1:
                axs[i,1].plot(tspec,fcH,color='magenta')
                axs[i,1].plot(tspec,2*fcH,color='magenta')
                axs[i,1].plot(tspec,3*fcH,color='magenta')
                axs[i,1].plot(tspec,4*fcH,color='magenta')
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='white',s=12)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
            axs[i,j].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)



#plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)







xr = [410,430]
xr = [460,470]
#xr = [452,455]
#xr = [450,460]
#xr = [400,420]
vr = [-80,-60]
fig, axs = plt.subplots(3,2, figsize=(16,9))
plt.title('GIRAFF ' + pld + ' Bernstein (gir_analysis_Bernstein_spectra.py)')
yr = [2000,6000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,0])
axs[0,0].set_xticklabels([])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1,0])
axs[1,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[0,180],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",ax=axs[2,0])

ps.plot_spectrogram(times_elip,freqs_elip,data_hel,vr=[0,0.6],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="helicity (V12,V34)\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,1])
axs[0,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,vr=[-1,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="ellipticity (V12,V34)\nfreq(Hz)",ax=axs[1,1])
axs[1,1].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_pol,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="polarization (V12,V34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,1])
fig.tight_layout(pad=0)

minvert = 0
for i in range(3):
      for j in range(2):
            axs[i,j].plot(tspec,fcH,color='magenta')
            axs[i,j].plot(tspec,2*fcH,color='magenta')
            axs[i,j].plot(tspec,3*fcH,color='magenta')
            axs[i,j].plot(tspec,4*fcH,color='magenta')
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='white',s=12)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
            axs[i,j].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)




