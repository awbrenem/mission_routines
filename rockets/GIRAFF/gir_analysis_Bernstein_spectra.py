#Make nice looking spectral plots for GIRAFF Bernstein waves on 380

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from scipy import signal
import numpy as np 
import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import plasma_params_get_flhr_freq as dflh
import pyIGRF
import pickle
#import scipy.io as sio
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
from scipy.interpolate import interp1d


pld = '380'



#--------------------------------------------------------------------------------------------------
#Read in polarization data of Bernstein waves from IDL save file
pathz = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/polarization_from_idl/'

idldat1 = readsav(pathz + 'Giraff_380_pol_from_idl_fft=4096_400-500sec–EDC_not_gp_calibrated.sav')


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


#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))[0]



#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()


#Reduce waveforms to times b/t 400-500 sec. 
goodt = np.where(tdat > 400)[0]
wf12 = wf12[goodt]
wf34 = wf34[goodt]
tdat = tdat[goodt]





"""
#Get lower hybrid freq from IRI
alts = np.arange(100,300,10)
flhr_z = np.zeros(len(alts))
timesz = np.zeros(len(alts))
for i in range(len(alts)):     
      Boz = pyIGRF.igrf_value(glat, glon, alts[i], 2022)[6]
      iriz = end_data_loader.load_iri(alt=alts[i])
      densz = iriz['Ne(m-3)']/(100**3)
      fcez = 28*Boz
      flhr_z[i] = dflh.flhr_IonMassFractions(densz, fcez, 0.01*iriz['H_ions'], 0.01*iriz['O_ions'])
      timesz[i] = iriz['times_downleg']
#plt.plot(alts,flhr_z)
"""




#magnetic field 
magv = GFL(pld,'mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]




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




