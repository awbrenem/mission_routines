#Make nice looking spectral plots for GIRAFF

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/PFISR/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from scipy import signal
import numpy as np 
import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import pickle
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
from scipy.interpolate import interp1d
import gir_load_data as gld
import pfisr_load_data as pfld
from datetime import datetime, timezone



pld = '381'
plot_harmonics = True



#Load IGRF data 
Bot,Bo = gld.load_igrf(pld)


#--------------------------------------------------------------------------------------------------
#Read in polarization data of GIRAFF waves from IDL save file
pathz = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/polarization_from_idl/'

if pld == '381':
      idldat1 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=4096_100-300sec.sav')
      idldat2 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=4096_300-550sec.sav')
if pld == '380':
      idldat1 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=8192_100-300sec.sav')
      idldat2 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=8192_300-550sec.sav')



times_pow = np.concatenate((idldat1['times_pow'],idldat2['times_pow']))
freqs_pow = np.concatenate((idldat1['freqs_pow'],idldat2['freqs_pow']))
data_pow = np.concatenate((idldat1['data_pow'],idldat2['data_pow']),axis=1)

times_pol = np.concatenate((idldat1['times_pol'],idldat2['times_pol']))
freqs_pol = np.concatenate((idldat1['freqs_pol'],idldat2['freqs_pol']))
data_pol = np.concatenate((idldat1['data_pol'],idldat2['data_pol']),axis=1)

times_elip = np.concatenate((idldat1['times_elip'],idldat2['times_elip']))
freqs_elip = np.concatenate((idldat1['freqs_elip'],idldat2['freqs_elip']))
data_elip = np.concatenate((idldat1['data_elip'],idldat2['data_elip']),axis=1)

times_hel = np.concatenate((idldat1['times_hel'],idldat2['times_hel']))
freqs_hel = np.concatenate((idldat1['freqs_hel'],idldat2['freqs_hel']))
data_hel = np.concatenate((idldat1['data_hel'],idldat2['data_hel']),axis=1)



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






#Load PFISR data for comparison
pftype = '20min'

if pld == '380':
    fn = '20250209.002_lp_'+pftype+'-fitcal.h5'
    tslice = datetime(2025,2,9,8,35,00, tzinfo=timezone.utc)  #GIRAFF 36.380 launch time
if pld == '381':
    fn = '20250202.002_lp_'+pftype+'-fitcal.h5'
    tslice = datetime(2025,2,2,7,7,00, tzinfo=timezone.utc)  #GIRAFF 36.381 launch time

colors = ['blue','orange','green','red','purple','maroon','magenta','navy','gold','cyan','darkblue','orangered','teal']


filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/pfisr/36.'+pld+'/' + fn

type = 'Ne_Mod'  #Density in m-3
pfisrTimes, pfisrAlts, pfisrNeMod = pfld.pfisr_load_data(filename, type, plot_test=0)

pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 1, plot_test=0)
pvstBeam1 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 2, plot_test=0)
pvstBeam2 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 3, plot_test=0)
pvstBeam3 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 4, plot_test=0)
pvstBeam4 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 5, plot_test=0)
pvstBeam5 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 6, plot_test=0)
pvstBeam6 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 7, plot_test=0)
pvstBeam7 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 8, plot_test=0)
pvstBeam8 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 9, plot_test=0)
pvstBeam9 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 10, plot_test=0)
pvstBeam10 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 11, plot_test=0)
pvstBeam11 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 12, plot_test=0)
pvstBeam12 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 13, plot_test=0)
pvstBeam13 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)





#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))[0]


#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()


#---------------------------------
#Overall spectral plot (Fig 1)
#---------------------------------

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=8)
fspec, tspec, powerc34, fs = fftspec.fft_spectrum_piecewise(tdat34, wf34, fs_thres=0.1, nfft=nfft, noverlap=8)

nfft = 4096
freqs_fin, tcenter_fin, csd_fin, coh_fin, phase_fin, fs_fin, spec_fin1, spec_fin2 =  ca.csd_spectrum_piecewise(tdat[tdat > 100], wf34[tdat > 100], wf12[tdat > 100], nfft=nfft, noverlap=8, fs_thres=0.2)

phase_fin = phase_fin*(180/3.14)

pow = np.abs(spec_fin1)
minpow = 5e-10
for i in range(len(tcenter_fin)):
      goo = np.where(pow[:,i] < minpow)[0]
      phase_fin[goo,i] = np.nan







#Load Marilia HF data
timesHF, freqsHF, powerHF = gld.load_marilia_HF(pld)




#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

#Gyrofrequencies (Hz)
fce = 28 * Bo   
#Interpolate high-cadence fce to low-cadence fpe
fce = np.interp(tspec, Bot, fce)
fcH = fce / 1836
fcO = fcH / 16



#-----------------------------------------
#create sonogram of Langmuir Probe data
#-----------------------------------------

#Load SLP data
lp = gld.load_slp(pld)
lptimes = lp[0]
lpvals = lp[1]
#plt.plot(lptimes, lpvals)

nfft=1024
fspecLP, tspecLP, powerLP, fsLP = fftspec.fft_spectrum_piecewise(lptimes, lpvals, fs_thres=0.1, nfft=nfft, noverlap=8)


#----------------------------------------------------
#----------------------------------------------------
#-----PLOTS-----
#----------------------------------------------------
#----------------------------------------------------

#-------------------------------------------------------------
#Overview plot for determining alt-dependent density changes
#-------------------------------------------------------------

xr = [100,520]
yr = [100,50000]
yrHF = [freqsHF.min(), freqsHF.max()]
#vr = [-72,-60] #strong power only
vrHF = [-95,-90]
vr = [-90,-50]
if pld == '380':
      vr_lp = [-250,-150] #380
elif pld == '381':
      vr_lp = [-230,-220] #381



fig, axs = plt.subplots(5,1, figsize=(10,9))
plt.title('GIRAFF ' + pld + ' - gir_analysis_plot_spectra.py')
ps.plot_spectrogram(timesHF, freqsHF, powerHF, vr=vrHF, yscale='linear', xr=xr, yr=yrHF, ylabel="power spectrum HF\nfreq(Hz); dB of (mV/m)^2/Hz", ax=axs[0],colorbar=0)
ps.plot_spectrogram(tspecLP, fspecLP, powerLP, vr=vr_lp, yscale='log', yr=yr, xr=xr, ylabel="Langmuir Probe power\nfreq(Hz); dB", ax=axs[1],colorbar=0)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='log',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[2],colorbar=0)
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='log',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[3],colorbar=0)
axs[4].plot(pvstBeam2['times'],8980*np.sqrt(pvstBeam2['pfisrdat_interp']/1e6), label='Beam 2 (VERTICAL)', color=colors[1], linewidth=4)
axs[4].plot(pvstBeam1['times'],8980*np.sqrt(pvstBeam1['pfisrdat_interp']/1e6), label='Beam 1', color=colors[0], linestyle='-.',linewidth=2)
axs[4].set_xlim(xr)
axs[4].set_yscale('linear')
axs[4].set_ylim(1e6,5e6)
axs[4].set_ylabel('fpe (Hz)')

minvert = 0
for i in range(5):
      axs[i].tick_params(axis='both', labelsize=6)
      axs[i].xaxis.label.set_fontsize(6)
      axs[i].yaxis.label.set_fontsize(6)
      axs[i].scatter(vertices[minvert:,0], vertices[minvert:,1],color='white',s=12)
      axs[i].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
      axs[i].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
      axs[i].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
      axs[i].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)

fig.tight_layout(pad=.8)
fig.subplots_adjust(left=0.15)

minvert = 0
for i in range(5):
      ax2 = axs[i].twinx()
      ax2.plot(time_ephem, alts, color='blue', linewidth=2)
      ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
      ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
      ax2.yaxis.set_label_position("left")
      ax2.tick_params(axis='y', labelleft=True, labelright=False)
      ax2.spines['left'].set_position(('outward', 60))







#------------------------
#Plots of sonograms, phase, and coherence
#------------------------

#xr = [100,520]
#xr = [128,188]
#xr = [188,248]
#xr = [248,308]
#xr = [308,368]
#xr = [368,428]
xr = [350,500]
#xr = [428,488]
#xr = [488,548]
#vr = [-72,-60] #strong power only
vr = [-90,-50]
#minzval = -80
#fig, axs = plt.subplots(2, figsize=(9,9), gridspec_kw={'height_ratios':[2,1]})
fig, axs = plt.subplots(3,2, figsize=(16,9))
plt.title('GIRAFF ' + pld)
yr = [0,10000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[0,0])
axs[0,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[1,0])
axs[1,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,0])
yr = [300,10000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='log',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[0,1])
axs[0,1].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='log',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[1,1])
axs[1,1].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin,vr=[0.6,1],zscale='linear',yscale='log',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,1])
fig.tight_layout(pad=0)

minvert = 0
for i in range(3):
      for j in range(2):
            axs[i,j].tick_params(axis='both', labelsize=6)
            axs[i,j].xaxis.label.set_fontsize(6)
            axs[i,j].yaxis.label.set_fontsize(6)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='white',s=12)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
            axs[i,j].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)
            if plot_harmonics:
                  axs[i,j].plot(tspec, 2*fcH, color='magenta',linewidth=1.5, label='2*fcH+')
                  axs[i,j].plot(tspec, 3*fcH, color='magenta',linewidth=1.5, label='3*fcH+')
                  axs[i,j].plot(tspec, 4*fcH, color='magenta',linewidth=1.5, label='4*fcH+')
                  axs[i,j].plot(tspec, 5*fcH, color='magenta',linewidth=1.5, label='5*fcH+')
                  axs[i,j].plot(tspec, 6*fcH, color='magenta',linewidth=1.5, label='6*fcH+')

fig.tight_layout(pad=2.5)
fig.subplots_adjust(left=0.15)

minvert = 0
for i in range(3):
      for j in range(2):
            ax2 = axs[i,j].twinx()
            ax2.plot(time_ephem, alts, color='blue', linewidth=2)
            ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
            ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
            ax2.yaxis.set_label_position("left")
            ax2.tick_params(axis='y', labelleft=True, labelright=False)
            ax2.spines['left'].set_position(('outward', 60))





print('h')




#---------------------
#Plot of different sonogram freq bands for comparison
#---------------------

#xr = [250,290] #VLF blob with DC structures
#xr = [290,370]  #Banded ELF w/ top band terminating just under fcH+
t0z = 200
t1z = t0z+10

xr = [t0z,t1z]
vr = [-90,-50]
fig, axs = plt.subplots(3, figsize=(8,9))
plt.title('GIRAFF ' + pld)
yr = [2000,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[0])
axs[0].set_xticklabels([])
yr = [500,2000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[1])
axs[1].set_xticklabels([])
yr = [10,500]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz); dB of (mV/m)^2/Hz",ax=axs[2])
fig.tight_layout(pad=0.01)

minvert = 0
for i in range(3):
      axs[i].scatter(vertices[0][minvert:,0], vertices[0][minvert:,1],color='white',s=12)
      axs[i].scatter(vertices[0][minvert:,0], vertices[0][minvert:,1],color='magenta',s=3, label='fLH')
      axs[i].plot(tspec, fce, color='black',linewidth=1.5, label='fce')
      axs[i].plot(tspec, fcH, color='magenta',linewidth=1.5, label='fcH+')
      axs[i].plot(tspec, fcO, color='blue', linewidth=1.5, label='fcO+')
      if plot_harmonics:
            axs[i].plot(tspec, 2*fcH, color='magenta',linewidth=1.5, label='2*fcH+')
            axs[i].plot(tspec, 3*fcH, color='magenta',linewidth=1.5, label='3*fcH+')
            axs[i].plot(tspec, 4*fcH, color='magenta',linewidth=1.5, label='4*fcH+')
            axs[i].plot(tspec, 5*fcH, color='magenta',linewidth=1.5, label='5*fcH+')
            axs[i].plot(tspec, 6*fcH, color='magenta',linewidth=1.5, label='6*fcH+')

      axs[i].legend()

fig.tight_layout(pad=2.5)
fig.subplots_adjust(left=0.25)

minvert = 0
for i in range(3):
      ax2 = axs[i].twinx()
      ax2.plot(time_ephem, alts, color='blue', linewidth=2)
      ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
      ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
      ax2.yaxis.set_label_position("left")
      ax2.tick_params(axis='y', labelleft=True, labelright=False)
      ax2.spines['left'].set_position(('outward', 60))







#plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)
