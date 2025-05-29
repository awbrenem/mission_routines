#Compare modulation of VLF hiss to precipitating electrons (reflected from much larger altitudes by field-aligned potentials (HARPS; Glocer+24))
#These might be the free energy source. 

#Also compare the hiss to density modulation.

#Hiss is either:
#     -modulated due to density
#     -modulated due to structure of reflection potentials
#     -modulated due to both


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
from scipy import signal
import numpy as np 
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import filter_wave_frequency as filt
import plasma_params_get_flhr_freq as dflh
import pyIGRF
import pickle
#import scipy.io as sio
from scipy.io import readsav
import signal_smoothing


iri = end_data_loader.load_iri()
ig = end_data_loader.load_ig()
slp = end_data_loader.load_slp()
df, goodscienceS, goodscienceE, badscienceS, badscienceE = end_data_loader.load_timeline()
ephem1, ephem2 = end_data_loader.load_ephemeris()

plt.plot(ephem1['Time'],ephem1['ECEF-Z-VEL'])

slp_times = slp['ToF [s]']
slp_alt = slp['Alt [km]']



#-------------------------------------------------------------------
#Read in the PES data Glyn made for me with precipitating electrons
#-------------------------------------------------------------------

pesdat, tsec = end_data_loader.load_pes_electron_precipitation()


#tsec values are defined as the beginning of each timestep. Shift to middle
delta_tsec = tsec[1] - tsec[0]

tsec += delta_tsec


#------------------
#Densities
#------------------

times,freqslow,datalow,freqshig,datahig = end_data_loader.load_slp_sonogram()

#grab lowest few freqs for DC density
denstimes = times 

nfreqs = 3 #sum over this many frequencies
#densDC, powarr, tarr = ps.slice_spectrogram(0, freqslow, np.transpose(datalow), nsec=nfreqs)
densDC, powarr, tarr = ps.slice_spectrogram(0, freqslow, datalow, nsec=0) #nfreqs)

plt.plot(denstimes,densDC)


nneut = ig['Dens_neutral']
ni = slp['SLP Ni [/m3]'] / 100**3
ni = np.asarray(ni)
ne = ni #cm-3




#Read in VLF wave data
v34 = EFL('VLF34D')
fs = v34.chnspecs['fs']
wf34, tdat = v34.load_data_gainphase_corrected()

nps = 4096
fspec, tspec, powerc = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
#Returns mV/m**2/Hz



#load the by-eye lower hybrid determination
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))


#Interpolate the lower hybrid id to the time cadence of the spectral data
vl2 = np.interp(tspec,verticesLOW[:,0],verticesLOW[:,1])

#Remove PES data during maneuver times 
    #badscienceS
    #[0,   192, 272.5, 352.5, 432.5, 512.5, 592.5, 692.5, 772.5, 852.5, 900.6]
    #badscienceE
    #[125, 203, 283,   362.9, 442.8, 522.5, 620.5, 702.5, 782.8, 862.5, 1000]
    #tsec - 5
    #array([129, 139, 149, 159, 169, 179, 189, 199, 209, 219, 229, 239, 249,  12
    #       259, 269, 279, 289, 299, 309, 319, 329, 339, 349, 359, 369, 379,  25
    #       389, 399, 409, 419, 429, 439, 449, 459, 469, 479, 489, 499, 509,  38
    #       519, 529, 539, 549, 559, 569, 579, 589, 599, 609, 619, 629, 639,  51
    #       649, 659, 669, 679, 689, 699, 709, 719, 729, 739, 749, 759, 769,  64
    #       779, 789, 799, 809, 819, 829, 839, 849, 859, 869, 879, 889],
    #      dtype='timedelta64[s]')

#badtsec = [46,47] #flip maneuver only
#badtsec = [6,14,22,30,38,46,56,64,72] #first data point that overlaps w/ maneuver
badtsec = [6,7,14,15,22,23,30,31,38,38,46,47,56,57,64,65,72,73] #all data points overlapping with maneuver

pesdat['el_influx_total'][badtsec] = np.nan
pesdat['el_influx_below_55eV'][badtsec] = np.nan
pesdat['total_influx_above_reflection_potential'][badtsec] = np.nan
pesdat['polar_rain'][badtsec] = np.nan




#Remove spectral data below the lower hybrid frequency 
powerc = np.abs(powerc)

powerc2 = np.copy(powerc) #version with maximally removed data
powerc3 = np.copy(powerc) #version with only times during maneuvers removed
powerc4 = np.copy(powerc) #version only DC values remaining
fudgefactor = 1000 #Hz (make sure we're not getting any lower hybrid continuous power into the integration)
for i in range(len(tspec)):
    powerc2[fspec < vl2[i]+fudgefactor, i] = np.nan

#Remove spectral data above DC freqs
for i in range(len(tspec)):
    powerc4[(fspec < 20) | (fspec > 500), i] = np.nan


#remove spectral power during attitude maneuvers
for i in range(len(badscienceS)):
    cond = (tspec >= badscienceS[i]) & (tspec <= badscienceE[i])
    powerc2[:,cond] = np.nan
    powerc3[:,cond] = np.nan
    powerc4[:,cond] = np.nan




#remove power at beginning and end of mission
powerc2[:,tspec > 890] = np.nan
powerc4[:,tspec > 890] = np.nan


#Integrate hiss power for f > flhr
powersum = np.zeros(len(tspec))
for i in range(len(powersum)):
    powersum[i] = np.nansum(powerc2[:,i])

#Integrate DC power
powersumDC = np.zeros(len(tspec))
for i in range(len(powersumDC)):
    powersumDC[i] = np.nansum(powerc4[:,i])
powersumDC[tspec > 880] = np.nan


#Remove hiss and DC power during attitude maneuvers
for i in range(len(badscienceS)):
    cond = (tspec >= badscienceS[i]) & (tspec <= badscienceE[i])
    powersum[cond] = np.nan
    powersumDC[cond] = np.nan

#smooth integrated hiss power 
smooveb = signal_smoothing.sliding_box_avg(powersum, box_pts=80)
smooveb2 = signal_smoothing.savgol_filter(powersum,window=50,poly=3)

smoovebDC = signal_smoothing.sliding_box_avg(powersumDC, box_pts=80)
smooveb2DC = signal_smoothing.savgol_filter(powersumDC,window=50,poly=3)

smooveb[np.where(smooveb == 0)[0]] = np.nan
smooveb2[np.where(smooveb2 == 0)[0]] = np.nan
smoovebDC[np.where(smoovebDC == 0)[0]] = np.nan
smooveb2DC[np.where(smooveb2DC == 0)[0]] = np.nan
 


#Bin the hiss data to the cadence of the precipitating electron flux (10 sec)
#Try (1) fractional occurrence; (2) avg/med/max amplitudes
hisstot = np.zeros(len(tsec))
hissmed = np.zeros(len(tsec))
hissmax = np.zeros(len(tsec))
DCtot = np.zeros(len(tsec))
DCmed = np.zeros(len(tsec))
DCmax = np.zeros(len(tsec))

for i in range(len(tsec)):
    tgood = np.where((tspec >= tsec[i]/np.timedelta64(1,'s')) & (tspec < tsec[i]/np.timedelta64(1,'s')+10))[0]
    hisstot[i] = np.nansum(powersum[tgood])
    hissmed[i] = np.median(powersum[tgood])
    hissmax[i] = np.nanmax(powersum[tgood])
    DCtot[i] = np.nansum(powersumDC[tgood])
    DCmed[i] = np.median(powersumDC[tgood])
    DCmax[i] = np.nanmax(powersumDC[tgood])


#remove zero values
hisstot[hisstot == 0] = np.nan
hissmed[hissmed == 0] = np.nan
hissmax[hissmax == 0] = np.nan
DCtot[DCtot == 0] = np.nan
DCmed[DCmed == 0] = np.nan
DCmax[DCmax == 0] = np.nan

for i in range(len(badscienceS)):
    cond = (tsec/np.timedelta64(1,'s') >= badscienceS[i]) & (tsec/np.timedelta64(1,'s') <= badscienceE[i])
    hisstot[cond] = np.nan
    hissmed[cond] = np.nan
    hissmax[cond] = np.nan
    DCtot[cond] = np.nan
    DCmed[cond] = np.nan
    DCmax[cond] = np.nan
  


#normalize all data binned at cadence of PES electrons 
hisstot = hisstot / np.nanmax(hisstot)
hissmed = hissmed / np.nanmax(hissmed)
hissmax = hissmax / np.nanmax(hissmax)
DCtot = DCtot / np.nanmax(DCtot)
DCmed = DCmed / np.nanmax(DCmed)
DCmax = DCmax / np.nanmax(DCmax)
pesdat['el_influx_below_55eV'] = pesdat['el_influx_below_55eV'] / np.nanmax(pesdat['el_influx_below_55eV'])


#Load SLP density spectrogram 
tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()


vr = [-40,-30]
yr = [4000,15000]
ys = 'linear'
xr = [100,900]


#j=0
#xr = [goodscienceS[j],goodscienceE[j]]

fig, axs = plt.subplots(6, figsize=(14,8))
ps.plot_spectrogram(tspec,fspec,powerc,vr=vr,yr=yr,xr=xr,yscale=ys,xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[0],colorbar=0)
axs[0].plot(tspec,vl2)
ps.plot_spectrogram(tspec,fspec,powerc2,vr=vr,yr=yr,xr=xr,yscale=ys,xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[1],colorbar=0)
axs[0].plot(tspec,vl2+fudgefactor)
axs[2].plot(tspec,powersum)
axs[2].plot(tspec,smooveb)
#axs[2].plot(tspec,smooveb2)
axs[2].set_ylim(0,0.5)
axs[2].set_xlim(xr)
#axs[3].plot(tsec,pesdat['el_influx_total']) #total influx
axs[3].plot(tsec,pesdat['el_influx_below_55eV'],'.') #reflected photoelectron precipitation
#axs[3].plot(tsec,pesdat['total_influx_above_reflection_potential'])
#axs[3].plot(tsec,pesdat['polar_rain'])
axs[3].set_xlim(xr)
#axs[3].set_ylim(0,10)

ps.plot_spectrogram(tspec,fspec,powerc4,vr=vr,yr=[10,8000],xr=xr,yscale='log',xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[4],colorbar=0)
axs[5].plot(tspec,powersumDC)
axs[5].plot(tspec,smoovebDC)
axs[5].set_xlim(xr)

print('h')







fig, axs = plt.subplots(5, figsize=(14,8))
ps.plot_spectrogram(tspec,fspec,powerc2,vr=vr,yr=yr,xr=xr,yscale=ys,xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[0],colorbar=0)
axs[1].plot(tspec,powersum)
axs[1].plot(tspec,smooveb)
axs[1].set_ylim(0,0.5)
axs[2].plot(tspec,powersumDC)
axs[2].plot(tspec,smoovebDC)
axs[3].plot(tsec,hisstot,'.',tsec,pesdat['el_influx_below_55eV'])
axs[4].plot(tsec,hisstot,'.',tsec,DCtot) #reflected photoelectron precipitation
#axs[5].plot(tsec,hissmed,'.',tsec,pesdat['el_influx_below_55eV'])
#axs[6].plot(tsec,hissmed,'.',tsec,DCmed) #reflected photoelectron precipitation
#axs[7].plot(tsec,hissmax,'.',tsec,pesdat['el_influx_below_55eV'])
#axs[8].plot(tsec,hissmax,'.',tsec,DCmax) #reflected photoelectron precipitation
axs[1].set_ylabel('Hiss power')
axs[2].set_ylabel('DC power')
axs[3].set_ylabel('binned hiss\ntotal power(blue)\nvs reflected\ne- flux')
axs[4].set_ylabel('binned hiss\ntotal power(blue)\nvs binned DC\ntotal power')

for i in range(5):
    axs[i].set_xlim(xr)



#Plot comparisons b/t quantities

#all data
plt.plot(pesdat['el_influx_below_55eV'], hisstot,'.')
plt.plot(DCtot, hisstot,'.')


#upleg only 
tgood = np.where(tsec/np.timedelta64(1,'s') < 592.5)[0]

fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(2,1,1, adjustable='box',aspect=1)
ax2 = fig.add_subplot(2,1,2, adjustable='box',aspect=1)

ax1.plot(pesdat['el_influx_below_55eV'][tgood], hisstot[tgood],'.')
ax2.plot(DCtot[tgood], hisstot[tgood],'.')

ax1.set_xlim(0,1.1)
ax2.set_xlim(0,1.1)
ax1.set_ylim(0,1.1)
ax2.set_ylim(0,1.1)
ax1.set_xlabel('e- flux')
ax1.set_ylabel('hiss power')
ax2.set_xlabel('DC power')
ax2.set_ylabel('hiss power')













vr = [-40,-30]
yr = [4000,15000]
ys = 'linear'

xr = [180,220]
xr = [260,290]
xr = [350,370]
xr = [430,450]
xr = [450,550]
xr = [510,530]
xr = [580,640]
xr = [690,710]
xr = [850,870]
xr = [760,800]

fig,axs = plt.subplots(2, figsize=(14,9))
ps.plot_spectrogram(tspec,fspec,powerc,vr=vr,yr=[10,8000],xr=xr,yscale='log',xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[0],colorbar=0)
ps.plot_spectrogram(tspec,fspec,powerc3,vr=vr,yr=[10,8000],xr=xr,yscale='log',xlabel='time(sec)',
                    ylabel="power spectrum x'\nfreq(Hz)",ax=axs[1],colorbar=0)


#plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)
