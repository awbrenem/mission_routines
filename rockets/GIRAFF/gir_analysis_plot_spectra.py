#Make nice looking spectral plots for GIRAFF

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
#import end_data_loader
from scipy import signal
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
from scipy.interpolate import interp1d

#load the crossover and cutoff freq estimation from IRI (assuming H+ and O+ plasma)
#crossover = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/crossover_cutoff_freq_iri.pkl'
#tco,fcr,fco = pickle.load(open(crossover,'rb'))



pld = '381'



##load SLP fixed bias density for the plot at end of mission
#goo = end_data_loader.load_slp_fixedbias()
#slpfb_dens = goo['Normalized Density [\m3]']
#slpfb_time = goo['ToF [s]']


#--------------------------------------------------------------------------------------------------
#Read in polarization data of GIRAFF waves from IDL save file
pathz = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/polarization_from_idl/'

idldat1 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=4096_100-300sec.sav')
idldat2 = readsav(pathz + 'Giraff_'+pld+'_pol_from_idl_fft=4096_300-550sec.sav')

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


#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_381_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))[0]



#Load E-fields data
v12 = GFL(pld,'VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL(pld,'VLF34D')
wf34, tdat34 = v34.load_data()


##There are some artificial signals in the Bernstein waves at higher freqs that are picked up by 12 but not 34. 
##Let's find where these are by finding high values of the below ratio. Then, we'll eliminate polarization data at these
##high values
#powerrat = powerc12 / powerc
#powerrat_lim = 10.



#ephem = end_data_loader.load_ephemeris()
#alt = ephem[0]['Altitude'] 

#iri = end_data_loader.load_iri()
#ig = end_data_loader.load_ig()
#slp = end_data_loader.load_slp()
#tl = end_data_loader.load_timeline()

#slp_times = slp['ToF [s]']
#slp_alt = slp['Alt [km]']





##plot height spectra of IRI to see where my assumption of O+ and H+ is OK. 
#plt.plot(iri['Height(km)'], iri['O_ions'])
#plt.plot(iri['Height(km)'], iri['H_ions'])
#plt.plot(iri['Height(km)'], iri['He_ions'])
#plt.plot(iri['Height(km)'], iri['O2_ions'])
#plt.plot(iri['Height(km)'], iri['NO_ions'])
#plt.plot(iri['Height(km)'], iri['N_ions'])



#Use this as the plotted time cadence
#plot_times = np.asarray(ephem[0]['Time'])
#plot_alt = np.asarray(ephem[0]['Altitude'])
#plot_times = plot_times[0::20]
#plot_alt = plot_alt[0::20]


#glat = 78 + (55/60)
#glon = 11 + (55/60)
#BoIGRF = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])





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




"""
#science collection times
sstart = tl[1]
send = tl[2]
bsstart = tl[3]
bsend = tl[4]
"""

#magnetic field 
magv = GFL(pld,'mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]


"""
#-------------------
#composition
# (comparison of flh to SLP densities)
#-------------------
fOp2 = 0.98
fHp2 = 0.02

fOp = iri['O_ions']/100
fHp = iri['H_ions']/100

#Interpolate IRI alt to SLP alt
irialt = iri['Height(km)']
fOp = np.interp(plot_alt, irialt, fOp)
fHp = np.interp(plot_alt, irialt, fHp)

#plt.plot(plot_times,fOp)
"""

"""
#------------------
#Densities
#------------------

nneut = ig['Dens_neutral']
ni = slp['SLP Ni [/m3]'] / 100**3
ni = np.asarray(ni)
ne = ni #cm-3

times_interp = slp_times

ne = np.interp(plot_times, times_interp, ne)
ni = np.interp(plot_times, times_interp, ni)

nOp = fOp * ne 
nHp = fHp * ne

"""





#Load SLP density spectrogram 
#tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()




#---------------------------------
#Overall spectral plot (Fig 1)
#---------------------------------

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=8)
fspec, tspec, powerc34, fs = fftspec.fft_spectrum_piecewise(tdat34, wf34, fs_thres=0.1, nfft=nfft, noverlap=8)


nfft = 4096
freqs_fin, tcenter_fin, csd_fin, coh_fin, phase_fin, fs_fin, spec_fin1, spec_fin2 =  ca.csd_spectrum_piecewise(tdat[tdat > 100], wf34[tdat > 100], wf12[tdat > 100], nfft=nfft, noverlap=8, fs_thres=0.2)


#freqs_fin = freqs_fin[:,0]
#freqs_coh = freqs_coh[:,0]
phase_fin = phase_fin*(180/3.14)

pow = np.abs(spec_fin1)
minpow = 5e-10
for i in range(len(tcenter_fin)):
      goo = np.where(pow[:,i] < minpow)[0]
      phase_fin[goo,i] = np.nan



#special run with cohmin = 0
#timechunk = 1.0
#cohz, phasez, tchunksz, ffz = ca.cross_spectral_density_spectrogram(wf34[tdat > 100],wf12[tdat > 100],tdat[tdat > 100],799999,timechunk,nperseg=2048,plot=False,coh_min=0)
#phasez = phasez*(180/3.14)



#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

#Gyrofrequencies (Hz)
fce = 28 * Bo   
#Interpolate high-cadence fce to low-cadence fpe
fce = np.interp(tspec, Bot, fce)
fcH = fce / 1836
fcO = fcH / 16


#fceIGRF = 28 * BoIGRF
#fcHIGRF = fceIGRF / 1836
#fcOIGRF = fcHIGRF / 16

#fceIGRF = np.interp(times_interp, plot_times,fceIGRF)
#fcHIGRF = np.interp(times_interp, plot_times,fcHIGRF)
#fcOIGRF = np.interp(times_interp, plot_times,fcOIGRF)


#plasma freqs (Hz)
#fpe = 8980 * np.sqrt(ne) 
#fpH = 104219
#fpO = 26054

#beta 
#beta_e = fpe / fce

#x = np.linspace(0, 2*np.pi, 10)
#y = np.sin(x)
#xvals = np.linspace(0, 2*np.pi, 50)
#yinterp = np.interp(xvals, x, y)


##Lower hybrid freq calculated based on fractional composition (f >> fH)
#flh = dflh.flhr_IonMassFractions(ne, fceIGRF, fHp, fOp)
#flh2 = dflh.flhr_IonMassFractions(ne, fceIGRF, fHp2, fOp2)



#Create a second version of flh that has times only up to 825 sec.
#The accuracy of my flh determination decreases after this time and the overplotted flh line on the
#zoomed in Bernstein plots looks bad. 

#flh_gootimes = np.where((plot_times > 150) & (plot_times < 825))





#------------------------
#Plot of sonograms, phase, and coherence
#------------------------


xr = [100,550]
#vr = [-72,-60] #strong power only
vr = [-90,-50]
#minzval = -80
#fig, axs = plt.subplots(2, figsize=(9,9), gridspec_kw={'height_ratios':[2,1]})
fig, axs = plt.subplots(3,2, figsize=(16,9))
plt.title('GIRAFF ' + pld)
yr = [0,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,0])
axs[0,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[1,0])
axs[1,0].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,0])
yr = [300,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='log',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0,1])
axs[0,1].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='log',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[1,1])
axs[1,1].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin,vr=[0.6,1],zscale='linear',yscale='log',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2,1])
fig.tight_layout(pad=0)

minvert = 0
for i in range(3):
      for j in range(2):
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='white',s=12)
            axs[i,j].scatter(vertices[minvert:,0], vertices[minvert:,1],color='magenta',s=3)
            axs[i,j].plot(tspec, fce, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcH, color='white',linestyle='--',linewidth=0.9)
            axs[i,j].plot(tspec, fcO, color='white',linestyle='--',linewidth=0.9)



print('h')




#---------------------
#Plot of different sonogram freq bands for comparison
#---------------------

#xr = [250,290] #VLF blob with DC structures
#xr = [290,370]  #Banded ELF w/ top band terminating just under fcH+

xr = [100,500]
vr = [-90,-50]
fig, axs = plt.subplots(3, figsize=(16,9))
plt.title('GIRAFF ' + pld)
yr = [2000,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
axs[0].set_xticklabels([])
yr = [500,2000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1])
axs[1].set_xticklabels([])
yr = [10,500]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[2])
fig.tight_layout(pad=0.01)

minvert = 0
for i in range(3):
      axs[i].scatter(vertices[0][minvert:,0], vertices[0][minvert:,1],color='white',s=12)
      axs[i].scatter(vertices[0][minvert:,0], vertices[0][minvert:,1],color='magenta',s=3, label='fLH')
      axs[i].plot(tspec, fce, color='black',linewidth=1.5, label='fce')
      axs[i].plot(tspec, fcH, color='magenta',linewidth=1.5, label='fcH+')
      axs[i].plot(tspec, fcO, color='blue', linewidth=1.5, label='fcO+')
      axs[i].legend()








#plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)
