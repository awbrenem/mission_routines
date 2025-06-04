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
import plasma_params_get_flhr_freq as dflh
import pyIGRF
import pickle
#import scipy.io as sio
from scipy.io import readsav
import fft_spectrum_piecewise as fftspec
from scipy.interpolate import interp1d

#load the crossover and cutoff freq estimation from IRI (assuming H+ and O+ plasma)
#crossover = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/crossover_cutoff_freq_iri.pkl'
#tco,fcr,fco = pickle.load(open(crossover,'rb'))





##load SLP fixed bias density for the plot at end of mission
#goo = end_data_loader.load_slp_fixedbias()
#slpfb_dens = goo['Normalized Density [\m3]']
#slpfb_time = goo['ToF [s]']

"""
#--------------------------------------------------------------------------------------------------
#Read in polarization data of Bernstein waves from IDL save file
idldat = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/polarization_from_idl/pol_from_idl_fft=4096.sav')

times_pow = idldat['times_pow']
freqs_pow = idldat['freqs_pow']
data_pow = idldat['data_pow']

times_pol = idldat['times_pol']
freqs_pol = idldat['freqs_pol']
data_pol = idldat['data_pol']

times_elip = idldat['times_elip']
freqs_elip = idldat['freqs_elip']
data_elip = idldat['data_elip']

times_hel = idldat['times_hel']
freqs_hel = idldat['freqs_hel']
data_hel = idldat['data_hel']


#reduce ellipticity values to only those with a high deg polarization
data_elip2 = np.copy(data_elip)
for i in range(np.shape(data_pol)[0]):
      goodv = np.where(data_pol[i,:] < 0.7)[0]
      data_elip2[i,goodv] = np.nan

"""



#--------------------------------------------------------------------------------------------------


##load the by-eye lower hybrid values
#flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
#verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))



#Load E-fields data
v12 = GFL('381','VLF12D')
wf12, tdat = v12.load_data()
v34 = GFL('381','VLF34D')
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





##plot height spectra if IRI to see where my assumption of O+ and H+ is OK. 
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
#...NOTE: the onboard magnetometer gets wonky after the flip maneuver, likely due to its proximity to sc body. 
#Better to use IGRF data

"""
magv = GFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]
"""


"""
#-------------------
#composition
#98% O+ at 270 km (comparison of flh to SLP densities)
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

"""
#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

##Gyrofrequencies (Hz)
#fce = 28 * Bo   
##Interpolate high-cadence fce to low-cadence fpe
#fce = np.interp(slp_times, Bot, fce)
#fcH = fce / 1836
#fcO = fcH / 16


fceIGRF = 28 * BoIGRF
fcHIGRF = fceIGRF / 1836
fcOIGRF = fcHIGRF / 16


#fceIGRF = np.interp(times_interp, plot_times,fceIGRF)
#fcHIGRF = np.interp(times_interp, plot_times,fcHIGRF)
#fcOIGRF = np.interp(times_interp, plot_times,fcOIGRF)


#plasma freqs (Hz)
fpe = 8980 * np.sqrt(ne) 
#fpH = 104219
#fpO = 26054


#beta 
#beta_e = fpe / fce
"""


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




#load the by-eye lower hybrid determination
#flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
#verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))




#Load SLP density spectrogram 
#tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()










#---------------------------------
#Overall spectral plot (Fig 1)
#---------------------------------

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=8)
fspec, tspec, powerc34, fs = fftspec.fft_spectrum_piecewise(tdat34, wf34, fs_thres=0.1, nfft=nfft, noverlap=8)

fspec = fspec[:,0]


nfft = 4096
freqs_fin, freqs_coh, tcenter_fin, csd_fin, coh_fin, phase_fin, spec_fin1, spec_fin2, fs_fin =  ca.csd_spectrum_piecewise(tdat[tdat > 100], wf34[tdat > 100], wf12[tdat > 100], nfft=nfft, noverlap=8, fs_thres=0.2)
freqs_fin = freqs_fin[:,0]
freqs_coh = freqs_coh[:,0]
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



xr = [100,550]
#xr = [500,510]
#yr = [100,40000]
yr = [0,50000]
#vr = [-72,-60] #strong power only
vr = [-90,-50]
#minzval = -80
#fig, axs = plt.subplots(2, figsize=(9,9), gridspec_kw={'height_ratios':[2,1]})
fig, axs = plt.subplots(4, figsize=(9,9))
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
#ps.plot_spectrogram(sono12['atimesfft'],sono12['afreq'],sono12['afftpow'],vr=vr,zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])#,minzval=minzval)
axs[0].set_xticklabels([])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=vr,yscale='log',yr=[300,100000],xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1])
axs[1].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_fin,vr=[-120,120],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[2])
axs[2].set_xticklabels([])
ps.plot_spectrogram(tcenter_fin,freqs_coh,coh_fin,vr=[0.6,1],zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="coh (VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[3])
fig.tight_layout(pad=0)

minvert = 0
#axs[0].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='white',s=12)
#axs[0].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='magenta',s=3)
#for j in range(6,11):
#      axs[0].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)
#axs[1].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='white',s=12)
#axs[1].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='magenta',s=3)
#for j in range(6,15):
#      axs[1].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)



#plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)



#-----------------------------------------------------
#Low alt on downleg to show the connection b/t Bernstein and two-stream

#These are plotted separately b/c the log scale doesn't work with multi panel plots
#-----------------------------------------------------


fig,axs = plt.subplots(5, figsize=(12,9))
xr = [850,900]
yr = [0,9000]
#yr = [200,10000]
p0 = axs[0].scatter(ephem[0]['Time'], alt, cmap='turbo',s=1)
axs[0].grid()
axs[0].set_ylabel('Altitude (km)')
axs[0].set_xticklabels([])
axs[0].set_xlim(xr)
axs[0].set_ylim(70,240)
cbar = plt.colorbar(p0)

ax2 = axs[0].twinx()

p1 = ax2.scatter(slpfb_time,slpfb_dens,s=1,color='magenta')
ax2.set_ylim(0,2e11)
ax2.set_ylabel('SLP density (m-3)')

ps.plot_spectrogram(tspec,fspec12,np.abs(powerc12),vr=[-75,-65],yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1])
axs[1].set_xticklabels([])
ps.plot_spectrogram(tchunks,ff,coh,zscale='linear',vr=[0.4,1],ax=axs[2],xr=xr,yr=yr,yscale='linear',ylabel='Coherence\n(VLF12, VLF34)\nfreq(Hz)')
axs[2].set_xticklabels([])
ps.plot_spectrogram(times_pol,freqs_pol,data_pol,zscale='linear',vr=[0.7,1],xr=xr,yr=yr,yscale='linear',ax=axs[3],ylabel='deg polarization\nfreq(Hz)')
axs[3].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,zscale='linear',vr=[-0.8,0.8],xr=xr,yr=yr,yscale='linear',ax=axs[4],ylabel='ellipticity\nfreq(Hz)',xlabel='time(sec)')
fig.tight_layout(pad=0)


#axs[3].set_xticklabels([])
#ps.plot_spectrogram(tchunks,ff,phase,zscale='linear',vr=[-3.14,3.14],ax=axs[4],xr=xr,yr=yr,yscale='log')
#fig.tight_layout(pad=1)
#plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)

for i in range(5): 
      if i > 0:
            #axs[i].plot(tco,fcr,color='magenta')
            if xr[0] < 600:
                  if i == 1 or i == 3:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='white',s=22)
                  else: 
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='black',s=22)
                  axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='magenta',s=10)
            if xr[0] > 600:
                  if i == 1 or i == 3:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='white',s=22)
                  else:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='black',s=22)
                  axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='magenta',s=10)
            #axs[i].plot(plot_times, fcOIGRF, color='white')
            #for j in range(11):
            #      if i == 1 or i == 3:
            #            axs[i].plot(plot_times,j*fcHIGRF,color='white',linestyle='dotted',linewidth=0.9)
            #      else:
            #            axs[i].plot(plot_times,j*fcHIGRF,color='black',linestyle='dotted',linewidth=0.9)

            axs[i].set_ylim(yr)

      axs[i].set_xlim(xr)

plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)


#plt.plot(alts,flhr_z)





#------------------------------------------------------------------
#Full timerange (removed density plot and adjusted zlimits)
#------------------------------------------------------------------


xr = [100,900]
yr = [2000,9000]
pk = {'cmap':'turbo'}


fig,axs = plt.subplots(4, figsize=(8,9))
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-75,-60],yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0],plot_kwargs=pk)
axs[0].set_xticklabels([])
ps.plot_spectrogram(tdens,fdensH,np.transpose(sdensH),vr=[-14,-11],zscale='linear',yr=yr,xr=xr,ylabel='Density power spectrum\nfreq(Hz)\nrelative units',ax=axs[1],plot_kwargs=pk)
axs[1].set_xticklabels([])
ps.plot_spectrogram(times_pol,freqs_pol,data_pol,zscale='linear',vr=[0.7,1],xr=xr,yr=yr,ax=axs[2],ylabel='degree polarization\nfreq(Hz)',plot_kwargs=pk)
axs[2].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,zscale='linear',vr=[-0.5,1],xr=xr,yr=yr,ax=axs[3],xlabel='time(sec)',ylabel='ellipticity\nfreq(Hz)',plot_kwargs=pk)
fig.tight_layout(pad=0)


for i in range(4): 
      axs[i].plot(tco,fcr,color='magenta')
      if xr[0] < 600:
            axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='white',s=12)
            axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='magenta',s=3)
      if xr[0] > 600:
            axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='white',s=12)
            axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='magenta',s=3)
      axs[i].plot(plot_times, fcOIGRF, color='white')
      for j in range(11):
           if i < 3:
                 axs[i].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)
           else:
                 axs[i].plot(plot_times,j*fcHIGRF,color='black',linestyle='--',linewidth=0.9)
      axs[i].set_xlim(xr)
      axs[i].set_ylim(yr)



plt.savefig("/Users/abrenema/Desktop/tst1.pdf", dpi=350)










#--------------------------------
#Separate upleg and downleg plots
#------------------------------------------------------------------

#pk = {'cmap':'viridis'}
#pk = {'cmap':'cividis'}
pk = {'cmap':'turbo'}

#xr = [100,240]
xr = [650,900]
#xr = [450,600]
#xr = [80,240]
yr = [4000,9000]
fig,axs = plt.subplots(6, figsize=(8,9))
p1 = axs[0].scatter(ephem[0]['Time'], alt, cmap='turbo',s=1)
axs[0].grid()
axs[0].set_ylabel('Altitude(km)')
axs[0].set_xticklabels([])
axs[0].set_xlim(xr)
cbar = plt.colorbar(p1)
#cbar.set_label('freq(Hz)',rotation=270,labelpad=15)

#axs[0].set_position([0.12,0.83,0.705,0.84])
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-75,-60],yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0],plot_kwargs=pk)
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-75,-60],yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1],plot_kwargs=pk)
ps.plot_spectrogram(tspec,fspec,powerc12,vr=[-75,-60],yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[1],plot_kwargs=pk)
axs[1].set_xticklabels([])
ps.plot_spectrogram(tchunks,ff,coh,zscale='linear',vr=[0.4,1],ax=axs[2],xr=xr,yr=yr,yscale='linear',ylabel='Coherence\n(VLF12, VLF34)\nfreq(Hz)')
axs[2].set_xticklabels([])
ps.plot_spectrogram(times_pol,freqs_pol,data_pol,zscale='linear',vr=[0.7,1],xr=xr,yr=yr,ax=axs[3],ylabel='degree polarization\nfreq(Hz)',plot_kwargs=pk)
axs[3].set_xticklabels([])
ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,zscale='linear',vr=[-1,1],xr=xr,yr=yr,ax=axs[4],ylabel='ellipticity\nfreq(Hz)',plot_kwargs=pk)
axs[4].set_xticklabels([])
ps.plot_spectrogram(tdens,fdensH,np.transpose(sdensH),vr=[-12,-10],zscale='linear',yr=yr,xr=xr,xlabel='time(sec)',ylabel='Density power spectrum\nfreq(Hz)\nrelative units',ax=axs[5],plot_kwargs=pk)

fig.tight_layout(pad=0)
#ps.plot_spectrogram(times_elip,freqs_elip,data_elip2,zscale='linear',vr=[0.2,0.8],xr=xr,yr=yr,ax=axs[4],xlabel='time(sec)',ylabel='ellipticity\nfreq(Hz)',plot_kwargs=pk)
#fig.tight_layout(pad=0)

for i in range(6): 
      if i > 0:
            axs[i].plot(tco,fcr,color='magenta')
            if xr[0] < 600:
                  if i == 1 or i == 3 or i == 5:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='white',s=12)
                  else: 
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='black',s=12)
                  axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='magenta',s=3)
            if xr[0] > 600:
                  if i == 1 or i == 3 or i == 5:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='white',s=12)
                  else:
                        axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='black',s=12)
                  axs[i].scatter(verticesLOW[:,0], verticesLOW[:,1],color='magenta',s=3)
            axs[i].plot(plot_times, fcOIGRF, color='white')
            for j in range(11):
                  if i == 1 or i == 3 or i == 5:
                        axs[i].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)
                  else:
                        axs[i].plot(plot_times,j*fcHIGRF,color='black',linestyle='--',linewidth=0.9)

            axs[i].set_ylim(yr)

      axs[i].set_xlim(xr)

plt.savefig("/Users/abrenema/Desktop/tst1.pdf", dpi=350)






#---------------------------------
#Spectra focusing on harmonics within the VLF hiss
#---------------------------------

#special run with cohmin = 0
timechunk = 0.2
cohz, phasez, tchunksz, ffz = ca.cross_spectral_density_spectrogram(wf34[tdat > 100],wf12[tdat > 100],tdat[tdat > 100],30000,timechunk,nperseg=32768,plot=False,coh_min=0)
phasez = phasez*(180/3.14)

xr = [600,850]
yr = [4000,10000]
vr = [-74,-63]
fig, axs = plt.subplots(2, figsize=(9,9), gridspec_kw={'height_ratios':[2,1]})
ps.plot_spectrogram(sono12['atimesfft'],sono12['afreq'],sono12['afftpow'],vr=vr,zscale='linear',yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])#,minzval=minzval)
axs[0].set_xticklabels([])
ps.plot_spectrogram(tchunksz,ffz,phasez,vr=[-120,120],zscale='linear',yscale='linear',yr=[3000,9000],xr=xr,ylabel="phase (deg; VLF12,VLF34)\nfreq(Hz)",xlabel='time since launch (sec)',ax=axs[1])
fig.tight_layout(pad=0)

minvert = 0
axs[0].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='white',s=12)
axs[0].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='magenta',s=3)
for j in range(6,12):
      axs[0].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)
axs[1].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='white',s=12)
axs[1].scatter(verticesLOW[minvert:,0], verticesLOW[minvert:,1],color='magenta',s=3)
for j in range(6,15):
      axs[1].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)





plt.savefig("/Users/abrenema/Desktop/tst3.pdf", dpi=350)


















#verticesLOW[9:11,1] += 400




#ps.plot_spectrogram(tdens,fdensH,np.transpose(sdensH),vr=[-20,-5],zscale='linear',yr=[1400,15000],xr=xr,ylabel='Density power spectru  x\nfreq(Hz)',xlabel='time(sec)')
#ps.plot_spectrogram(tdens,fdensH,np.transpose(sdensH),vr=[-12,-10],zscale='linear',yr=yr,xr=xr,ylabel='Density power spectrum\nfreq(Hz)',xlabel='time(sec)')
#plt.savefig("/Users/abrenema/Desktop/tst4.pdf", dpi=350)



#Look at low freq SLP density sonogram
yr=[1,2000]
xr=[100,240]
ps.plot_spectrogram(tdens,fdensL,np.transpose(sdensL),vr=[-12,-5],yscale='linear',zscale='linear',yr=yr,xr=xr,ylabel='Density power spectrum\nfreq(Hz)',xlabel='time(sec)')















xr = [200,600]
yr = [4000,12000]
ps.plot_spectrogram(tspec,fspec,powerc_v2m,vr=[-100,-90],yr=yr,xr=xr,yscale='linear',ylabel="power spectrum wf34\nfreq(Hz)\ndB (mV/m)^2/Hz")


fig,axs = plt.subplots(2, figsize=(8,9))
ps.plot_spectrogram(tspec,fspec,powerc_v2m,vr=[-105,-95],yr=yr,xr=xr,yscale='linear',ylabel="power spectrum wf34\nfreq(Hz)\ndB (mV/m)^2/Hz",ax=axs[0])
ps.plot_spectrogram(tspec,fspec,powerc_v2m,vr=[-95,-90],yr=yr,xr=xr,yscale='linear',ylabel="power spectrum wf34\nfreq(Hz)\ndB (mV/m)^2/Hz",ax=axs[1])
axs[0].set_xticklabels([])

axs[0].plot(plot_times[flh_gootimes[0]], flh[flh_gootimes[0]],color='white')
axs[0].plot(plot_times[flh_gootimes[0]], flh2[flh_gootimes[0]], color='white')
axs[0].plot(plot_times, fcOIGRF, color='white')
for j in range(11):
      axs[0].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)








#Test of SLP sonogram time tagging

xr = [100,900]
yr = [2000,9000]
pk = {'cmap':'turbo'}


fig,axs = plt.subplots(3, figsize=(8,9))
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-75,-60],yr=yr,xr=xr,ylabel="power spectrum VLF34\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0],plot_kwargs=pk)
axs[0].set_xticklabels([])
ps.plot_spectrogram(tdens+20,fdensH,np.transpose(sdensH),vr=[-16,-11],zscale='linear',yr=yr,xr=xr,ylabel='Density power spectrum\nfreq(Hz)\nTIMESHIFT',ax=axs[1],plot_kwargs=pk)
axs[1].set_xticklabels([])
ps.plot_spectrogram(tdens,fdensH,np.transpose(sdensH),vr=[-16,-11],zscale='linear',yr=yr,xr=xr,ylabel='Density power spectrum\nfreq(Hz)\nNO TIMESHIFT',ax=axs[2],plot_kwargs=pk)
axs[2].set_xticklabels([])






for i in range(4): 
      axs[i].plot(tco,fcr,color='magenta')
      if xr[0] < 600:
            axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='white',s=12)
            axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='magenta',s=3)
      if xr[0] > 600:
            axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='white',s=12)
            axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='magenta',s=3)
      axs[i].plot(plot_times, fcOIGRF, color='white')
      for j in range(11):
           if i < 3:
                 axs[i].plot(plot_times,j*fcHIGRF,color='white',linestyle='--',linewidth=0.9)
           else:
                 axs[i].plot(plot_times,j*fcHIGRF,color='black',linestyle='--',linewidth=0.9)
      axs[i].set_xlim(xr)
      axs[i].set_ylim(yr)

