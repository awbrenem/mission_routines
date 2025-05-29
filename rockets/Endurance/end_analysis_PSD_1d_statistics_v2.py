#Build up statistics of the 1D power spectral density vs time for various plots. 
#Note that it takes many minutes to create the slice files, so I've saved data as a pickle file. 

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import numpy as np 
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import pyIGRF
from scipy.io import readsav
import copy 
import pickle
from copy import deepcopy


#-------------------------------------------------------------
#-------------------------------------------------------------
#Load a pickle file with the data or create new file with this program (takes a while)
#load_pickle = True

#leg = 'upleg'
leg = 'downleg'

#Define a minimum threshold power for consideration
if leg == 'upleg':
    #powthres = 1.5e-13
    powthres = 3.6e-14
    tr = [115,220]   #timerange for plotting
    fr = [5400,7500] #freq range for plotting
    er = [-0.8,0.3]  #ellipticity range for plotting
    fcolorrange = [5500,7400] #for squeezing the colorbar
if leg == 'downleg':
    #powthres = 1.5e-12
    powthres = 2e-13
    tr = [710,850]
    fr = [5000,7500]
    er = [-0.8,0.6]
    fcolorrange = [5000,7400] #for squeezing the colorbar
    #tr = [650,900]

#-------------------------------------------------------------
#-------------------------------------------------------------



##if loading a pickle file, then load this file
#if leg == 'upleg':
#    picklefile = 'pol_slice_values_upleg.pkl'
#if leg == 'downleg':
#    picklefile = 'pol_slice_values_downleg.pkl'

#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))


#load the flh and fcH, fcO lines to overplot
tmpfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/IGRF/IGRF_fcH_flh_lines_for_plotting.pkl'
plot_times,flh,fcHIGRF,fcOIGRF = pickle.load(open(tmpfile,'rb'))

##only use flh values in this range (outside of it they're not very accurate)
#flh_gootimes = np.where((plot_times > 150) & (plot_times < 825))


#--------------------------------------------------------------------------------------------------
#Read in polarization data of Bernstein waves from IDL save file
idldat = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/polarization_from_idl/pol_from_idl_fft=4096.sav')

timesIDL = idldat['times_pow']
freqs_pow = idldat['freqs_pow']
data_pow = idldat['data_pow']
freqs_pol = idldat['freqs_pol']
data_pol = idldat['data_pol']
freqs_elip = idldat['freqs_elip']
data_elip = idldat['data_elip']
freqs_hel = idldat['freqs_hel']
data_hel = idldat['data_hel']


#Reduce these to limited freq range (DON'T CHANGE - NOT RELATED TO THE PLOTTING)
good = np.where((freqs_pow > fcolorrange[0]) & (freqs_pow < fcolorrange[1]))[0]
data_pow = data_pow[good,:]
data_pol = data_pol[good,:]
data_elip = data_elip[good,:]
data_hel = data_hel[good,:]
freqs_pow = freqs_pow[good]

#--------------------------------------------------------------------------------------------------

data_pow2 = np.copy(data_pow)
data_elip2 = np.copy(data_elip)
data_pol2 = np.copy(data_pol)

#For each time, grab the distribution of wave polarization and power. 
#goo = np.where(timesIDL > 140)



for i in range(len(timesIDL)):
    bad = np.where(data_pow[:,i] < powthres)[0]
    data_pow2[bad,i] = np.nan
    data_elip2[bad,i] = np.nan
    data_pol2[bad,i] = np.nan




ttmp = np.zeros(np.shape(data_pow))
for i in range(len(freqs_pow)):
    ttmp[i,:] = timesIDL
ctmp = np.zeros(np.shape(data_pow))
for i in range(len(timesIDL)):
    ctmp[:,i] = freqs_pow
    bad = np.where(data_pow[:,i] < powthres)[0]
    ctmp[bad,i] = np.nan



goot = np.where((timesIDL > tr[0]) & (timesIDL < tr[1]))[0]


#cmap = 'turbo'
cmap = 'viridis'
#cmap = 'gist_heat'
fig,axs = plt.subplots(3,figsize=(7,9))
for i in range(3):
    axs[i].set_xlim(tr[0],tr[1])
axs[0].set_ylim(fr[0],fr[1])
axs[1].set_ylim(0.7,1)
axs[2].set_ylim(er[0],er[1])
scatter = axs[0].scatter(ttmp[:,goot],ctmp[:,goot],c=ctmp[:,goot], cmap=cmap,s=4, alpha=1)
cbar = plt.colorbar(scatter)
cbar.set_label('freq(Hz)',rotation=270,labelpad=15)
scatter.set_alpha(0.03)
axs[0].set_ylabel('spectrum of VLF power')
axs[0].set_xticklabels([])

scatter = axs[1].scatter(ttmp[:,goot],data_pol2[:,goot],c=ctmp[:,goot], cmap=cmap,s=4, alpha=1)
cbar = plt.colorbar(scatter)
cbar.set_label('freq(Hz)',rotation=270,labelpad=15)
scatter.set_alpha(0.03)
axs[1].set_ylabel('degree polarization')
axs[1].set_xticklabels([])
scatter = axs[2].scatter(ttmp[:,goot],data_elip2[:,goot],c=ctmp[:,goot], cmap=cmap,s=4, alpha=1)
cbar = plt.colorbar(scatter)
cbar.set_label('freq(Hz)',rotation=270,labelpad=15)
scatter.set_alpha(0.03)
axs[2].set_ylabel('ellipticity')
axs[2].set_xlabel('time(sec)')
fig.tight_layout(pad=1)


for i in range(3): 
    if leg == 'upleg':
        axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='black',s=12)
        axs[i].scatter(verticesLOW[9:,0], verticesLOW[9:,1],color='magenta',s=3)
    if leg == 'downleg':
        axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='black',s=12)
        axs[i].scatter(verticesHIG[5:-13:,0], verticesHIG[5:-13:,1],color='magenta',s=3)
    
    #axs[i].plot(plot_times, fcOIGRF, color='white')
    for j in range(11):
        axs[i].plot(plot_times,j*fcHIGRF,color='black',linestyle='--',linewidth=0.9) 
      


plt.savefig("/Users/abrenema/Desktop/tst1.pdf", dpi=350)




print('h')











#for i in range(1000):
#    ttmp[:] = timesIDL[i]
#    dtmp = np.squeeze(data_pow2[:,i])
#    #etmp = np.squeeze(data_elip2[:,i])
#    axs.scatter(ttmp,,c=cval,cmap='gist_heat',s=16)




axs[0].scatter(timesIDL[gootspec],freqsE,c=cval,s=4,cmap='gist_heat',zorder=2)


axs[1].scatter(timesIDL[gootspec],planE,c=cval,s=4,cmap='gist_heat',zorder=2)



plt.scatter(dtmp,etmp,c=cval,cmap='gist_heat',s=16)
plt.xlim(1e-16,1e-12)
plt.xscale('log')



freqscale = deepcopy(freqsE)
cval = freqscale
#cval = np.log(powscale)

#Plot time profile of values
fig, axs = plt.subplots(4)
axs[0].scatter(timesIDL[gootspec],freqsE,c=cval,s=16,cmap='gist_heat')















ephem = end_data_loader.load_ephemeris()
alts = ephem[0]['Altitude'] 
times = ephem[0]['Time']



v34 = EFL('VLF34D')
fs = v34.chnspecs['fs']
wf34, tdat = v34.load_data_gainphase_corrected()
#mV/m units

if leg == 'upleg':
    tr = [100,240]
if leg == 'downleg':
    tr = [700,898]

fr = [4000,9000]

#altitude at the desired time
gooalt = np.where(times >= tr[0])[0][0]
altz = alts[gooalt]

#flh (by-eye) at the desired time 
gooflh = np.where(verticesLOW[:,0] >= tr[0])[0][0]
flhz = verticesLOW[gooflh,1]
#--------------------------------------
#Get cyclotron frequencies

glat = 78 + (55/60)
glon = 11 + (55/60)
#Bo = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])
Bo = np.asarray(pyIGRF.igrf_value(glat, glon, altz, 2022)[6])

fce = 28 * Bo
fcH = fce / 1836
fcO = fcH / 16
fcHe = fcH / 4


#--------------------------------------

#Load SLP density spectrogram 
tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()

tdens = np.asarray(tdens)
fdensH = np.asarray(fdensH)
sdensH = np.asarray(sdensH)

#Get relevant time arrays
goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))[0]
gootspec = np.where((timesIDL > tr[0]) & (timesIDL < tr[1]))[0]
gootdens = np.where((tdens > tr[0]) & (tdens < tr[1]))[0]

#get relevant frequency locations
goofspec = np.where((freqs_pow > fr[0]) & (freqs_pow < fr[1]))[0]


##nfft = 4096
#nfft = 2048
##nfft = 8192
#psd12, psdf = correlation_analysis.psd(wf12[goot], tdat[goot], fs, tr, nft=nfft)
##units of mV/m**2 / sqrt(Hz)
##Change to dB of V/m**2 / sqrt(Hz) to be more consistent with most publications
##psd12 = psd12 / 1000000


#number of seconds to average over in slice
#nsec = tr[1]-tr[0]
nsec = 0

#Test run to get array sizes
elipA, elipARR, elipT = ps.slice_spectrogram(tr[0], timesIDL, data_elip, nsec=nsec)
densA, densARR, densT = ps.slice_spectrogram(tr[0], tdens, np.transpose(sdensH), nsec=nsec)


elipA = np.zeros([len(elipA), len(gootspec)])
powA = copy.deepcopy(elipA)
planA = copy.deepcopy(elipA)

#densA = np.zeros([len(densA), len(gootdens)])


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
if load_pickle == True:
    #Load pickle values for either upstream or downstream
    fnload = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/polarization_sliced_data/'+picklefile
    loadvals = pickle.load(open(fnload, 'rb'))
    elipA = loadvals[0]
    planA = loadvals[1]
    powA = loadvals[2]

else:
    #slice through each time to get polarization properties (entire freq range)
    #NOTE: this loop takes a while to run, so I'll save data to pickle file
    for t in range(len(gootspec)):
        elipA[:,t], elipARR, elipT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_elip, nsec=nsec)
        planA[:,t], planARR, planT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_pol, nsec=nsec)
        powA[:,t], powARR, powT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_pow, nsec=nsec)
    fnsav = '/Users/abrenema/Desktop/tmp'
    pickle.dump([elipA,planA,powA], open(fnsav + '.pkl','wb'))
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------



#Extract values corresponding to frequency of max power (within freq range of interest) at each time
powE = np.zeros(len(gootspec))
elipE = np.zeros(len(gootspec))
planE = np.zeros(len(gootspec))
freqsE = np.zeros(len(gootspec))
#densE = np.zeros(len(gootdens))

for t in range(len(gootspec)):
    powloc = np.where(powA[goofspec,t] == np.max(powA[goofspec,t]))[0]  
    powE[t] = powA[goofspec[powloc],t]
    elipE[t] = elipA[goofspec[powloc],t]
    planE[t] = planA[goofspec[powloc],t]
    freqsE[t] = freqs_pow[goofspec[powloc]]

badv = powE < powthres
powE[badv] = np.nan
elipE[badv] = np.nan
planE[badv] = np.nan
freqsE[badv] = np.nan

#powE[powE < powthres] = np.nan
#elipE[powE < powthres] = np.nan
#planE[powE < powthres] = np.nan
#freqsE[powE < powthres] = np.nan


#Define color scale based on amplitude while avoiding high amplitude noise at certain times 
#NOTE: this does not remove any data. 
if leg == 'upleg':
    goop = np.where((timesIDL[gootspec] > 112) & (timesIDL[gootspec] < 180))[0]
if leg == 'downleg':
    goop = np.where((timesIDL[gootspec] > 725) & (timesIDL[gootspec] < 750))[0]

powmax = np.nanmax(powE[goop])
powscale = deepcopy(powE)
powscale[powscale > powmax] = powmax


freqscale = deepcopy(freqsE)
cval = freqscale
#cval = np.log(powscale)

#Plot time profile of values
fig, axs = plt.subplots(4)
axs[0].scatter(timesIDL[gootspec],freqsE,c=cval,s=16,cmap='gist_heat')
axs[1].scatter(timesIDL[gootspec],powE,c=cval,s=16,cmap='gist_heat')
axs[2].scatter(timesIDL[gootspec],planE,c=cval,s=16,cmap='gist_heat')
axs[3].scatter(timesIDL[gootspec],elipE,c=cval,s=16,cmap='gist_heat')
axs[0].set_title('end_analysis_PSD_1d_statistics.py')
axs[0].set_ylim(4000,9000)
axs[0].set_ylabel('freq (Hz)\nof peak power')
axs[1].set_yscale('log')
axs[1].set_ylim(1e-14,0.25e-10)
axs[1].set_ylabel('Power')
axs[2].set_ylim(0.6,1.05)
axs[2].set_ylabel('deg pol')
axs[3].set_ylabel('ellipticity')
axs[3].set_xlabel('time (sec)')

#Plot freq vs ellipticity
plt.figure(figsize=(8,8))
plt.scatter(freqsE,elipE,c=cval,cmap='gist_heat')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz\nend_analysis_PSD_1d_statistics.py')
plt.xlabel('freq (Hz)\nof peak power')
plt.ylabel('ellipticity')
#Plot power vs ellipticity
plt.figure(figsize=(8,8))
plt.scatter(powE,elipE,c=cval,cmap='gist_heat')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz')
plt.xlabel('power')
plt.ylabel('ellipticity')
plt.xlim(0,0.25e-10)
#Plot deg polarization vs ellipticity
plt.figure(figsize=(8,8))
plt.scatter(planE,elipE,c=cval,cmap='gist_heat')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz')
plt.xlabel('deg pol')
plt.ylabel('ellipticity')

plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)





#For each time calculate statistics of quantities for all values where the power exceeds a threshold (within freq range of interest)
#powE = np.zeros(len(gootspec))
elipEmin = np.zeros(len(gootspec))
elipEmax = np.zeros(len(gootspec))
elipEmed = np.zeros(len(gootspec))
planEmin = np.zeros(len(gootspec))
planEmax = np.zeros(len(gootspec))
planEmed = np.zeros(len(gootspec))
freqEmin = np.zeros(len(gootspec))
freqEmax = np.zeros(len(gootspec))



for t in range(len(gootspec)):
    #powthres = 500*np.min(powA[goofspec,t])
    powloc = 0
    powloc = np.where((powA[goofspec,t] >= powthres) & (planA[goofspec,t] > 0.6))[0]
    if len(powloc) > 1:
        elipEmin[t] = np.min(elipA[goofspec[powloc],t])
        elipEmax[t] = np.max(elipA[goofspec[powloc],t])
        elipEmed[t] = np.median(elipA[goofspec[powloc],t])
        planEmin[t] = np.min(planA[goofspec[powloc],t])
        planEmax[t] = np.max(planA[goofspec[powloc],t])
        planEmed[t] = np.median(planA[goofspec[powloc],t])
        freqEmin[t] = freqs_pow[goofspec[np.min(powloc)]]
        freqEmax[t] = freqs_pow[goofspec[np.max(powloc)]]

elipEmin[elipEmin == 0] = np.nan
elipEmed[elipEmed == 0] = np.nan
elipEmax[elipEmax == 0] = np.nan
planEmin[planEmin == 0] = np.nan
planEmed[planEmed == 0] = np.nan
planEmax[planEmax == 0] = np.nan
freqEmin[freqEmin == 0] = np.nan
freqEmax[freqEmax == 0] = np.nan




fig, axs = plt.subplots(3,figsize=(8,8))
#fig.tight_layout(pad=0)
axs[0].set_title('Values above power threshold='+str(powthres)+' plotted\nbottomfreq(purple), topfreq(blue), dots(maxvalue)\ncoh>0.6 only\nend_analysis_PSD_1d_statistics.py')
axs[0].plot(timesIDL[gootspec],freqEmin,"+",color='mediumpurple',markersize=8)
axs[0].plot(timesIDL[gootspec],freqEmax,"+",color='cornflowerblue',markersize=8)
axs[0].scatter(timesIDL[gootspec],freqsE,c=cval,s=4,cmap='gist_heat',zorder=2)
axs[0].set_ylabel('freq(Hz)')
axs[0].set_ylim(fr[0],fr[1])
axs[0].plot(verticesLOW[9:,0], verticesLOW[9:,1],color='black',linewidth=3)
axs[0].plot(verticesHIG[5:,0], verticesHIG[5:,1],color='black',linewidth=3)
axs[0].plot(plot_times, fcHIGRF,color='black')
axs[0].plot(plot_times, fcOIGRF,color='black')
axs[0].plot(plot_times, 2*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 3*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 4*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 5*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 6*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 7*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 8*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 9*fcHIGRF,linestyle='--',color='black')
axs[0].plot(plot_times, 10*fcHIGRF,linestyle='--',color='black')
#axs[0].set_xticklabels([])

axs[1].plot(timesIDL[gootspec],planEmin,'+',color='mediumpurple',markersize=8)
axs[1].plot(timesIDL[gootspec],planEmax,'+',color='cornflowerblue',markersize=8)
axs[1].scatter(timesIDL[gootspec],planE,c=cval,s=4,cmap='gist_heat',zorder=2)
axs[1].set_ylim(0.6,1)
axs[1].set_ylabel('deg pol')
#axs[1].set_xticklabels([])

axs[2].plot(timesIDL[gootspec],elipEmin,'+',color='mediumpurple',markersize=8)
axs[2].plot(timesIDL[gootspec],elipEmax,'+',color='cornflowerblue',markersize=8)
axs[2].scatter(timesIDL[gootspec],elipE,c=cval,s=4,cmap='gist_heat',zorder=2)
axs[2].set_ylabel('ellipticity')
fig.tight_layout(pad=2)

if leg == 'upleg':
    for i in range(3):
        axs[i].set_xlim(100,240)
else:
    for i in range(3):
        axs[i].set_xlim(700,900)

# define color map
cmap = cm.get_cmap("gist_heat")
norm = colors.Normalize(4000, 9000)
for i in range(3):
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),ax=axs[i],label='frequency (Hz)')




plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)


