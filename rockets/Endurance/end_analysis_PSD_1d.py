#Plot the 1D power spectral density
#and make plot for paper

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
from matplotlib.ticker import MultipleLocator
from scipy.signal import find_peaks
import pyIGRF
from scipy.io import readsav
import pickle
import correlation_analysis as ca
import filter_wave_frequency

#Load by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))




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
#--------------------------------------------------------------------------------------------------


ephem = end_data_loader.load_ephemeris()
alts = np.asarray(ephem[0]['Altitude'])

#interpolate alts to times of spectra
alts2 = np.interp(timesIDL,np.asarray(ephem[0]['Time']),alts)

times = ephem[0]['Time']



v34 = EFL('VLF34D')
fs = v34.chnspecs['fs']
wf34, tdat = v34.load_data_gainphase_corrected()

v12 = EFL('VLF12D')
wf12, tdat = v12.load_data_gainphase_corrected()


#mV/m units

#tr = [130,135]
#tr = [130,170]
#tr = [140,145]
#tr = [150,155]
#tr = [160,165]
#tr = [160,161]
#tr = [170,175]
#tr = [210,215]
tr = [320,340]
#tr = [710,715]
#tr = [750,755]
#tr = [830,835]
#tr = [810,815]
#tr = [803,805]
goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
#window = np.hanning(len(goot[0]))

#altitude at the desired time
gooalt = np.where(times >= tr[0])[0][0]
altz = alts[gooalt]

#all altitudes during the desired time range
gooalts = np.where((timesIDL >= tr[0]) & (timesIDL <= tr[1]))[0]
#altv = np.asarray(alts[gooalts])
altv = np.asarray(alts2[gooalts])



#flh (by-eye) at the desired time 
gooflh = np.where(verticesHIG[:,0] >= tr[0])[0][0]
flhz = verticesHIG[gooflh,1]
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


#find fce values for all altitudes during timerange 
Bo_all = np.empty(len(altv))
for i in range(len(altv)):
       Bo_all[i] = np.asarray(pyIGRF.igrf_value(glat, glon, altv[i], 2022)[6])
fce_all = Bo_all * 28
fcH_all = fce_all / 1836


#----------------------
#----------------------
#----------------------
#interpolate the fcH_all array size to the number of elements of the final freq array. 
#This will allow us to plot f/fcH instead of just f
#freqs_pow
#dt = tr[1] - tr[0]
#oldtimes = np.arange(tr[0],tr[1],1/len(altv))

f_fce = np.zeros((len(freqs_pow),len(gooalts)))
for i in range(len(gooalts)):
       f_fce[:,i] = freqs_pow / fcH_all[i]

#alts2 




#----------------------
#----------------------
#----------------------

#--------------------------------------
#Calculate coherence
timechunk = 0.2
coh, phase, tchunks, ff = ca.cross_spectral_density_spectrogram(wf34[tdat > 100],wf12[tdat > 100],tdat[tdat > 100],30000,timechunk,nperseg=512,plot=False,coh_min=0)



#Load SLP density spectrogram 
tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()
tdens = np.squeeze(np.asarray(tdens))
fdensH = np.squeeze(np.asarray(fdensH))
sdensH = np.squeeze(np.asarray(sdensH))


#nfft = 4096
nfft = 2048
#nfft = 8192
psd12, psdf = correlation_analysis.psd(wf12[goot], tdat[goot], fs, tr, nft=nfft)
#units of mV/m**2 / sqrt(Hz)
#Change to dB of V/m**2 / sqrt(Hz) to be more consistent with most publications
#psd12 = psd12 / 1000000
psd12 = np.interp(freqs_pow,psdf,psd12)


#number of seconds to average over in slice
nsec = tr[1]-tr[0]
#nsec = 0

#Slice through the ellipticity and deg polarization arrays
elipA, elipARR, elipT = ps.slice_spectrogram(tr[0], timesIDL, data_elip, nsec=nsec)
planA, planARR, planT = ps.slice_spectrogram(tr[0], timesIDL, data_pol, nsec=nsec)
#f_fceA, f_fceARR, f_fceT = ps.slice_spectrogram(tr[0], timesIDL[gooalts], f_fce, nsec=nsec)

#slice through SLP density sonogram and interpolate to times of the polarization data
densA, densARR, densT = ps.slice_spectrogram(tr[0], tdens, np.transpose(sdensH), nsec=nsec)



#High pass the density data to get rid of DC offset 
#densA = filter_wave_frequency.butter_highpass_filter(densA, 1390, fdensH, order=5)
gooloc = np.where(fdensH > 4000)[0]
densA -= densA[gooloc[0]]




cohA, cohARR, cohT = ps.slice_spectrogram(tr[0], tchunks, coh, nsec=nsec)
#interpolate coherence to times of the polarization data. 
cohA2 = np.interp(freqs_pow,ff,cohA)



#Change to np arrays
elipA = np.asarray(elipA)
psd12 = np.asarray(psd12)

#Highlight data where deg polarization > 0.7 and coherence is > 0.4
planA = np.asarray(planA)

#goo = np.where((planA < 0.7) | (cohA2 < 0.4))[0]
goo = np.where(planA < 0.7)[0]
psd12z = np.copy(psd12)
psd12z[goo] = np.nan

planAz = np.copy(planA)
planAz[goo] = np.nan

elipAz = np.copy(elipA)
elipAz[goo] = np.nan

cohA2z = np.copy(cohA2)
cohA2z[goo] = np.nan


#Remove abnormally high values near DC
cohA2[0:20] = np.nan
cohA2z[0:20] = np.nan

psd12[0:20] = np.nan
planA[0:20] = np.nan

elipA[0:20] = np.nan
elipAz[0:20] = np.nan



fig, axs = plt.subplots(4,figsize=(6,7))

axs[0].set_title('end_analysis_PSD_1d.py for t=['+str(tr[0])+ '-' + str(tr[1]) + '] sec; nfft=' + str(nfft) + '\n fcO='+str(np.floor(fcO))+'; fcHe='+str(np.floor(fcHe))+'; fcH='+str(np.floor(fcH))+' Hz')
axs[0].plot(freqs_pow, 10*np.log(psd12),color='gray')
axs[0].plot(freqs_pow, 10*np.log(psd12z),color='red')
axs[0].set_ylim([-100,-50])
plt.minorticks_on()
axs[0].set_ylabel('PSD\ndB of (mV/m)^2/Hz')
axs[0].set_xticklabels([])

#axs[1].plot(freqs_pow, cohA2, color='gray')
#axs[1].plot(freqs_pow, cohA2z, color='red')
#axs[1].set_ylabel('Coherence\n(VLF12, VLF34)')
#axs[1].set_xticklabels([])
#axs[1].set_ylim(0,1)


axs[1].plot(freqs_pol,planA,color='gray')
axs[1].plot(freqs_pol,planAz,color='red')
axs[1].set_yscale('linear')
axs[1].set_ylim([0.5,1])
plt.minorticks_on()
axs[1].set_ylabel('deg polarization')
axs[1].set_xticklabels([])

axs[2].plot(freqs_elip,elipA,color='gray')
axs[2].plot(freqs_elip,elipAz,color='red')
axs[2].set_yscale('linear')
axs[2].set_ylim([-0.8,0.4])
axs[2].set_ylabel('Ellipticity')
axs[2].set_xticklabels([])

axs[2].axhline(y=0.0, color='black',linestyle='--')


axs[3].plot(fdensH,densA,color='black')
axs[3].set_yscale('linear')
axs[3].set_ylim(np.nanmin(densA),np.nanmax(densA))
axs[3].set_ylabel('SLP density\nrelative units')


plt.minorticks_on()
plt.xlabel('freq (Hz)')


for i in range(4):
       axs[i].vlines(fcH*6, -400, 10, color='black',linestyle='--')
       axs[i].vlines(fcH*7, -400, 10, color='black',linestyle='--')
       axs[i].vlines(fcH*8, -400, 10, color='black',linestyle='--')
       axs[i].vlines(fcH*9, -400, 10, color='black',linestyle='--')
       axs[i].vlines(fcH*10, -400, 10, color='black',linestyle='--')
       if tr[0] < 800:
              axs[i].vlines(flhz, -400, 10, color='magenta',linestyle='--')
       axs[i].set_xlim(4000,10000)
       axs[i].xaxis.set_minor_locator(MultipleLocator(100))   # minor ticks every 10

fig.tight_layout(pad=2)
plt.subplots_adjust(hspace=0.1)

plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)




plt.plot(fdensH,densA)









"""
for tr=[750-755]
plt.ylim([1e-8, 2e-6])
plt.xlim(4000,8000)
plt.axvline(4500,color='orange')
plt.axvline(4500 + fcO,color='orange')
plt.axvline(4500 + fcHe,color='orange')
plt.axvline(5174)
plt.axvline(5299)
plt.axvline(5458)
plt.axvline(5618)
plt.axvline(5794)
plt.axvline(5942)
plt.axvline(6182)
plt.axvline(6298)
plt.axvline(6414)
plt.axvline(6622)
plt.axvline(6752)
plt.axvline(6871)
plt.axvline(7018)
plt.axvline(7080)
plt.axvline(7295)
plt.axvline(fcH*5,color='green')
plt.axvline(fcH*6,color='green')
plt.axvline(fcH*7,color='green')
plt.axvline(fcH*8,color='green')
plt.axvline(fcH*9,color='green')
plt.axvline(fcH*10,color='green')
"""



"""
for tr=[170-175]
plt.ylim([1e-8, 1e-5])
plt.xlim(4000,7000)
plt.axvline(4700,color='orange')
plt.axvline(4700 + fcO,color='orange')
plt.axvline(4700 + fcHe,color='orange')
plt.axvline(5743)
plt.axvline(6083)
plt.axvline(6306)
plt.axvline(fcH*6,color='green')
plt.axvline(fcH*7,color='green')
plt.axvline(fcH*8,color='green')
plt.axvline(fcH*9,color='green')
plt.axvline(fcH*10,color='green')
"""






plt.plot(psdf, psd12**2)
plt.title('end_analysis_PSD_1d.py')
plt.yscale('log')
plt.ylim([1e-8, 1e-5])
plt.xlim(4000,8000)
plt.ylabel('PSD (mV/m**2/Hz)')
peaks, _ = find_peaks(psd12**2, height=1e-6, distance=5)
plt.plot(psdf[peaks], psd12[peaks]**2, "x")
#plt.plot(np.zeros_like(x), "--", color="gray")
plt.show()


print('h')




#Vm_rootHz = psd12 / 1000 
#Vm2_Hz = Vm_rootHz**2
#plt.plot(psdf, Vm2_Hz)
#plt.ylim([1e-14, 1e-11])

