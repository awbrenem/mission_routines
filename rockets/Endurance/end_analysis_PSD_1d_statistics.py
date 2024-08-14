#Build up statistics of the 1D power spectral density vs time for various plots. 


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
from matplotlib.ticker import MultipleLocator
from scipy.signal import find_peaks
import pyIGRF
import scipy.io as sio
from scipy.io import readsav
import copy 

vertices2 = [[ 118.55265075, 7964.38255622],
       [ 126.59693574, 8254.90835695],
       [ 134.64122073, 8530.14332606],
       [ 145.90321972, 8514.85249444],
       [ 156.36079021, 8438.39833636],
       [ 162.7962182 , 8316.07168342],
       [ 170.8405032 , 8102.00004078],
       [ 179.68921669, 8010.25505107],
       [ 186.92907318, 7903.21922975],
       [ 199.79992917, 7796.18340843],
       [ 208.64864266, 7643.27509225],
       [ 223.12835565, 7551.53010255],
       [ 233.58592614, 7475.07594446],
       [ 243.23906813, 7199.84097535],
       [ 255.30549562, 7459.78511285],
       [ 266.56749461, 7398.62178638],
       [ 270.58963711, 7169.25931212],
       [ 281.8516361 , 7199.84097535],
       [ 293.11363508, 7062.22349079],
       [ 302.76677708, 7352.74929152],
       [ 311.61549057, 7291.58596506],
       [ 320.46420406, 7092.80515403],
       [ 330.92177455, 6955.18766947],
       [ 338.16163104, 7016.35099594],
       [ 348.61920153, 6939.89683786],
       [ 358.27234352, 6955.18766947],
       [ 367.12105702, 6832.86101653],
       [ 377.57862751, 6909.31517462],
       [ 388.036198  , 6894.024343  ],
       [ 399.29819699, 7138.67764888],
       [ 404.92919648, 6848.15184815],
       [ 409.75576748, 7092.80515403],
       [ 425.03990896, 6955.18766947],
       [ 434.69305095, 6848.15184815],
       [ 453.99933494, 6832.86101653],
       [ 468.47904792, 6817.57018492],
       [ 474.91447592, 6786.98852168],
       [ 499.0473309 , 6802.2793533 ],
       [ 511.91818688, 6863.44267977],
       [ 520.76690038, 6924.60600624],
       [ 532.02889936, 6955.18766947],
       [ 542.48646986, 6955.18766947],
       [ 565.81489633, 6924.60600624],
       [ 577.88132382, 7092.80515403],
       [ 588.33889431, 6939.89683786],
       [ 626.95146228, 7184.55014373],
       [ 641.43117526, 7230.42263859],
       [ 651.88874575, 7337.45845991],
       [ 664.75960174, 7429.20344961],
       [ 680.04374323, 7429.20344961],
       [ 693.71902772, 7520.94843932],
       [ 706.5898837 , 7444.49428123],
       [ 725.09173919, 7475.07594446],
       [ 736.35373818, 7475.07594446],
       [ 746.81130867, 7612.69342902],
       [ 763.70430715, 7673.85675549],
       [ 775.77073464, 7811.47424005],
       [ 786.22830513, 7903.21922975],
       [ 797.49030412, 8086.70920916],
       [ 807.14344611, 8224.32669371],
       [ 819.2098736 , 8407.81667312],
       [ 832.88515809, 8499.56166283],
       [ 841.73387158, 8468.97999959],
       [ 852.99587057, 8331.36251504],
       [ 857.01801307, 8147.87253563],
       [ 859.43129856, 7887.92839813],
       [ 863.45344106, 7536.23927093],
       [ 869.08444056, 7306.87679667],
       [ 871.49772605, 7016.35099594],
       [ 877.12872555, 6710.5343636 ],
       [ 877.93315405, 6312.97274154],
       [ 881.95529654, 5961.28361434]]

vertices2 = np.asarray(vertices2)


#--------------------------------------------------------------------------------------------------
#Read in polarization data of Bernstein waves from IDL save file
idldat = readsav('/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/polarization_from_idl/pol_from_idl.sav')

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
alts = ephem[0]['Altitude'] 
times = ephem[0]['Time']



v12 = EFL('VLF34D')
fs = v12.chnspecs['fs']
wf12, tdat = v12.load_data_gainphase_corrected()
#mV/m units

tr = [100,240]
#tr = [150,155]
#tr = [160,165]
#tr = [160,161]
#tr = [170,175]
#tr = [180,185]
#tr = [800,805]
#tr = [830,835]
#tr = [750,755]
#tr = [725,730]

fr = [4000,9000]

#altitude at the desired time
gooalt = np.where(times >= tr[0])[0][0]
altz = alts[gooalt]

#flh (by-eye) at the desired time 
gooflh = np.where(vertices2[:,0] >= tr[0])[0][0]
flhz = vertices2[gooflh,1]
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


#nfft = 4096
nfft = 2048
#nfft = 8192
psd12, psdf = correlation_analysis.psd(wf12[goot], tdat[goot], fs, tr, nft=nfft)
#units of mV/m**2 / sqrt(Hz)
#Change to dB of V/m**2 / sqrt(Hz) to be more consistent with most publications
psd12 = psd12 / 1000000


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

#slice through each time to get polarization properties (entire freq range)
for t in range(len(gootspec)):
    elipA[:,t], elipARR, elipT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_elip, nsec=nsec)
    planA[:,t], planARR, planT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_pol, nsec=nsec)
    powA[:,t], powARR, powT = ps.slice_spectrogram(timesIDL[gootspec[t]], timesIDL, data_pow, nsec=nsec)


#for t in range(len(gootdens)):
#    densA[:,t], densARR, densT = ps.slice_spectrogram(tdens[gootdens[t]], tdens, np.transpose(sdensH), nsec=nsec)



plt.plot(freqs_pol,elipA[:,10])
#plt.plot(fdensH,densA[:,10])



#Extract values at each time for frequency range of interest
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

    
#Plot time profile of values
fig, axs = plt.subplots(4)
axs[0].plot(timesIDL[gootspec],freqsE,'.')
axs[1].plot(timesIDL[gootspec],powE,'.')
axs[2].plot(timesIDL[gootspec],planE,'.')
axs[3].plot(timesIDL[gootspec],elipE,'.')
axs[0].set_title('end_analysis_PSD_1d_statistics.py')
axs[0].set_ylim(5000,7000)
axs[0].set_ylabel('freq (Hz)\nof peak power')
axs[1].set_ylim(0,0.25e-10)
axs[1].set_ylabel('Power')
axs[2].set_ylim(0.6,1.05)
axs[2].set_ylabel('deg pol')
axs[3].set_ylabel('ellipticity')
axs[3].set_xlabel('time (sec)')

#Plot freq vs ellipticity
plt.figure(figsize=(8,8))
plt.plot(freqsE,elipE,'.')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz\nend_analysis_PSD_1d_statistics.py')
plt.xlabel('freq (Hz)\nof peak power')
plt.ylabel('ellipticity')
#Plot power vs ellipticity
plt.figure(figsize=(8,8))
plt.plot(powE,elipE,'.')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz')
plt.xlabel('power')
plt.ylabel('ellipticity')
plt.xlim(0,0.25e-10)
#Plot deg polarization vs ellipticity
plt.figure(figsize=(8,8))
plt.plot(planE,elipE,'.')
plt.title('values at peak power for each time\n'+str(tr[0])+'-'+str(tr[1])+' sec\nfreqs=['+str(fr[0])+'-'+str(fr[1])+'] Hz')
plt.xlabel('deg pol')
plt.ylabel('ellipticity')




plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)

