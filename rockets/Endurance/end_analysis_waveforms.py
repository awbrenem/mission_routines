#Plot the 1D power spectral density
#and make plot for paper

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
from matplotlib.ticker import MultipleLocator
from scipy.signal import find_peaks
import filter_wave_frequency


v12 = EFL('VLF12D')
wf12, tdat = v12.load_data_gainphase_corrected()

v34 = EFL('VLF34D')
wf34, tdat = v34.load_data_gainphase_corrected()
#mV/m units


fs = 1/(tdat[1]-tdat[0])
wf12  = filter_wave_frequency.butter_highpass_filter(wf12,2000,fs,order=5)
wf34  = filter_wave_frequency.butter_highpass_filter(wf34,2000,fs,order=5)

#waveforms
tst = 141.
dt = 0.02
goot = np.where((tdat >= tst) & (tdat <= tst+dt))[0]
wf12z = wf12[goot]
wf34z = wf34[goot]
tz = tdat[goot]
tz -= tz[0]
tz *= 1000

ttitle = 'time [msec since ' + str(int(tdat[goot][0])) + ' sec]'

fig, axs = plt.subplots(1)

#plt.figure(figsize=(6,3))
axs.plot(tz,wf34z)
axs.set_ylabel('VLF34\n[mV/m]')
axs.set_xlabel(ttitle)


plt.ylabel('VLF34\n[mV/m]')
plt.xlabel(ttitle)
plt.margins(x=0.1,y=8)

#hodograms
tst = 141.
dth = 0.001
goot = np.where((tdat >= tst) & (tdat <= tst+dth))[0]
wf12z = wf12[goot]
wf34z = wf34[goot]
tz = tdat[goot]
maxv = np.max([wf12z,wf34z])

plt.plot(wf12z,wf34z)
plt.gca().set_aspect('equal')
plt.ylim(-1*maxv,maxv)
plt.xlim(-1*maxv,maxv)




print('h')


#def auto_correlation(wf):






























#altitude at the desired time
gooalt = np.where(times >= tr[0])[0][0]
altz = alts[gooalt]

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


#--------------------------------------

#Load SLP density spectrogram 
tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()

tdens = np.asarray(tdens)
fdensH = np.asarray(fdensH)
sdensH = np.asarray(sdensH)


#nfft = 4096
nfft = 2048
#nfft = 8192
psd12, psdf = correlation_analysis.psd(wf34[goot], tdat[goot], fs, tr, nft=nfft)
#units of mV/m**2 / sqrt(Hz)
#Change to dB of V/m**2 / sqrt(Hz) to be more consistent with most publications
#psd12 = psd12 / 1000000

#number of seconds to average over in slice
nsec = tr[1]-tr[0]
#nsec = 0

#Slice through the ellipticity and deg polarization arrays
elipA, elipARR, elipT = ps.slice_spectrogram(tr[0], timesIDL, data_elip, nsec=nsec)
planA, planARR, planT = ps.slice_spectrogram(tr[0], timesIDL, data_pol, nsec=nsec)
#slice through SLP density sonogram
densA, densARR, densT = ps.slice_spectrogram(tr[0], tdens, np.transpose(sdensH), nsec=nsec)

#tdens, fdensL, sdensL, fdensH, sdensH = end_data_loader.load_slp_sonogram()


#Highlight data where deg polarization > 0.7 
planA = np.asarray(planA)

goo = np.where(planA < 0.7)[0]
psd12z = np.copy(psd12)
psd12z[goo] = np.nan

planAz = np.copy(planA)
planAz[goo] = np.nan

elipAz = np.copy(elipA)
elipAz[goo] = np.nan


fig, axs = plt.subplots(4,figsize=(6,7))

axs[0].set_title('end_analysis_PSD_1d.py for t=['+str(tr[0])+ '-' + str(tr[1]) + '] sec; nfft=' + str(nfft) + '\n fcO='+str(np.floor(fcO))+'; fcHe='+str(np.floor(fcHe))+'; fcH='+str(np.floor(fcH))+' Hz')
axs[0].plot(psdf, 10*np.log(psd12),color='gray')
axs[0].plot(psdf, 10*np.log(psd12z),color='red')
#axs[0].set_yscale('log')
#axs[0].set_ylim([5e-11, 1e-8])
#axs[0].set_ylim([-230,-180])
axs[0].set_ylim([-100,-50])
plt.minorticks_on()
axs[0].set_ylabel('PSD\ndB of (mV/m)^2/Hz')
axs[0].set_xticklabels([])

axs[1].plot(fdensH,densA,color='black')
axs[1].set_yscale('linear')
axs[1].set_ylim([-15,-8])
#axs[1].set_ylim([-500,500])
axs[1].set_ylabel('SLP density\nrelative units')
plt.minorticks_on()
plt.xlabel('freq (Hz)')
axs[1].set_xticklabels([])

axs[2].plot(freqs_pol,planA,color='gray')
axs[2].plot(freqs_pol,planAz,color='red')
axs[2].set_yscale('linear')
axs[2].set_ylim([0.5,1])
plt.minorticks_on()
axs[2].set_ylabel('deg polarization')
axs[2].set_xticklabels([])

axs[3].plot(freqs_elip,elipA,color='gray')
axs[3].plot(freqs_elip,elipAz,color='red')
axs[3].set_yscale('linear')
axs[3].set_ylim([-0.4,0.8])
plt.minorticks_on()
plt.ylabel('Ellipticity')
plt.xlabel('freq (Hz)')
plt.axhline(y=0.0, color='black',linestyle='--')

for i in range(4):
       axs[i].vlines(fcH*6, -400, 1, color='black',linestyle='--')
       axs[i].vlines(fcH*7, -400, 1, color='black',linestyle='--')
       axs[i].vlines(fcH*8, -400, 1, color='black',linestyle='--')
       axs[i].vlines(fcH*9, -400, 1, color='black',linestyle='--')
       axs[i].vlines(fcH*10, -400, 1, color='black',linestyle='--')
       if tr[0] < 800:
              axs[i].vlines(flhz, -400, 1, color='magenta',linestyle='--')
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

