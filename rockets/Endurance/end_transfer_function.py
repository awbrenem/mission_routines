"""
Calculate transfer function for Endurance channels

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians
(3) Apply transfer function correction as: fft_data_corrected = fft_data/H


(see https://resources.pcb.cadence.com/blog/2021-understanding-a-circuit-transfer-function-from-a-bode-plot)
"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import matplotlib.pyplot as plt
import numpy as np
import end_load_data as end
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
import end_load_gainphase as gainphase
from scipy import signal
import filter_wave_frequency as filt



#-------------------------------------------------------
#Select data type (DC or VLF) and channel (12, 34, etc) for calibration 
#Load data...
#-------------------------------------------------------


type = 'VLF'
chn = '12'
filename = "Endurance_Analog 1_" + type + chn + "D_6-30000-100.txt"
pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_VLF/'
wavegoo = end.efield_vlf()
tdat = wavegoo.tvlf  #times
if chn == '12': wavedat = wavegoo.dvlf12_mvm
elif chn == '34': wavedat = wavegoo.dvlf34_mvm
elif chn == '24': wavedat = wavegoo.dvlf24_mvm
elif chn == '32': wavedat = wavegoo.dvlf32_mvm

#----------------------------------------------------
"""
type = 'DC'
chn = '12'
filename = "Endurance_Analog 1_V"+chn+"D_10-10000-100.txt"
pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/efield_DC/'
wavegoo = end.efield_dc()
tdat = wavegoo.times
if chn == '12': wavedat = wavegoo.dv12_mvm
elif chn == '34': wavedat = wavegoo.dv34_mvm
elif chn == '13': wavedat = wavegoo.dv13_mvm
elif chn == '32': wavedat = wavegoo.dv32_mvm
elif chn == '24': wavedat = wavegoo.dv24_mvm
elif chn == '41': wavedat = wavegoo.dv41_mvm
"""

fs = wavegoo.samplerate




#-----------------------------------------------------------------------------------------
#Load gain/phase data for selected channel. 
#-----------------------------------------------------------------------------------------


prad, Hmag, f = gainphase.end_load_gainphase(filename)


#Chop off calibration values above the Nyquist freq of data.
#Not doing this can mess up the interpolations below...

f = np.asarray(f)
prad = np.asarray(prad)
Hmag = np.asarray(Hmag)

f2 = np.asarray(f < fs)
f = f[f2]
prad = prad[f2]
Hmag = Hmag[f2]




#At this point we have an issue with the DC data. The interpolated frequency array (freq) does go down to 0 Hz, however
#the gain and phase have NaN values at < 10 Hz (or similar) b/c the AC testing setup can't go to 0 Hz in freq. 
#Therefore we need to assign these NaN values a finite value based on the values at ~10 Hz. 
#NOTE: need to do multiple values <10 Hz or else the interpolation to higher resolution (later) gets wonky

if type == 'DC':
    newfgoo = np.arange(10)
    Hgoo = np.full(10,Hmag[0]) 
    Pgoo = np.full(10,prad[0])
    f = np.insert(f,0,newfgoo)
    Hmag = np.insert(Hmag,0,Hgoo)
    prad = np.insert(prad,0,Pgoo)

if type == 'VLF':
    #average value b/t 100-1000 Hz 
    #goo = list(map(tuple, np.where((f > 100) & (f < 1000))))
    #Hmag_avg = np.mean(Hmag[goo])
    #badgain_ind = list(map(tuple, np.where(f < 30)))
    #Hmag[badgain_ind] = Hmag_avg

    #Change any gain values < 1 to 1 so that low freqs are not amplified.     
    Hmag[Hmag < 1] = 1


    #Remove values < 10 Hz b/c I don't trust these. The signal/noise ratio is very high in the 
    #cal tests, and these gain/phase values can be all over the place. 
    fmin_trust = 30
    tt = f < fmin_trust
    f = np.delete(f, f < fmin_trust)
    Hmag = np.delete(Hmag, tt)
    prad = np.delete(prad, tt)

    #Change the < fmin_trust (Hz) values to the gain/phase just above fmin_trust
    #NOTE: best to have these low values as continuous as possible with higher freq values otherwise the 
    #interpolation can have unrealistic peaks. I'll bandpass the data later to remove any remaining low freq power.
    newfgoo = np.arange(fmin_trust)
    Hnewval = Hmag[0]
    Pnewval = prad[0]
    Hgoo = np.full(fmin_trust,Hnewval) 
    Pgoo = np.full(fmin_trust,Pnewval)
    f = np.insert(f,0,newfgoo)
    Hmag = np.insert(Hmag,0,Hgoo)
    prad = np.insert(prad,0,Pgoo)
    







#-------------------------------------------------------------
#FFT Endurance waveform data to apply transfer function
#-------------------------------------------------------------


wavedatFFT = rfft(wavedat)


N = len(wavedatFFT)
n = np.arange(N)
T = N/fs
freq = n/T/2 

#fig, axs = plt.subplots(3)
#axs[0].plot(freq,np.abs(wavedatFFT))
#axs[1].plot(freq,np.real(wavedatFFT))
#axs[2].plot(freq,np.imag(wavedatFFT))
#for i in range(3): axs[i].set_xscale('linear')
#for i in range(3): axs[i].set_xlim(0,100)





#------------------------------------------------------------------
#Transfer function H(w) = |H|*exp(i*theta)
#------------------------------------------------------------------


# The interpolated phase has an over/undershoot at 1 kHz. This is due to interpolating over a sudden change.
# To get around this, need to unwrap values, interpolate, then rewrap
prad_unwrapped = np.unwrap(prad)
#plt.plot(f,prad, f,prad_unwrapped)
#plt.xscale('log')


#Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(f,Hmag,kind='cubic', bounds_error=False)
Hmag2 = interp(freq)
interp2 = interp1d(f,prad_unwrapped,kind='cubic', bounds_error=False)
prad2_unwrapped = interp2(freq)


#plt.plot(freq,Hmag2)
#plt.plot(freq,prad2)
#plt.xlim(0,1000)


#-------------------------
#rewrap angles from -2pi to 2pi 
#Note that there are some minor interpolation issues here due to granularity of data (shoots past 2*pi)
#THIS ISN'T AN ISSUE. 
prad2 = (prad2_unwrapped + np.pi) % (2 * np.pi) - np.pi
plt.plot(f, prad, freq,prad2_unwrapped, freq, prad2)
plt.xscale('log')
plt.ylim(-4,4)
print('check wrapping')




fig, axs = plt.subplots(2)
axs[0].plot(f, Hmag, freq, Hmag2)
axs[1].plot(f, prad, freq, prad2)
for i in range(2): axs[i].set_xscale('log')
for i in range(2): axs[i].set_xlim(1,30000)
axs[1].set_ylim(-4,4)

print('Check modified gain/phase curves')

#---------------------------------------------------
#Define the transfer function.
#Two versions:
#1) gain + phase
#2) unity gain + phase 

#...note that the files Steve sent me for Endurance already have the DC gain correction applied
#...for both the DC and VLF channels (Paulo must have done this). So when applying the transfer function
#...to these data use the unity gain version. 
#---------------------------------------------------

H = [Hmag2[i] * np.exp(1j*prad2[i]) for i in range(len(prad2))]

#fig2, axs2 = plt.subplots(3)
#axs2[0].plot(freq,np.abs(H))
#axs2[1].plot(freq,np.real(H))
#axs2[2].plot(freq,np.imag(H))
#for i in range(3): axs2[i].set_xscale('log')
#print('Check transfer function')


maxv = np.nanmax(Hmag2)
Hmag2_unity = [i/maxv for i in Hmag2]
H_unitygain = [Hmag2_unity[i] * np.exp(1j*prad2[i]) for i in range(len(prad2))]



"""
#Bode plots
fig2, axs2 = plt.subplots(3)
axs2[0].plot(freq,np.abs(H))
axs2[1].plot(freq,np.real(H))
axs2[2].plot(freq,np.imag(H))
for i in range(3): axs2[i].set_xscale('log')
for i in range(3): axs2[i].set_xlim(1,20000)
for i in range(3): axs2[i].set_xlabel('freq(Hz)')
axs2[0].set_title('Bode plots \n ' + filename)
axs2[0].set_ylabel('Magnitude |H(jw)|')
axs2[1].set_ylabel('Re(H(jw))')
axs2[2].set_ylabel('Phase Im(H(jw))')
"""


print('Check unity gain transfer function')



#---------------------------------------------------
#Apply transfer function to FFT'd data (amplitude, not power)
#---------------------------------------------------

wavedatFFTc = [wavedatFFT[i]/H_unitygain[i] for i in range(len(H_unitygain))]


#----------
#----------
#----------
#UNDER CONSTRUCTION - need to test this
#Remove very low freqs due to remaining artificial power
#tst = freq < 20
#wavedatFFTc2 = wavedatFFTc
wavedatFFTc[freq < 20] = 0-0j
#----------
#----------
#----------



freq2 = list(freq)
v1 = np.abs(wavedatFFT)
v2 = np.abs(wavedatFFTc)
rat = [v1[i]/v2[i] for i in range(len(v1))]

fig3,axs3 = plt.subplots(4)
axs3[0].plot(freq2,np.abs(wavedatFFT))
axs3[1].plot(freq2,np.abs(wavedatFFTc))
axs3[0].set_ylim(0,3e6)
axs3[1].set_ylim(0,3e6)
axs3[2].plot(freq2,rat)
axs3[2].set_yscale('linear')
axs3[3].plot(freq2,np.abs(H))
axs3[3].set_yscale('linear')
for i in range(4): axs3[i].set_xlim(200,201)
print('Quick look at calibrated vs original data')





#---------------------------------------------------
#Inverse FFT to get back to corrected waveform
#---------------------------------------------------

wf = irfft(wavedatFFT, n=len(tdat))
wf_corr = irfft(wavedatFFTc, n=len(tdat))



##Remove low freq power in the VLF channels. This can result in unrealistic sine wave at low freqs. 
#wf_corr = filt.butter_highpass_filter(wf_corr, 20, fs, order=5)




fig,axs = plt.subplots(3)
axs[0].plot(tdat,wavedat)
axs[1].plot(tdat,wf)
axs[2].plot(tdat,wf_corr)
#for i in range(3): axs[i].set_ylim(-100,100)
for i in range(3): axs[i].set_xlim(400,402)
for i in range(3): axs[i].set_ylim(-1,1)



import plot_spectrogram as ps
fspec, tspec, powerc = signal.spectrogram(wavedatFFT, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
#ps.plot_spectrogram(tspec,fspec,powerc,vr=[-50,-35],yr=[1,1000],xr=[100,800], yscale='log')
#Before I inverse transform I need to get rid of NaN values at the beginning of new data arrays
bad = list(map(tuple, np.where(np.isnan(wavedatFFTc))))
wavedatFFTc = np.asarray(wavedatFFTc)
wavedatFFTc[bad] = 0





#print('quick look at calibrated vs original waveforms')



#----------------------------------------------------
#Save corrected waveform data
#----------------------------------------------------

import pickle

dict_fin = {'tvals':tdat, 'wf':wf_corr}

fnsav = filename[:-4] + '_gainphase_corrected'
pickle.dump(dict_fin, open(pathoutput + fnsav + '.pkl','wb'))
#wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))


