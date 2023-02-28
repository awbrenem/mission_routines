"""
Calculate transfer function for Endurance channels

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians

(see https://resources.pcb.cadence.com/blog/2021-understanding-a-circuit-transfer-function-from-a-bode-plot)

"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import end_load_data as end
from scipy.interpolate import interp1d
import plot_spectrogram as ps


#-------------------------------
#Select desired gain/phase files (from Paulo) from the calibration testing.
#From these files Paulo derives a fit (ax + b) that allows calibration from counts to volts.
#These values (a, b) are found in the Endurance channel list document as the yellow boxes in far right column
#-------------------------------

def end_load_gainphase(fn):

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'


    with open(path + fn) as f:
        lines = f.readlines()


    f = lines[0].split()  #freq in Hz
    p = lines[1].split()  #phase in deg
    g = lines[2].split()  #gain in dB
    f = [float(i) for i in f]
    p = [float(i) for i in p]
    g = [float(i) for i in g]

    #change to radians
    prad = [np.deg2rad(i) for i in p]

    #change gain from dB to linear scale for calculation of transfer function
    #From Steve Martin email on Nov 7, 2022: 
    #Gain=10^(0.05 * (opchan+gainoffset))
    offset = 0.
    Hmag = [10**(0.05*i + offset) for i in g]

    fig, axs = plt.subplots(3)
    axs[0].plot(f,g)
    axs[1].plot(f,Hmag)
    axs[2].plot(f,p)
    #axs[2].plot(f,prad)
    axs[0].set_title('gain/phase; \n fn='+ fn)
    axs[0].set_xscale('log')
    axs[1].set_xscale('log')
    axs[2].set_xscale('log')
    axs[0].set_yscale('linear')
    axs[1].set_ylim(0,20)
    axs[0].set(ylabel='gain(dB)',xlabel='freq(kHz)')
    axs[1].set(ylabel='gain(linear)',xlabel='freq(kHz)')
    axs[2].set(ylabel='phase(deg)',xlabel='freq(kHz)')
    #axs[:].set_xlim(-40,10)
    plt.show()


    return prad, Hmag, f





#-----------------------------------------------
#Load gain/phase curves for desired channel
#-----------------------------------------------
 
prad, Hmag, f = end_load_gainphase("Endurance_Analog 1_VLF12D_6-30000-100.txt")


#Chop off values above 15 kHz, the Nyquist freq. The calibration data goes out to 30 kHz. 
#Not doing this can mess up the interpolation b/t lowres and highres data

f = f[:91]
prad = prad[:91]
Hmag = Hmag[:91]



#-------------------------------
#Load Endurance data for testing
#-------------------------------

vlf = end.efield_vlf()
tvlf = vlf.tvlf
vlf12 = vlf.dvlf12_mvm

fs = vlf.samplerate


#-------------------------------
#FFT Endurance VLF data
#-------------------------------

vr = [-80,-60]
#vr = [-0.1,0.1]
xr = [100,500]
yr = [0,12000]

#FFT to get power spectral density (V^2/Hz). 
#Mode defaults to PSD (V^2/Hz). 
#Complex consists of amplitude (where magnitude = abs(complex))
fspec, tspec, powerm = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
fspec, tspec, powerc = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
#fspec, tspec, powera = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='angle')
#ft, tt, pt = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True)

#-------------------------------------------------------------
#FFT entire Endurance waveform to apply transfer function
#-------------------------------------------------------------


from scipy.fft import rfft, irfft


vlf12FFT = rfft(vlf12)
power12FFT = vlf12FFT**2



N = len(vlf12FFT)
n = np.arange(N)
T = N/fs
freq = n/T/2 

fig, axs = plt.subplots(3)
axs[0].plot(freq,np.abs(power12FFT))
axs[1].plot(freq,np.real(power12FFT))
axs[2].plot(freq,np.imag(power12FFT))
for i in range(3): axs[i].set_xscale('linear')
for i in range(3): axs[i].set_xlim(0,100)



print('here')

#------------------------------------------------------------------
#Transfer function H(w) = |H|*exp(i*theta)
#------------------------------------------------------------------


#Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(f,Hmag,kind='cubic', bounds_error=False)
Hmag2 = interp(freq)
interp2 = interp1d(f,prad,kind='cubic', bounds_error=False)
prad2 = interp2(freq)

#********
#NoTE: Two issues to deal with here before we define the transfer function:
# 1) the gains at the lowest freqs (<10 Hz) are artificial and come from the calibration procedure. 
#   If I leave these as they are (gains < 1) then this ends up hugely boosting calibrated data resulting in a 10 Hz continuous signal that isn't real. 
#   I need to just get rid of this. Set any gain <30 Hz to the value of the gain at 100 Hz (~17.3) 
# 2) The interpolated phase has an over/undershoot at 1 kHz. This is due to interpolating over a sudden change.
#   Need to unwrap values, interpolate, then rewrap
#********

#Remove artifical gains below about 30 Hz. 
badgain_ind = list(map(tuple, np.where(freq < 30)))
Hmag2[badgain_ind] = 17.6



fig, axs = plt.subplots(2)
axs[0].plot(f, Hmag, freq, Hmag2)
axs[1].plot(f, prad, freq, prad2)
for i in range(2): axs[i].set_xscale('log')
for i in range(2): axs[i].set_xlim(1,20000)
axs[1].set_ylim(-4,4)


H = [Hmag2[i] * np.exp(1j*prad2[i]) for i in range(len(prad2))]

fig2, axs2 = plt.subplots(3)
axs2[0].plot(freq,np.abs(H))
axs2[1].plot(freq,np.real(H))
axs2[2].plot(freq,np.imag(H))
for i in range(3): axs2[i].set_xscale('log')




#Apply transfer function to FFT'd data (amplitude, not power)
vlf12FFTc = [vlf12FFT[i]/H[i] for i in range(len(H))]

freq2 = list(freq)




v1 = np.abs(vlf12FFT)
v2 = np.abs(vlf12FFTc)
rat = [v1[i]/v2[i] for i in range(len(v1))]

fig3,axs3 = plt.subplots(4)
axs3[0].plot(freq2,np.abs(vlf12FFT))
axs3[1].plot(freq2,np.abs(vlf12FFTc))
axs3[2].plot(freq2,rat)
axs3[2].set_yscale('linear')
axs3[3].plot(freq2,np.abs(H))
axs3[3].set_yscale('linear')
for i in range(4): axs3[i].set_xlim(0,50)




#Before I inverse transform I need to get rid of NaN values at the beginning of new data arrays
#good = list(map(tuple, np.where(np.isfinite(vlf12FFTc))))
bad = list(map(tuple, np.where(np.isnan(vlf12FFTc))))

vlf12FFTc = np.asarray(vlf12FFTc)
vlf12FFTc[bad] = 0




#Inverse FFT to get back to corrected waveform
wf = irfft(vlf12FFT)
wf_corr = irfft(vlf12FFTc)



fig,axs = plt.subplots(3)
axs[0].plot(tvlf,vlf12)
axs[1].plot(tvlf,wf)
for i in range(3): axs[i].set_ylim([-0.5,0.5])
axs[2].plot(tvlf,wf_corr)
axs[2].set_ylim(-0.02,0.02)
for i in range(3): axs[i].set_xlim([393.1+0.0625,393.1+0.065])


#Phase shift test using normalized gain versions.
goomax = np.max(wf)
wfnorm = wf/goomax
goomax = np.max(wf_corr)
wf_corrnorm = wf_corr/goomax

plt.plot(tvlf,wfnorm, tvlf,wf_corrnorm)
plt.xlim([393.1+0.0625,393.1+0.065])
plt.ylim(-0.02,0.02)







"""
#HF channels 
Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt
Endurance_Analog 1_HF12_1000-20000000-100.txt
Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt
Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt


#VLF digitial files
Endurance_Analog 1_VLF12D_6-30000-100.txt
Endurance_Analog 1_VLF13D_6-30000-100.txt
Endurance_Analog 1_VLF41D_6-30000-100.txt
Endurance_Analog 1_VLF24D_6-30000-100.txt
Endurance_Analog 1_VLF42D_6-30000-100.txt
Endurance_Analog 1_VLF32D_6-30000-100.txt
Endurance_Analog 1_VLF34D_6-30000-100.txt


Endurance_Analog 1_VLF12A_6-100000-100.txt

Endurance_Analog 1_V42D_10-10000-100.txt
Endurance_Analog 1_V41D_10-10000-100.txt
Endurance_Analog 1_V34D_10-10000-100.txt
Endurance_Analog 1_V34A_10-10000-100.txt
Endurance_Analog 1_V32D_10-10000-100.txt
Endurance_Analog 1_V24D_10-10000-100.txt
Endurance_Analog 1_V13D_10-10000-100.txt
Endurance_Analog 1_V12D_10-10000-100.txt
Endurance_Analog 1_V12A_10-10000-100.txt
Endurance_Analog 1_V4SD_10-10000-100.txt
Endurance_Analog 1_V4SA_10-10000-100.txt
Endurance_Analog 1_V3SD_10-10000-100.txt
Endurance_Analog 1_V3SA_10-10000-100.txt
Endurance_Analog 1_V2SD_10-10000-100.txt
Endurance_Analog 1_V2SA_10-10000-100.txt
Endurance_Analog 1_V1SD_10-10000-100.txt
Endurance_Analog 1_V1SA_10-10000-100.txt

"""


