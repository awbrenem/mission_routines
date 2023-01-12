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

    #change gain from dB to linear scaler for calculation of transfer function
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
    axs[1].set(ylabel='gain',xlabel='freq(kHz)')
    axs[2].set(ylabel='phase(deg)',xlabel='freq(kHz)')
    #axs[:].set_xlim(-40,10)
    plt.show()


    return prad, Hmag, f





#-----------------------------------------------
#Load gain/phase curves for desired channel
#-----------------------------------------------
 
prad, Hmag, f = end_load_gainphase("Endurance_Analog 1_VLF12D_6-30000-100.txt")




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
fspec, tspec, powerc = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerm = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
fspec, tspec, powera = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='angle')
#ft, tt, pt = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True)




#------------------------------------------------------------------
#Transfer function  H(w) = |H|*exp(i*theta)
#------------------------------------------------------------------


#Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(f,Hmag,kind='cubic', bounds_error=False)
Hmag = interp(fspec)
interp2 = interp1d(f,prad,kind='cubic', bounds_error=False)
prad = interp2(fspec)

H = [Hmag[i] * np.exp(1j*prad[i]) for i in range(len(prad))]


"""
fig, axs = plt.subplots(2)
axs[0].plot(f,p)
axs[1].plot(fspec, prad)
axs[0].set_xlim(0,15000)
axs[1].set_xlim(0,15000)
plt.show()
"""

#----------------------------------------------------------------------------
#Compare distribution of uncorrected vs corrected power FOR A SINGLE TEST TIME
#----------------------------------------------------------------------------


#Apply transfer function to the complex FFT data.
#(1) apply it to the complex FFT 
ttime = 2000



#The transfer function correction applies to the POWER, not the amplitude
corrected = powerc[:,ttime]/H
#corrected = powerc[:,ttime]/Htst

#Check to make sure the gain has been applied correctly
ratio = [np.abs(powerc[i,ttime])/np.abs(corrected[i]) for i in range(len(H))]


fig, axs = plt.subplots(2)
axs[0].plot(fspec,np.abs(powerc[:,ttime]),'.',fspec,np.abs(corrected),'x')
axs[1].plot(fspec,ratio)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlim(1,100000)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xlim(1,100000)
axs[0].set_ylim(1e-7,1e-3)
plt.show()


#----------------------------------------------------------------------------
#Gain/phase correct all the spectral data. 
#----------------------------------------------------------------------------

powerc_corr = powerc.copy()
for i in range(len(powerc[0,:])):   #Loop through each time
    powerc_corr[:,i] = powerc[:,i]/H




ptmp = np.abs(powerc_corr)

vr = [-80,-20]
xr = [100,500]
yr = [0,12000]

fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspec, fspec, powerm, vr=vr, yscale='linear', zscale='log', xr=xr, yr=yr, ax=axs[0], show=False)
axs[0].set_title = 'test'
ps.plot_spectrogram(tspec, fspec, ptmp, vr=vr, yscale='linear', zscale='log', xr=xr, yr=yr, ax=axs[1])
axs[1].set_title = 'test2'
plt.show()


plt.plot(fspec,powerm[:,3000])
plt.plot(fspec,powerc[:,3000])
plt.plot(fspec,powerc_corr[:,3000])


#------------------------------------------------------
#invert the corrected spectral data to waveforms
#NOTE: IT SEEMS THAT I NEED THE NEGATIVE FREQUENCIES TO USE IFFT
#------------------------------------------------------


#array with positive followed by negative frequencies (as required by ifft)
fspecc = np.concatenate((fspec,-1*np.flip(fspec)))

#Loop through each time and add the negative freq data to the positive freq data to create array with both
pgoo = []
pgoo_corr = []
for i in range(len(tspec)):
    ftst = np.concatenate((powerc[:,i],np.flip(powerc[:,i])))
    pgoo.append(ftst)
    ftst = np.concatenate((powerc_corr[:,i],np.flip(powerc_corr[:,i])))
    pgoo_corr.append(ftst)

powerc2 = np.transpose(pgoo)
powerc_corr2 = np.transpose(pgoo_corr)


#plt.plot(fspecc,pgoo[:,100])
plt.plot(fspecc,powerc2[:,0])
plt.plot(fspecc,powerc_corr2[:,0])
plt.ylim(-0.01,0.01)
plt.xlim(-1000,1000)


#get rid of NaN values
#ttime = 2000
#ptmp = powerc[5:,ttime]
#ptmp_corr = powerc_corr[5:,ttime]


wf = np.fft.ifft(powerc2[:,10])
wf_corr = np.fft.ifft(powerc_corr2[:,5900])

wf2 = np.ndarray.flatten(wf)
wf_corr2 = np.ndarray.flatten(wf_corr)

wf2 = wf2[:len(wf2)//2]
wf_corr2 = wf_corr2[:len(wf_corr2)//2]





wf = np.fft.ifft(powerc2)
wf_corr = np.fft.ifft(powerc_corr2)

wf2 = np.ndarray.flatten(wf)
wf_corr2 = np.ndarray.flatten(wf_corr)

wf2 = wf2[:len(wf2)//2]
wf_corr2 = wf_corr2[:len(wf_corr2)//2]

#wf = np.fft.ifft(np.abs(powerc[:,ttime]))
plt.plot(np.real(wf2[100:200]))
plt.plot(np.real(wf_corr2[100:200]))

plt.plot(np.abs(wf_corr2[0:len(wf_corr2)//10]))
plt.show()


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


