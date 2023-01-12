"""
Calibrate the AC data to exact mV/m based on the DC data in the overlap region around 10 Hz.

From Paulo: The DC data have a very accurate calibration from counts to mV/m. The AC data are not accurate.
Here we use the fact that the two channels overlap in frequency around 10 Hz to better calibrate the AC data.

"""


import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
import end_load_data as end
from scipy.interpolate import interp1d
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
import filter_wave_frequency



path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'


#Load Endurance data for testing
vlf = end.efield_vlf()
tvlf = vlf.tvlf
vlf12 = vlf.dvlf12_mvm

fs = vlf.samplerate

edc = end.efield_dc()
edcE = edc['dv12_mvm']
edcEc = edc['dv12_raw']
tedc = edc['times']
fsDC = edc.samplerate

plt.plot(edcEc)
#dv{12/34/13/32/24/41}_mvm = mV/m
#dv{12/34/13/32/24/41}_volts = Volts
#dv{12/34/13/32/24/41}_raw = Counts
#dv{12/34/13/32/24/41}_rawnormal = [-1,1]
#times = Seconds since T-0 for all channels'



#FFT Endurance VLF data
vr = [-80,-60]
#vr = [-0.1,0.1]
xr = [100,500]
yr = [0,12000]

#FFT to get power spectral density (V^2/Hz). 
#fig, ax = plt.subplots(2)


fspec, tspec, powerm = signal.spectrogram(vlf12, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
fspecDC, tspecDC, powermDC = signal.spectrogram(edcE, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
ps.plot_spectrogram(tspec,fspec,powerm,vr=[-80,-10],yr=[1,1000],xr=[150,170], yscale='log')
ps.plot_spectrogram(tspecDC,fspecDC,powermDC,vr=[-80,-10],yr=[1,1000],xr=[150,170], yscale='log')



#Filter wave data to b/t 100-1000 Hz for DC to AC comparison

lowcut = 1000
highcut = 2000
vlf12bp = filter_wave_frequency.butter_bandpass_filter(vlf12, lowcut, highcut,fs,order=12)
edcEbp = filter_wave_frequency.butter_bandpass_filter(edcE, lowcut, highcut,fsDC,order=12)

fspecbp, tspecbp, powermbp = signal.spectrogram(vlf12bp, fs, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
fspecDCbp, tspecDCbp, powermDCbp = signal.spectrogram(edcEbp, fsDC, nperseg=1024,noverlap=1024/2,window='hann',return_onesided=True,mode='magnitude')
ps.plot_spectrogram(tspecbp,fspecbp,powermbp,vr=[-80,-10],yr=[1,1000],xr=[150,170], yscale='log')
ps.plot_spectrogram(tspecDCbp,fspecDCbp,powermDCbp,vr=[-60,0],yr=[1,1000],xr=[150,170], yscale='log')


plt.plot(tvlf,vlf12bp,tedc,edcEbp)
#plt.xlim(150,170)
plt.xlim(165.50,165.54)
plt.ylim(-0.02,0.02)
plt.ylabel("Efield\n[mV/m]")
plt.xlabel("time (sec)")
plt.title("Endurance Efield comparison (1000-2000 Hz)\nOrange=EDC12; Blue=VLF12\nSteve's Calibrations")


#Test
f, t, p = signal.spectrogram(wfbp, fs, nperseg=16384,noverlap=16384/2,window='hann')
ps.plot_spectrogram(t,f,p,vr=[-100,-40],yr=yrspec, xr=trspec, yscale='linear')




#fspec, tspec, powera = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='angle')
#ft, tt, pt = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True)

print("here")





#Test calibrate b/t counts to volts to see if I'm doing this correctly. 
#Compare with Steve Martin's data
#Vout = G * Ncounts * (2.5/bit_depth)
#G = Vin/Vadc 

#VLF12 
#G = 19.7
#bit_depth = 2**18

#Steve's calibrations for EDC:
    #Counts -> Volts: 
    #dv12_volts = A * dv12_rawnormal + B, 
    #Volts -> mV/m = dv12_mvm=1000.0*dv12_volts/bl12 

#Steve's calibrations for VLF12
    #Counts -> Volts: 
    #dvlf12_volts = A * dvlf12_normal + B 
    #A = 1.265e-1; B=  5.5604e-5,  boomlength=3.212'


#EDC 12 
#G = 2 
#From Paulo's "conversion equation" parameters a, b on the AC cal document
G = 1.2555
B = 0.0011  #offset
bit_depth = 2**18  #so range goes from +/-131072

edcE_compare = edc['dv12_mvm']
edcEv_compare = edc['dv12_volts']
edcEc = edc['dv12_raw']
edcEcn = edc['dv12_rawnormal']
edcEv = edc['dv12_volts']

#counts to volts
#Aaron method
#edcEv_test = G * edcEc * (2.5/(bit_depth/2)) + B
#Steve method
edcEv_test = G * edcEcn

fig, axs = plt.subplots(3)

axs[0].plot(edc['times'],edcEc)
axs[1].plot(edc['times'],edcEv_test)
axs[2].plot(edc['times'],edcEv_compare)

plt.plot(edc['times'], edcE_test, edc['times'],edcE_compare)

"""

#---------------------------
#Load AC calibration data
#---------------------------

fn = "Endurance_Analog 1_VLF12D_6-30000-100.txt"
with open(path + fn) as f:
    lines = f.readlines()

f = lines[0].split()  #freq in Hz
p = lines[1].split()  #phase in deg
g = lines[2].split()  #gain in dB
f_ac = [float(i) for i in f]
p_ac = [float(i) for i in p]
g_ac = [float(i) for i in g]
#change to radians
prad = [np.deg2rad(i) for i in p_ac]
#change gain from dB to linear scaler for calculation of transfer function
#From Steve Martin email on Nov 7, 2022: 
#Gain=10^(0.05 * (opchan+gainoffset))
offset = 0.
Hmag_ac = [10**(0.05*i + offset) for i in g_ac]


#---------------------------
#Load DC calibration data
#---------------------------

fn = "Endurance_Analog 1_V12D_10-10000-100.txt"
with open(path + fn) as f:
    lines = f.readlines()

f = lines[0].split()  #freq in Hz
p = lines[1].split()  #phase in deg
g = lines[2].split()  #gain in dB
f_dc = [float(i) for i in f]
p_dc = [float(i) for i in p]
g_dc = [float(i) for i in g]
#change to radians
prad = [np.deg2rad(i) for i in p_dc]
offset = 0.
Hmag_dc = [10**(0.05*i + offset) for i in g_dc]




fig, axs = plt.subplots(3)
axs[0].plot(f_dc,g_dc,f_ac,g_ac)
axs[1].plot(f_dc,Hmag_dc,f_ac,Hmag_ac)
axs[2].plot(f_dc,p_dc,f_ac,p_ac)
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
for i in range(3): axs[i].set_xlim(0,100)
#axs[:].set_xlim(-40,10)
plt.show()

"""












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

