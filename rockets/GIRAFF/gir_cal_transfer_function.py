"""
Calculate transfer function for GIRAFF channels. This is used to gain/phase correct the data.

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians
(3) Apply transfer function correction as: fft_data_corrected = fft_data/H


(see https://resources.pcb.cadence.com/blog/2021-understanding-a-circuit-transfer-function-from-a-bode-plot)
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
from math import remainder
import pickle


#-------------------------------------------------------
#Select data channel for calibration 
#-------------------------------------------------------

pld = '381'
#ch = 'V34D'
ch = 'VLF34D'
#ch = 'V4SD'
#ch = 'VLF41D'
v = GFL(pld, ch)

v.plot_gainphase()

vdat = v.load_data()
wavedat = vdat[0]
tdat = vdat[1]


#sample rate
fs = v.chnspecs['fs']

#FOR VLF CHANNELS ONLY: Max freq of the downsampled VLF channels (e.g. VLF12D vs VLF12DF)
fmax_vlf_downsampled = 50000



#-----------------------------------------------------------------------------------------
#Load gain/phase data for selected channel. 
#-----------------------------------------------------------------------------------------


phase_lowres = np.asarray(v.phase)
gain_lowres = np.asarray(v.gain)
freq_lowres = np.asarray(v.freq_gainphase)

#For the normal VLF data loaded, I've desampled to 500 kHz (otherwise too much data to load)


nyquist = fs/2


#-----------------------------------------------------------------------------------------
#Chop off calibration values above 7/8 of the Nyquist freq of data.
#Not doing this can mess up the interpolations below...
#-----------------------------------------------------------------------------------------

good = np.where(freq_lowres <= nyquist * (7/8))
phase_lowres = phase_lowres[np.where(freq_lowres <= nyquist * (7/8))]
gain_lowres = gain_lowres[np.where(freq_lowres <= nyquist * (7/8))]
freq_lowres = freq_lowres[np.where(freq_lowres <= nyquist * (7/8))]



#-----------------------------------------------------------------------------------------
#At this point we have an issue with the DC data. The interpolated frequency array (freq) does go down to 0 Hz, however
#the gain and phase have NaN values at < 10 Hz (or similar) b/c the AC testing setup can't go to 0 Hz in freq. 
#Therefore we need to assign these NaN values a finite value based on the values near ~10 Hz. 
#NOTE: need to do multiple values <10 Hz or else the interpolation to higher resolution (later) gets wonky
#-----------------------------------------------------------------------------------------

if v.type == 'DC' or v.type == 'skins':

    goo = np.where(freq_lowres >= 20)    
    g_replace = gain_lowres[goo[0][0]]
    p_replace = phase_lowres[goo[0][0]]


    #define additional value at zero freq
    gain_lowres = np.append([g_replace]*10, gain_lowres)
    phase_lowres = np.append([p_replace]*10, phase_lowres)
    freq_lowres = np.append(list(range(10)), freq_lowres)


    fmin_trust = 0  #used later
    fmax_trust = (7/8)*v.chnspecs['fs'] / 2

if v.type == 'VLF':

    #--Change any gain values < 1 to 1 so that low freqs are not amplified. We will just consider these values to be outside
    # of the VLF channel range. Not doing so tends to produce large power at very low and high freqs.     
    gain_lowres[np.where(gain_lowres < 1)] = 1



    #--Replace VLF cal test gain/phase values < fmin_trust (Hz) b/c I don't trust these. 
    #--The signal/noise ratio is very high in the cal tests, and these gain/phase values can be all over the place. 
    fmin_trust = 2*v.chnspecs['hpf']
    fmax_trust = fmax_vlf_downsampled    #For the normal VLF file to load I downsample the data to freqs <= 50k. 

    #--define replacement values 
    whgoo = np.where(freq_lowres <= fmin_trust)
    g_replace = gain_lowres[whgoo[0][-1]]
    p_replace = phase_lowres[whgoo[0][-1]]
    gain_lowres[whgoo] = g_replace
    phase_lowres[whgoo] = p_replace

    #define additional value at zero freq
    gain_lowres = np.append(g_replace, gain_lowres)
    phase_lowres = np.append(p_replace, phase_lowres)
    freq_lowres = np.append(0, freq_lowres)




#-------------------------------------------------------------
#FFT GIRAFF waveform data to apply transfer function
#-------------------------------------------------------------


wavedatFFT = rfft(wavedat)


#--Construct a list of higher resolution frequencies (than those from the gain/phase files) that go up to half the Nyquist
N = len(wavedatFFT)
n = list(range(N))
T = N/fs
freq_hires = np.asarray([i/T/2 for i in n])


#------------------------------------------------------------------
#Transfer function = |gain|*exp(i*theta)
#------------------------------------------------------------------


#--To avoid interpolation issues due to phase jumps (-180 to 180), need to unwrap values, interpolate, then rewrap
phase_lowres_unwrapped = np.unwrap(phase_lowres)


#--Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(freq_lowres,gain_lowres,kind='cubic', bounds_error=False)
gain_hires = interp(freq_hires)
interp2 = interp1d(freq_lowres,phase_lowres_unwrapped,kind='cubic', bounds_error=False)
phase_hires_unwrapped = interp2(freq_hires)




#--rewrap angles from -2pi to 2pi 
phase_hires = np.asarray([remainder(phase_hires_unwrapped[i], 2*np.pi) for i in range(len(phase_hires_unwrapped))])


plt.plot(freq_lowres, phase_lowres, '*', freq_hires,phase_hires_unwrapped, freq_hires, phase_hires)
plt.axvline(x=nyquist, color='r', linestyle='--')
plt.xscale('log')
plt.xlim(1,1000000)
plt.ylim(-4,4)
print('check angle wrapping')
plt.close()



fig, axs = plt.subplots(2)
axs[0].plot(freq_lowres, gain_lowres, '*', freq_hires, gain_hires)
axs[1].plot(freq_lowres, phase_lowres,'*', freq_hires, phase_hires)
axs[0].axvline(nyquist, color='r', linestyle='--')
axs[1].axvline(nyquist, color='r', linestyle='--')
for i in range(2): axs[i].set_xscale('log')
for i in range(2): axs[i].set_xlim(1,1000000)
axs[1].set_ylim(-4,4)

print('Check modified gain/phase curves')
plt.close(fig)


#---------------------------------------------------
#Define the transfer function.
#Two versions:
#1) gain + phase
#2) unity gain + phase 

#--note that the files Steve sent me for Endurance already have the DC gain correction applied
#--for both the DC and VLF channels (Paulo must have done this using his yellow-highlighted fit values on the MEB). So when applying the transfer function
#--to these data use the unity gain version. 
#---------------------------------------------------

transfer_func = gain_hires * np.exp(1j*phase_hires)
goodf = np.where((freq_hires > fmin_trust) & (freq_hires < fmax_trust))[0]
maxv = np.nanmax(gain_hires[goodf])
transfer_func_unitygain = (gain_hires/maxv) * np.exp(1j*phase_hires)


#Remove any NaN values in the transfer function. These will mess up the inverse FFT. 
bad = np.where(np.isnan(transfer_func))
if len(bad[0]) != 0:
    transfer_func[bad] = transfer_func[bad[0][0]-1]
    transfer_func_unitygain[bad] = transfer_func_unitygain[bad[0][0]-1]



fig2, axs2 = plt.subplots(3)
axs2[0].plot(freq_hires,np.abs(transfer_func))
axs2[1].plot(freq_hires,np.real(transfer_func))
axs2[2].plot(freq_hires,np.imag(transfer_func))
for i in range(3): axs2[i].set_xscale('log')
for i in range(3): axs2[i].set_xlim(1,1000000)
for i in range(3): 
    axs2[i].axvline(nyquist, color='r', linestyle='--')
for i in range(3): axs2[i].set_xscale('log')

#For the downsampled VLF channels, indicate that the data are truncated above 50 kHz
if ch == 'VLF12D' or ch == 'VLF34D' or ch == 'VLF13D' or ch == 'VLF32D' or ch == 'VLF24D' or ch == 'VLF41D':
    for i in range(3): 
        axs2[i].axvline(fmax_vlf_downsampled, color='r')
print('Check transfer function')





#--Bode plots
fig2, axs2 = plt.subplots(3)
axs2[0].plot(freq_hires,np.abs(transfer_func_unitygain))
axs2[1].plot(freq_hires,np.real(transfer_func_unitygain))
axs2[2].plot(freq_hires,np.imag(transfer_func_unitygain))
for i in range(3): axs2[i].set_xscale('log')
for i in range(3): axs2[i].set_xlim(1,1000000)
for i in range(3): axs2[i].set_xlabel('freq(Hz)')
for i in range(3): 
    axs2[i].axvline(nyquist, color='r', linestyle='--')
axs2[0].set_title('Bode plots \n ' + v.chnspecs["gainphase_file"][:-4])
axs2[0].set_ylabel('Magnitude |H(jw)|')
axs2[1].set_ylabel('Re(H(jw))')
axs2[2].set_ylabel('Phase Im(H(jw))')

#For the downsampled VLF channels, indicate that the data are truncated above 50 kHz
if ch == 'VLF12D' or ch == 'VLF34D' or ch == 'VLF13D' or ch == 'VLF32D' or ch == 'VLF24D' or ch == 'VLF41D':
    for i in range(3): 
        axs2[i].axvline(fmax_vlf_downsampled, color='r')
print('Check unity gain transfer function')



#---------------------------------------------------
#Apply transfer function to FFT'd data (amplitude, not power)
#---------------------------------------------------

wavedatFFTc = wavedatFFT/transfer_func_unitygain



#---------------------------------------------------
#Inverse FFT to get back to corrected waveform
#---------------------------------------------------

wf = irfft(wavedatFFT, n=len(tdat))
wf_corr = irfft(wavedatFFTc, n=len(tdat))



#----------------------------------------------------
#Save corrected waveform data
#----------------------------------------------------

dict_fin = {'tvals':tdat, 'wf':wf_corr}
pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/' + v.chnspecs['folder'] + '/'


#----------------------------
#Set output file name   Giraff_381_Analog 2_VLF12D_100-100000-50_gainphase_corrected.pkl
#----------------------------
if pld == '381':
    agoo = '2_'
else: 
    agoo = '1_'


t1 = 'Giraff_'+pld+'_'+'Analog '+agoo+ch+'_'
flow = str(int(fmin_trust)) + '-'


goo = v.chnspecs["gainphase_file"].find('-')
stmp = v.chnspecs["gainphase_file"][goo:]
goo = stmp.find('-',1)
tnsteps = stmp[goo+1:-4]+'steps'
if ch == 'VLF12D' or ch == 'VLF34D' or ch == 'VLF13D' or ch == 'VLF32D' or ch == 'VLF24D' or ch == 'VLF41D':
    fhigh = str(100000) + 'Hz_'
else:
    fhigh = stmp[1:goo] + 'Hz_'

t2 = '_gainphase_corrected'


fn = t1 + flow + fhigh + tnsteps + t2


pickle.dump(dict_fin, open(pathoutput + fn + '.pkl','wb'))
#wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))



