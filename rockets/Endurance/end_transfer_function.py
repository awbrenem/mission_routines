"""
Calculate transfer function for Endurance channels. This is used to gain/phase correct the data.

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians
(3) Apply transfer function correction as: fft_data_corrected = fft_data/H


(see https://resources.pcb.cadence.com/blog/2021-understanding-a-circuit-transfer-function-from-a-bode-plot)
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
from math import remainder


#-------------------------------------------------------
#Select data channel for calibration 
#-------------------------------------------------------

ch = 'V13D'
#ch = 'V4SD'
#ch = 'VLF41D'
v = EFL(ch)

v.plot_gainphase()

vdat = v.load_data()
wavedat = vdat[0]
tdat = vdat[1]


#sample rate
fs = v.chnspecs['fs']



#-----------------------------------------------------------------------------------------
#Load gain/phase data for selected channel. 
#---NOTE: channels on the MEB with an identified negative polarity have had their phases flipped 
#---in the following load routine. This needs to be done.
#-----------------------------------------------------------------------------------------


phase_lowres = np.asarray(v.phase)
gain_lowres = np.asarray(v.gain)
freq_lowres = np.asarray(v.freq_gainphase)

nyquist = fs/2


#-----------------------------------------------------------------------------------------
#Chop off calibration values above the Nyquist freq of data.
#Not doing this can mess up the interpolations below...
#-----------------------------------------------------------------------------------------

good = np.where(freq_lowres <= nyquist)
phase_lowres = phase_lowres[np.where(freq_lowres <= nyquist)]
gain_lowres = gain_lowres[np.where(freq_lowres <= nyquist)]
freq_lowres = freq_lowres[np.where(freq_lowres <= nyquist)]



#-----------------------------------------------------------------------------------------
#At this point we have an issue with the DC data. The interpolated frequency array (freq) does go down to 0 Hz, however
#the gain and phase have NaN values at < 10 Hz (or similar) b/c the AC testing setup can't go to 0 Hz in freq. 
#Therefore we need to assign these NaN values a finite value based on the values near ~10 Hz. 
#NOTE: need to do multiple values <10 Hz or else the interpolation to higher resolution (later) gets wonky
#-----------------------------------------------------------------------------------------

if v.type == 'DC' or v.type == 'skins':

    #NOTE: this method doesn't work if there's a phase flip at low freqs!!!
        ##--Mean value b/t 20-100 Hz - to be used to populate <=10 Hz values
        #g_replace = np.mean([gain_lowres[i] for i in range(len(freq_lowres)) if (freq_lowres[i] > 20) & (freq_lowres[i] < 100)])
        #p_replace = np.mean([phase_lowres[i] for i in range(len(freq_lowres)) if (freq_lowres[i] > 20) & (freq_lowres[i] < 100)])

    goo = np.where(freq_lowres >= 20)    
    g_replace = gain_lowres[goo[0][0]]
    p_replace = phase_lowres[goo[0][0]]


    #define additional value at zero freq
    gain_lowres = np.append([g_replace]*10, gain_lowres)
    phase_lowres = np.append([p_replace]*10, phase_lowres)
    freq_lowres = np.append(list(range(10)), freq_lowres)


    fmin_trust = 0  #used later

if v.type == 'VLF':

    #--Change any gain values < 1 to 1 so that low freqs are not amplified.     
    gain_lowres[np.where(gain_lowres < 1)] = 1



    #--Replace VLF cal test gain/phase values < fmin_trust (Hz) b/c I don't trust these. 
    #--The signal/noise ratio is very high in the cal tests, and these gain/phase values can be all over the place. 
    fmin_trust = 30
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
#FFT Endurance waveform data to apply transfer function
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


plt.plot(freq_lowres, phase_lowres, freq_hires,phase_hires_unwrapped, freq_hires, phase_hires)
plt.axvline(x=nyquist, color='r', linestyle='--')
plt.xscale('log')
plt.xlim(1,50000)
plt.ylim(-4,4)
print('check wrapping')
plt.close()



fig, axs = plt.subplots(2)
axs[0].plot(freq_lowres, gain_lowres, freq_hires, gain_hires)
axs[1].plot(freq_lowres, phase_lowres, freq_hires, phase_hires)
axs[0].axvline(nyquist, color='r', linestyle='--')
axs[1].axvline(nyquist, color='r', linestyle='--')
for i in range(2): axs[i].set_xscale('log')
for i in range(2): axs[i].set_xlim(1,30000)
axs[1].set_ylim(-4,4)

print('Check modified gain/phase curves')
fig.clear()
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
maxv = np.nanmax(gain_hires)
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
for i in range(3): axs2[i].set_xlim(1,20000)
for i in range(3): 
    axs2[i].axvline(nyquist, color='r', linestyle='--')
for i in range(3): axs2[i].set_xscale('log')
print('Check transfer function')





#--Bode plots
fig2, axs2 = plt.subplots(3)
axs2[0].plot(freq_hires,np.abs(transfer_func_unitygain))
axs2[1].plot(freq_hires,np.real(transfer_func_unitygain))
axs2[2].plot(freq_hires,np.imag(transfer_func_unitygain))
for i in range(3): axs2[i].set_xscale('log')
for i in range(3): axs2[i].set_xlim(1,20000)
for i in range(3): axs2[i].set_xlabel('freq(Hz)')
for i in range(3): 
    axs2[i].axvline(nyquist, color='r', linestyle='--')
axs2[0].set_title('Bode plots \n ' + v.chnspecs["gainphase_file"][:-4])
axs2[0].set_ylabel('Magnitude |H(jw)|')
axs2[1].set_ylabel('Re(H(jw))')
axs2[2].set_ylabel('Phase Im(H(jw))')

print('Check unity gain transfer function')



#---------------------------------------------------
#Apply transfer function to FFT'd data (amplitude, not power)
#---------------------------------------------------

wavedatFFTc = wavedatFFT/transfer_func_unitygain



"""
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

fig3.clear()
plt.close(fig3)
"""


#---------------------------------------------------
#Inverse FFT to get back to corrected waveform
#---------------------------------------------------

wf = irfft(wavedatFFT, n=len(tdat))
wf_corr = irfft(wavedatFFTc, n=len(tdat))



fig,axs = plt.subplots(3, clear=True)

axs[0].plot(tdat,wf)
axs[1].plot(tdat,wf_corr)


#--DC signal test
for i in range(3): axs[i].set_xlim(100,180)
for i in range(3): axs[i].set_xlim(110,115)
for i in range(3): axs[i].set_ylim(-1,1)
#--Small amp test
for i in range(3): axs[i].set_xlim(535,536.2)
for i in range(3): axs[i].set_ylim(-0.25,0.25)
#--Maneuver test
for i in range(3): axs[i].set_xlim(524,525)
for i in range(3): axs[i].set_ylim(-1,1)



for i in range(3): axs[i].set_xlim(656,656.06)
for i in range(3): axs[i].set_xlim(656.02,656.03)
for i in range(3): axs[i].set_xlim(450,460)
for i in range(3): axs[i].set_ylim(-0.45,0.45)
axs[1].set_ylim(-0.25,0.25)
axs[2].set_ylim(-0.25,0.25)


fig.clear()
plt.close(fig)
 




#----------------------------------------------------
#Save corrected waveform data
#----------------------------------------------------

import pickle

dict_fin = {'tvals':tdat, 'wf':wf_corr}
pathoutput = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/' + v.chnspecs['folder'] + '/'

fnsav = v.chnspecs["gainphase_file"][:-4] + '_gainphase_corrected'
pickle.dump(dict_fin, open(pathoutput + fnsav + '.pkl','wb'))
#wf_corr_load = pickle.load(open(pathoutput + fnsav + ".pkl", 'rb'))



