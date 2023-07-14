
from scipy.signal import butter, sosfilt, sosfreqz
#sos bandpass
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    #b, a = butter(order, [low, high], btype='band')
    #This method (using sosfiltfilt) preserves signal phase!! (tested and verified)(seehttps://gist.github.com/junzis/e06eca03747fc194e322)
    sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    return sos
 
 
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sosfilt(sos, data)
    return y

import scipy.io
from matplotlib import pyplot as plt
import numpy as np

rs = scipy.io.readsav(r"C:\Users\Esther Kwon\Documents\Breneman docs\v_amp_data(valid)\35040_LFDSP_S1_VLF12_mvm.sav")
print(rs.units)
s12 = scipy.io.readsav(r"C:\Users\Esther Kwon\Documents\Breneman docs\v_amp_data(valid)\35040_LFDSP_S1_VLF12_mvm.sav")

rs = scipy.io.readsav(r"C:\Users\Esther Kwon\Documents\Breneman docs\v_amp_data(valid)\35040_LFDSP_S1_VLF34_mvm.sav")
print(rs.units)
s34 = scipy.io.readsav(r"C:\Users\Esther Kwon\Documents\Breneman docs\v_amp_data(valid)\35040_LFDSP_S1_VLF34_mvm.sav")

from scipy.io import wavfile
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col

#load vlf data
times34 = s34.tvlf34
times12 = s12.tvlf12
print(np.shape(times12), np.shape(times34))

amps34 = s34.dvlf34
amps12 = s12.dvlf12
print(np.shape(amps12), np.shape(amps34))

#sampling freq
sr_lf34 = [1/(times34[i+1]-times34[i]) for i in range(times34.size-1)]
sr_lf12 = [1/(times12[i+1]-times12[i]) for i in range(times12.size-1)]
print(np.shape(sr_lf12), np.shape(sr_lf34))

fsVLF34 = np.mean(sr_lf34)
fsVLF12 = np.mean(sr_lf12)
print(fsVLF12, fsVLF34)

#add image of non-bandpassed

#bandpassed to remove low and high freq power that we're not interested in
amps12 = butter_bandpass_filter(amps12, 200, 3000, fsVLF12, order=10)
amps34 = butter_bandpass_filter(amps34, 200, 3000, fsVLF34, order=10)



import sys
sys.path.append(r'C:\Users\Esther Kwon\PycharmProjects\estherpython')
import plot_spectrogram as ps
specfreqs12, spectimes12, power12 = signal.spectrogram(amps12, fsVLF12, nperseg=512,noverlap=None,window='hann') #, return_onesided=1)
specfreqs34, spectimes34, power34 = signal.spectrogram(amps34, fsVLF34, nperseg=512,noverlap=None,window='hann') #, return_onesided=1)
ps.plot_spectrogram(spectimes12,specfreqs12,power12,vr=[-80,-40],yr=[100,10000],xr=[100,800], yscale='log')
ps.plot_spectrogram(spectimes34,specfreqs34,power34,vr=[-40,-40],yr=[100,10000],xr=[100,800], yscale='log')




spectimes34 = spectimes34[0 : 41452]
power34 = power34[:, 0 : 41452]
print(np.shape(power12))
print(np.shape(power34))

#BBELF1(broadband extremely low freq wave?)

#ps.plot_spectrogram(spectimes12,specfreqs12,power12,vr=[-80,-20],yr=[300,1000],xr=[460,500], yscale='linear')
#ps.plot_spectrogram(spectimes34,specfreqs34,power34,vr=[-80,-20],yr=[300,1000],xr=[460,500], yscale='linear')

ps.plot_spectrogram(spectimes12,specfreqs12,power12,vr=[-80, -20],yr=[1000,4000],xr=[650,700], yscale='log')
ps.plot_spectrogram(spectimes34,specfreqs34,power34,vr=[-80, -20],yr=[1000,4000],xr=[650,700], yscale='log')


flo = 1300
fhi = 2500
order = 10
tmin = 670
tmax = 690

"""BBELF1(broadband extremely low freq wave?) zoom in parameters
polarization varies from circular to linear over t range "onion"
flo = 600
fhi = 700
order = 10
tmin = 481.055
tmax = 481.055+.018
"""

"""BBELF1(broadband extremely low freq wave?) zoom out parameters
#flo = 350
#fhi = 800
#tmin = 481.04
#tmax = 481.04+.018*10
"""

"""cleanps during collection
flo = 250
fhi = 310
order = 10
tmin = 382.68
tmax = 382.72

"""

"""phase shift of analog @ beginning
flo = 900
fhi = 1100
order = 10
tmin = 130.66
tmax = 130.80
"""

"""og
flo = 400
fhi = 700
tmin = 525.8570
tmax = 525.8640
order = 10
"""

amps12bp = butter_bandpass_filter(amps12, flo, fhi, fsVLF12, order=order)
amps34bp = butter_bandpass_filter(amps34, flo, fhi, fsVLF34, order=order)
print(np.shape(amps12bp))
print(np.shape(amps12))

nperseg = 1024*4
specfreqs12, spectimes12, power12bp = signal.spectrogram(amps12bp, fsVLF12, nperseg=nperseg, noverlap=nperseg//2, window='hann') #, return_onesided=1)
ps.plot_spectrogram(spectimes12,specfreqs12,power12bp,vr=[-80,-20],yr=[200,400],xr=[tmin,tmax], yscale='linear')

specfreqs34, spectimes34, power34bp = signal.spectrogram(amps34bp, fsVLF34, nperseg=nperseg, noverlap=nperseg//2, window='hann') #, return_onesided=1)
ps.plot_spectrogram(spectimes34,specfreqs34,power34bp,vr=[-80,-20],yr=[200,400],xr=[tmin,tmax], yscale='linear')

a = times12

indices = [idx for idx,a in enumerate(a) if (a > tmin) & (a < tmax)]
times12n = times12[indices]
amps12bpn = amps12bp[indices]
times34n = times34[indices]
amps34bpn = amps34bp[indices]
print(len(indices))


plt.plot(amps12bpn)
plt.plot(amps34bpn)

#plt.xlim(525.843, 525.843+.022)
plt.show()

rgoo = .25

plt.axes().set_aspect('equal')
plt.plot(amps12bpn, amps34bpn)
plt.xlim(-rgoo, rgoo)
plt.ylim(-rgoo, rgoo)

plt.show()
