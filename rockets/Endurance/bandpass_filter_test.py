#Crib sheet I wrote for Esther


import sys 
#sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
#import pandas as pd
import matplotlib.pyplot as plt
#from scipy import signal
import numpy as np
import plot_spectrogram as ps
import filter_wave_frequency as filt



time = np.arange(0,10,0.001)  #seconds, say...
fs = len(time)/10

#I'll define what the waveform period is b/c this is easier to visualize than the frequency 
#on a timescale ranging from 0-10 seconds.
p1 = 0.5  #period in sec
p2 = 4*p1

f1 = 1/p1   # = 2 Hz
f2 = 1/p2   # = 0.5 Hz

a1 = 1
a2 = 0.5

#waveforms are defined with angular frequency (2*np.pi*f), not frequency
waveform1 = a1*(np.sin(2*np.pi*f1*time))
waveform2 = a2*(np.sin(2*np.pi*f2*time))
waveform = waveform1 + waveform2
plt.plot(time, waveform1, time, waveform2, time, waveform)
plt.show()


#Filter out the 0.5 Hz wave (should be nearly identical to waveform1)
f1 = w1/(2*np.pi) # = 2 Hz
f2 = w2/(2*np.pi) # = 0.5 Hz
highpass = filt.butter_bandpass_filter(waveform, 1, 10, fs, order= 10)
plt.plot(time, highpass, time, waveform1)
plt.show()


#Filter out the 2 Hz wave (should be nearly identical to waveform2)
lowpass = filt.butter_bandpass_filter(waveform, 0.1, 0.7, fs, order= 5)
plt.plot(time, lowpass, time, waveform2)
plt.show()


