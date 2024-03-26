"""
See what effect an instrumental timing delay has on the measured polarization of a wave. 
For example, if the Bernstein waves are linearly polarized but are observed as elliptical due to a timing delay. 

Sampling frequency is 32 kHz, or about 31 microseconds b/t samples

Tests with 5000 Hz wave (200 microseconds/sample. 
A 5000 Hz wave has ~9 samples per wave period @32 kHz
1 microsecond - barely visible
10 microseconds - elliptical with < 0.5 ellipticity
30 microseconds - elliptical with about 0.7 ellipticity

Effect is more pronounced at lower frequencies. This is b/c the time shift is a smaller fraction of 
the wave phase 

"""

import numpy as np
import matplotlib.pyplot as plt


freq = 500.

#VLF channel sample rate
sr = 32000. 
ntimes = 100
t = np.asarray(range(ntimes))*(1/sr)

#timeshift
microsec = 30

tshift = t + microsec*1e-6
w = freq*2.*np.pi

#V2 waveform acts as if it's on times "tshift", not "t"
V1 = np.sin(w*t)
V2 = np.sin(w*tshift)

#As far as plotting waveforms, both V1 and V2 appear (on my end) to be on times "t" (not "tshift")
plt.plot(t,V1,t,V2)
plt.plot(t,V1,'*',t,V2,'*')
plt.show()


plt.plot(V1,V2)
plt.show()





