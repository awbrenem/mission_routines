#Look for deviations in the DC magnetometer indicative of signficant field-aligned currents 


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
#import end_data_loader
import numpy as np 
#import plot_spectrogram as ps
import matplotlib.pyplot as plt
#import plasma_params_get_flhr_freq
#import plasma_params_get_flhr_freq as dflh
from scipy import signal
from statsmodels.tsa.tsatools import detrend

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


#------------------
#magnetic field 
#------------------
magv = EFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)

#perp Bo
Bop = np.sqrt(mag[1]**2 + mag[2]**2)
Bot = mag[3]

#plt.plot(Bot,mag[1],Bot,mag[2])
#plt.plot(Bot,Bo)

sr = 1/(Bot[1]-Bot[0])



sz = 1   #seconds to smooth over
#sz = 0.1
sznum = int(np.floor(sz*sr))
rolling_mean = moving_average(Bop,sznum)

Bop2 = Bop[0:-(sznum-1)]
Bot2 = Bot[0:-(sznum-1)]
Bo_det = Bop2 - rolling_mean



plt.plot(Bot2,Bo_det)
#plt.xlim(400,420)
plt.ylim(-50,50)
plt.show()


