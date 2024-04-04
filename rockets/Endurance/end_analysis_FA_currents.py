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




#------------------
#magnetic field 
#------------------
magv = EFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]

plt.plot(Bot,mag[1],Bot,mag[2])
plt.plot(Bot,Bo)

Bo_det = signal.detrend(Bo,type='constant')
plt.plot(Bot,Bo_det)
plt.show()


