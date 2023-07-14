#Testing



import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
import end_load_data
import end_functions as end
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
import plasma_params_get_density_from_flhr_freq as dflh
import plasma_params_get_flhr_freq as dflh2
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
from scipy import signal

import numpy as np 
import matplotlib.pyplot as plt
import plasmapy
from astropy import units as u  



"""Load E-field VLF data"""
evlf = end_load_data.efield_vlf()

print('here')
#fig,axs = plt.subplots(2)


good = list(map(tuple,np.where((evlf['tvlf'] > 300) & (evlf['tvlf'] < 305))))

goo1 = evlf['dvlf12_mvm'][good]
t1 = evlf['tvlf'][good]
goo2 = evlf['dvlf12_mvm_gpcal'][good]

plt.plot(t1[0],goo1[0],t1[0],goo2[0])



