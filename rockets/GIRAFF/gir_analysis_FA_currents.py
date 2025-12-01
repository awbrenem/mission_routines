"""
Estimate magnitude of FA currents on GIRAFF using the curlometer technique

References:
Kintner+05 (SCIFER-2 proposal)
Kintner+96 (SCIFER)

muo*Jz = dBy/dx - dBx/dy = (1/Vx)*dBy/dt - (1/Vy)*dBx/dt


NOTE: Attempting to do a rough Jz calculating using the spinning data 
isn't working b/c the signal is completely dominated by the spin modulation
"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import gir_load_data as gld
#from scipy import signal
import numpy as np 
#import correlation_analysis as ca
import plot_spectrogram as ps
import matplotlib.pyplot as plt
#import plasma_params_get_flhr_freq as dflh
#import pyIGRF
#import pickle
#import scipy.io as sio
#from scipy.io import readsav
#import fft_spectrum_piecewise as fftspec
from scipy.interpolate import interp1d
import numpy as np
#from scipy.signal import find_peaks, peak_widths, peak_prominences
#from scipy import interpolate
#from scipy.signal import find_peaks
#from scipy.signal import peak_widths



#******************************
#NOTE: Need to despin the B-field data first so I can use it alongside the Velocity data
#******************************




#ephem, ephem2 = gld.load_ephemeris()



def remove_spikes(data, threshold=5, window=5):
    """
    Simple spike removal: replaces points that deviate from a moving median by more than threshold*std
    """
    from scipy.ndimage import median_filter
    med = median_filter(data, size=window)
    diff = np.abs(data - med)
    std = np.nanstd(data)
    spikes = np.where(diff > threshold * std)[0]
    data_despiked = data.copy()
    data_despiked[spikes] = med[spikes]
    return data_despiked


pld = '381'
chn = 'mag'


#Load B-fields data
Bchn = GFL(pld,chn)
bx,by,bz,tdat = Bchn.load_data()



bx_despiked = remove_spikes(bx, threshold=0.5, window=31)
by_despiked = remove_spikes(by, threshold=0.5, window=31)
bz_despiked = remove_spikes(bz, threshold=0.5, window=31)



# Plot original and despiked data
plt.figure()
plt.plot(tdat, bx, label='Original bx')
plt.plot(tdat, bx_despiked, label='Despiked bx')
plt.xlim(400, 500)
plt.legend()
plt.title('Spike removal on bx')
plt.show()



dt = np.gradient(tdat)
dbx_dt = np.gradient(bx_despiked) / dt
dby_dt = np.gradient(by_despiked) / dt



Jz_rough = dby_dt - dbx_dt


plt.figure()
plt.plot(tdat, Jz_rough, label='Rough Jz estimate')
plt.ylim(-2e5,2e5)
plt.show()
#muo*Jz = dBy/dx - dBx/dy = (1/Vx)*dBy/dt - (1/Vy)*dBx/dt






#bx_despiked = spike_removal(bx, 
#                            width_threshold = 3,         
#                            prominence_threshold = 5,   #20
#                            moving_average_window=40, 
#                            width_param_rel=1.2, 
#                            interp_type='linear')



#plt.plot(tdat,bx)
plt.plot(tdat,bx_despiked)
plt.plot(tdat,bx)
plt.xlim(400,500)





dbx = np.gradient(bx)
plt.plot(tdat,dbx)
plt.xlim(405,410)

bad = np.where(np.abs(dbx) > 500)[0]
bx2 = bx
bx2[bad] = np.nan

plt.plot(tdat,bx2)
plt.plot(tdat,bx)
plt.xlim(500,550)
plt.xlim(100,550)


dbx2 = np.gradient(bx2)
plt.plot(tdat,dbx2)
plt.xlim(100,550)
plt.ylim(-1000,1000)




