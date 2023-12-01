"""
Plot PSD of the skins to compare their response vs freq. 
"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import filter_wave_frequency as filt



#Calibrated
v1 = EFL('V1SD')
v1s, t, = v1.load_data_gainphase_corrected()
v2 = EFL('V2SD')
v2s, t, = v2.load_data_gainphase_corrected()
v3 = EFL('V3SD')
v3s, t, = v3.load_data_gainphase_corrected()
v4 = EFL('V4SD')
v4s, t, = v4.load_data_gainphase_corrected()

"""
#Raw
v1 = EFL('V1SD')
v1s, t, = v1.load_data()
v2 = EFL('V2SD')
v2s, t, = v2.load_data()
v3 = EFL('V3SD')
v3s, t, = v3.load_data()
v4 = EFL('V4SD')
v4s, t, = v4.load_data()
"""

fs = 1/(t[1]-t[0])

plt.plot(t,v1s)
plt.ylim(-0.3,-0.225)
plt.xlim(430,432)
plt.plot(t,v2s)
plt.plot(t,v3s)
plt.plot(t,v4s)


nfft = 1024
tr = [150,850]
psd1, psdf = correlation_analysis.psd(v1s, t, fs, tr, nft=nfft)
psd2, psdf = correlation_analysis.psd(v2s, t, fs, tr, nft=nfft)
psd3, psdf = correlation_analysis.psd(v3s, t, fs, tr, nft=nfft)
psd4, psdf = correlation_analysis.psd(v4s, t, fs, tr, nft=nfft)


fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)

axs[0].plot(psdf,np.sqrt(psd1))  #blue
axs[0].set_yscale('log')
#axs[0].set_yscale('linear')
#axs[0].set_ylim(0,1e-3)
#axs[0].set_ylim(1e-10,4e-6)
axs[0].set_xscale('log')
axs[0].set_xlim(1,1000)

axs[0].plot(psdf,np.sqrt(psd2))  #orange
axs[0].plot(psdf,np.sqrt(psd3))  #green
axs[0].plot(psdf,np.sqrt(psd4))  #red

plt.plot(psdf,psd1)

