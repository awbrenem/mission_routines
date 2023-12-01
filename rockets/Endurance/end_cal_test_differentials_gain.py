"""
Plot PSD of the differential channels to compare their response vs freq. 
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
from scipy.io import readsav 




#raw voltages converted into mV/m
folder = "efield_DC"
fn = "47001_TM1_LFDSP_S5DCE_DCES5_calibrated.sav"
path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'

vals = readsav(path + folder + '/' + fn)

bld = 2.26
bl = 3.2

t = vals.times
v12r = 1000*vals.dv12_volts / bl
v34r = 1000*vals.dv34_volts / bl
v13r = 1000*vals.dv13_volts / bld
v32r = 1000*vals.dv32_volts / bld
v24r = 1000*vals.dv24_volts / bld
v41r = 1000*vals.dv41_volts / bld
fs = 10000



#Gainphase Calibrated mV/m
vt = EFL('V12D')
v12c, t, = vt.load_data_gainphase_corrected()
vt = EFL('V34D')
v34c, t, = vt.load_data_gainphase_corrected()
vt = EFL('V13D')
v13c, t, = vt.load_data_gainphase_corrected()
vt = EFL('V32D')
v32c, t, = vt.load_data_gainphase_corrected()
vt = EFL('V24D')
v24c, t, = vt.load_data_gainphase_corrected()
vt = EFL('V41D')
v41c, t, = vt.load_data_gainphase_corrected()


#Raw mV/m
vt = EFL('V12D')
v12, t, = vt.load_data()
vt = EFL('V34D')
v34, t, = vt.load_data()
vt = EFL('V13D')
v13, t, = vt.load_data()
vt = EFL('V32D')
v32, t, = vt.load_data()
vt = EFL('V24D')
v24, t, = vt.load_data()
vt = EFL('V41D')
v41, t, = vt.load_data()




#---------------
#comparison of the different types (gp calibrated, non-calibrated, raw, from skins)
#for a single channel

nfft = 1024
tr = [150,850]
psd12r, psdfr = correlation_analysis.psd(v12r, t, fs, tr, nft=nfft)
psd12c, psdfc = correlation_analysis.psd(v12c, t, fs, tr, nft=nfft)
psd12, psdf = correlation_analysis.psd(v12, t, fs, tr, nft=nfft)


fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.4,hspace=0.4)

axs[0].set_yscale('linear')
axs[0].set_xscale('log')
axs[0].set_xlim(1,10000)
axs[0].plot(psdfr,psd12r,'*')  #blue
axs[0].plot(psdfc,psd12c,'*')  #orange
axs[0].plot(psdf,psd12,'*')  #green










nfft = 1024
#tr = [150,850]
tr = [610,620]
psd12c, psdf = correlation_analysis.psd(v12c, t, fs, tr, nft=nfft)
psd34c, psdf = correlation_analysis.psd(v34c, t, fs, tr, nft=nfft)
psd13c, psdf = correlation_analysis.psd(v13c, t, fs, tr, nft=nfft)
psd32c, psdf = correlation_analysis.psd(v32c, t, fs, tr, nft=nfft)
psd24c, psdf = correlation_analysis.psd(v24c, t, fs, tr, nft=nfft)
psd41c, psdf = correlation_analysis.psd(v41c, t, fs, tr, nft=nfft)


fig, axs = plt.subplots(2,figsize=(8,6))
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.4,hspace=0.4)

axs[0].plot(psdf,psd12c)  #blue
axs[0].set_yscale('linear')
axs[0].set_xscale('log')
axs[0].set_xlim(1,10000)

axs[0].plot(psdf,psd34c)  #orange
axs[0].plot(psdf,psd13c)  #green
axs[0].plot(psdf,psd32c)  #red
axs[0].plot(psdf,psd24c)  #purple
axs[0].plot(psdf,psd41c)  #brown






plt.plot(t,v12)
plt.ylim(-0.3,-0.225)
plt.xlim(430,432)
plt.plot(t,v34)
plt.plot(t,v13)
plt.plot(t,v32)
