"""
Test the Endurance transfer function on skins and DC data

"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
from end_fields_loader import Endurance_Fields_Loader as EFL



#---------------------------------------------
#Load corrected data
#---------------------------------------------

#For comparison
v1c = EFL('V12D')
wf1c, tdat1c = v1c.load_data_gainphase_corrected()
fs1c = v1c.chnspecs['fs']

#skin 1
v2c = EFL('V1SD')
wf2c, tdat2c = v2c.load_data_gainphase_corrected()
fs2c = v2c.chnspecs['fs']

#skin 2
v3c = EFL('V2SD')
wf3c, tdat3c = v3c.load_data_gainphase_corrected()
fs3c = v3c.chnspecs['fs']



#-------------------------------
#Construct E-fields from skin channels
#-------------------------------




#get rid of slight DC offset before I use skins to create Efield
hp = 2  #Hz
wf2c = filt.butter_highpass_filter(wf2c, hp, fs2c, order=8)
wf3c = filt.butter_highpass_filter(wf3c, hp, fs2c, order=8)




wfVskinsc = (wf2c - wf3c)


#Redefine the second channel as the Efield from skins
#This allows the rest of the program to function without change    
wf2c = wfVskinsc




#boom lengths (for reference)
#bl1234 = 3.212
#bl13413224 = 2.271227
wfEskinsc = 1000*(wfVskinsc)/3.212
#wfEskinsc = 1000*(wfVskinsc)/2.271227




#------------------------------
#Plot comparison
#------------------------------


tr = [400,405]
goo1 = np.where((tdat1c > tr[0]) & (tdat1c < tr[1]))
goo2 = np.where((tdat2c > tr[0]) & (tdat2c < tr[1]))

#fig,axs = plt.subplots(3)
#axs[0].plot(tdat1c[goo1], wf1c[goo1])
#axs[1].plot(tdat2c[goo2], wfVskinsc[goo2])
#axs[2].plot(tdat2c[goo2], wfEskinsc[goo2])
#for i in range(3): axs[i].set_xlim(tr[0],tr[1])


#get rid of DC offset and plot
wf1chp = filt.butter_highpass_filter(wf1c, hp, fs1c, order=8)

plt.plot(tdat2c[goo2], wfEskinsc[goo2], tdat1c[goo1], wf1chp[goo1])


print("here")
