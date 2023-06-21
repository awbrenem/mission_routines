"""
Select desired gain/phase files (from Paulo) from the calibration testing.
From these files Paulo derives a fit (ax + b) that allows calibration from counts to volts.
These values (a, b) are found in the Endurance channel list document as the yellow boxes in far right column

See end_gainphase_test.py for comparisons to theoretical behavior
"""


#import sys 
#sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
#sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import matplotlib.pyplot as plt
import numpy as np


def end_load_gainphase(fn):

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'


    with open(path + fn) as f:
        lines = f.readlines()


    f = lines[0].split()  #freq in Hz
    p = lines[1].split()  #phase in deg
    g = lines[2].split()  #gain in dB
    f = [float(i) for i in f]
    p = [float(i) for i in p]
    g = [float(i) for i in g]

    #change to radians
    prad = [np.deg2rad(i) for i in p]

    #change gain from dB to linear scale for calculation of transfer function
    #From Steve Martin email on Nov 7, 2022: 
    #Gain=10^(0.05 * (opchan+gainoffset))
    #NOTE: This gives the correct max gain values for all the channels EXCEPT for HF.

    offset = 0.
    Hmag = [10**(0.05*i + offset) for i in g]



    #-----------------------------------------------
    #Remove bad data points that can occur at very low frequencies. This happens b/c 
    #the signal/noise ratio can become high leading to artificial values at very low freqs
    #during the gain/phase tests. 
    #-----------------------------------------------

    if fn == 'Endurance_Analog 1_V13D_10-10000-100.txt':
        Hmag[0] = Hmag[3]
        Hmag[1] = Hmag[3]
        Hmag[2] = Hmag[3]
    if fn == 'Endurance_Analog 1_V34D_10-10000-100.txt':
        Hmag[0] = Hmag[1]
    if fn == 'Endurance_Analog 1_V32D_10-10000-100.txt':
        Hmag[0] = Hmag[1]
    if fn == 'Endurance_Analog 1_V24D_10-10000-100.txt':
        Hmag[0] = Hmag[1]




    fig, axs = plt.subplots(3)
    axs[0].plot(f,g)
    axs[1].plot(f,Hmag)
    axs[2].plot(f,p)
    #axs[2].plot(f,prad)
    axs[0].set_title('gain/phase; \n fn='+ fn)
    axs[0].set_xscale('log')
    axs[1].set_xscale('log')
    axs[2].set_xscale('log')
    axs[0].set_yscale('linear')
    #axs[1].set_ylim(0,20)
    axs[0].set(ylabel='gain(dB)',xlabel='freq(kHz)')
    axs[1].set(ylabel='gain(linear)',xlabel='freq(kHz)')
    axs[2].set(ylabel='phase(deg)',xlabel='freq(kHz)')
    #axs[:].set_xlim(-40,10)
    plt.show()

    
    return prad, Hmag, f




