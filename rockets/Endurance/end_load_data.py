#-------------------------------------------------------------------------
#  _______ ______  _____   _     _ ______         ______   ______ _______
# (_______|  ___ \(____ \ | |   | (_____ \   /\  |  ___ \ / _____(_______)
#  _____  | |   | |_   \ \| |   | |_____) ) /  \ | |   | | /      _____
# |  ___) | |   | | |   | | |   | (_____ ( / /\ \| |   | | |     |  ___)
# | |_____| |   | | |__/ /| |___| |     | | |__| | |   | | \_____| |_____
# |_______|_|   |_|_____/  \______|     |_|______|_|   |_|\______|_______)
#
#--------------------------------------------------------------------------
# 9/16/2022 - A.W. Breneman
#   - Loads Endurance data
#   
#   T-0: May 11, 2022  01:31:00.0 U.T.
#   Data array original sample rate 30kHz
#   Column 1, time (sec) since T-0, format:f13.7
#   Column 2, VLF12 (mV/m), format:f10.3
#--------------------------------------------------------------------------


import pandas as pd 
import numpy as np
from scipy import signal
from scipy.io import readsav 
import matplotlib.pyplot as plt
import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps

#from scipy.fft import fftshift

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'
folder = 'ephemeris'
fn = 'Endurance_GPS_velocity_position_altitude.csv'

header = ['time','xvel','yvel','zvel','lat','long','alt']
ephem =  pd.read_csv(path + folder + '/' + fn, skiprows=1, names=header)


#vmag = np.sqrt(ephem.xvel**2 + ephem.yvel**2 + ephem.zvel**2)

#fig, axs = plt.subplots(4)
#axs[0].plot(ephem.time,vmag)
#axs[1].plot(ephem.time,ephem.alt)
#axs[2].plot(ephem.time,ephem.lat)
#axs[3].plot(ephem.time,ephem.long)


def efield_dc():

    folder = "efield_DC"
    fn = "47001_TM1_LFDSP_S5DCE_DCES5_calibrated.sav"

    edc = readsav(path + folder + '/' + fn)
    edckeys = list(edc.keys())
    #dict_keys(['times', 'bl1234', 'bl13413224', 'dv12_mvm', 'dv12_volts', 'dv12_raw', 'dv12_rawnormal', 'a12', 'b12', 'dv34_mvm', 'dv34_volts', 'dv34_raw', 'dv34_rawnormal', 'a34', 'b34', 'dv13_mvm', 'dv13_volts', 'dv13_raw', 'dv13_rawnormal', 'a13', 'b13', 'dv32_mvm', 'dv32_volts', 'dv32_raw', 'dv32_rawnormal', 'a32', 'b32', 'dv24_mvm', 'dv24_volts', 'dv24_raw', 'dv24_rawnormal', 'a24', 'b24', 'dv41_mvm', 'dv41_volts', 'dv41_raw', 'dv41_rawnormal', 'a41', 'b41', 'author', 'timenote', 'calnote', 'dataunits', 'flight', 'format', 'filein', 'link', 'samplerate', 't0', 'timetagmethod', 'timeunits'])


    #Get rid of all times t<0
    good = np.squeeze(np.where(edc.times >= 0.))
    tstlen = len(edc['times'])

    for i in edckeys:
        if np.shape(edc[i]) != ():
            if len(edc[i]) == tstlen:
                edc[i] = edc[i][good]

    return edc



def efield_vlf():


    folder = 'efield_VLF'

    """
    fn = '47001_TM1_LFDSP_VLF12_Glyn_30kHz_v1.txt'
    header = ['tsec', 'amp']
    vlf12 = pd.read_csv(path + folder + '/' + fn, skiprows=17, names=header, delim_whitespace=True)

    
    #fs = 1/(np.mean(vlf12.tsec - vlf12.tsec.shift()))
    #freq12, tspec12, power12 = signal.spectrogram(vlf12.amp, fs, nperseg=512, return_onesided=1,window='hann')
    #ps.plot_spectrogram(tspec12,freq12,power12,vr=[0.4,0.65], xr=[0,900],yr=[1,10000],pl=1)
    #ps.plot_spectrogram(tspec12,freq12,power12,vr=[0.4,0.65], xr=[0,900],yr=[0,10000],pl=1)

    """


    fn = '47001_TM1_LFDSP_S5_VLF_mvm.sav'

    vlf12 = readsav(path + folder + '/' + fn)
    #dict_keys(['tvlf', 'dvlf12_mvm', 'dvlf34_mvm', 'dvlf24_mvm', 'dvlf32_mvm', 'author', 
    # 'calnote_vlf1234', 'calnote_vlf2432', 'dataunits', 'flight', 'format', 'filein', 'link', 
    # 'samplerate', 't0', 'timetagmethod', 'timeunits'])


    #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
    good = np.squeeze(np.where(vlf12.tvlf >= 0.))

    vlf12['tvlf'] = vlf12['tvlf'][good]
    vlf12['dvlf12_mvm'] = vlf12['dvlf12_mvm'][good]
    vlf12['dvlf34_mvm'] = vlf12['dvlf34_mvm'][good]
    vlf12['dvlf24_mvm'] = vlf12['dvlf24_mvm'][good]
    vlf12['dvlf32_mvm'] = vlf12['dvlf32_mvm'][good]

    #fs = evlf.samplerate
    #freq12, tspec12, power12 = signal.spectrogram(evlf.dvlf12_mvm, fs, nperseg=512, return_onesided=1,window='hann')
    #ps.plot_spectrogram(tspec12,freq12,power12,vr=[0.4,0.65], xr=[0,900],yr=[0,10000],pl=1)


    return vlf12


def efield_hf():

    """
    Some notes on reading the arrays AFFTPOW12, AFFTPOW34:
    
    AFFTPOW12 corresponds to HF12, and AFFTPOW34 to HF34
    There are 945 spectra for each channel, corresponding time stamps in ATIMESFFT
    Frequencies are arranged low to high, so AFFTPOW12[*,0] corresponds to DC, and AFFTPOW12[*,4200] corresponds to 4 Mhz.
    
    Oh, and the way I define power (as per Rob) is as follows:
    
    AFFTPOW = 10log10(Win*FFT^2*factor)
    
    Where Win=window function (typically Hanning) and factor is a constant.
    """ 


    folder = 'efield_HF'
    fn = '47001_TM2_0-LOS_S1_HFsnippets_HF1234_method2_mblkskips_4_FFT8400.sav'
    ehf = readsav(path + folder + '/' + fn)

    ehf.keys()
    #dict_keys(['afftpow12', 'afftpow34', 'afreq', 'atimesfft', 'fftsize', 'overlap', 'weight', 'samplerate', 'nfreq', 'nfftlines', 'in_file'])




    #hf12 = ehf.afftpow12
    #hf34 = ehf.afftpow34
    #hffreqs = ehf.afreq
    #hftimes = ehf.atimesfft

    #ps.plot_spectrogram(hftimes,hffreqs,hf12,pl=0,vr=[0.4,0.48],yr=[1.4e6,4.5e6],xr=[0,900])
    #ps.plot_spectrogram(hftimes,hffreqs,hf12,pl=0,yscale='log',vr=[0.35,0.48],yr=[2e3,4e6],xr=[0,900])
  

    return ehf
    #Slice at t=810


    """
    #fn = '47001_TM2_0-LOS_S1_HFsnippets_HF1234_method2.sav'
    ehf = readsav(path + folder + '/' + fn)

    v1 = ehf['dhf12']
    v2 = ehf['dhf34']


    hf = {'dhf12':[v1], 'dhf34':[v2]}
    hf = pd.DataFrame.from_dict(hf)

    plt.plot(v1[:,1])

    fs = 8e6
    freq12, tspec12, power12 = signal.spectrogram(hf.dhf12, fs, nperseg=512, return_onesided=1,window='hann')

    """


"""
Return DC magnetometer values. 

*******************************************
Contents of this file:
ENDURANCE 47.001
LFDSP TM1 DC Magnetometer Channels
T-0: May 11, 2022  01:31:00.0 U.T.
Data arrays original sample rate 2kHz
*** Mag X, Y, Z labels may NOT correspond
    to sensor coordinates ***
Column 1, time (sec) since T-0, format:f13.7
Column 2, Mag X (nT), format:f10.3
Column 3, Mag Y (nT), format:f10.3
Column 4, Mag Z (nT), format:f10.3
*******************************************

"""
def mag_dc():

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'
    folder = 'mag'
    fn = '47001_TM1_LFDSP_S5Mag_MagXYZS5_nominalcal_Glyn_2kHz_v1.txt'

    header = ['tsec','Bx','By','Bz']
    x = pd.read_csv(path+folder+'/'+fn, skiprows=21, names=header, delim_whitespace=True)

    bm = np.sqrt(x.Bx**2 + x.By**2 + x.Bz**2)
    bmag = pd.DataFrame(bm,columns=['Bmag'])
    x['Bmag'] = bmag


    #plt.plot(x.tsec, x.Bmag)

    return x



#allow to be run as a script
if __name__ == '__main__':
    print("<> Running end_load_data.py as a script! <>")
    #x = load_particle("eea")
    #print(x)
    x = efield_vlf()
    b = mag_dc()
    x = efield_hf()
    x.keys()
    print("here_fin")





