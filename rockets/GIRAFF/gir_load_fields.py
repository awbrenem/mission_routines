# 2025-05-22 - A.W. Breneman
#   - Loads GIRAFF E-fields data
#   

#**********FIX THE BELOW DATA**********
#   T-0: May 11, 2022  01:31:00.0 U.T.
#   Data array original sample rate 30kHz
#   Column 1, time (sec) since T-0, format:f13.7
#   Column 2, VLF12 (mV/m), format:f10.3


"""
Note on channel polarity:

Polarity is accounted for by Steve from the conversion equation ("a" coefficient) in the channel specs document. 

"""
#--------------------------------------------------------------------------


class GIRAFF_Fields_Loader:

    """
    Loadable channels:
    DC: 'V12D','V34D'
    skins: 'V1SD','V2SD','V3SD','V4SD'
    VLF: 'VLF12D','VLF34D','VLF13D','VLF32D','VLF24D','VLF41D'
    HF: 'HF12','HF34'
    analog: 'V12A','V34A','VLF12A','V1SA','V2SA','V3SA','V4SA'
    mag: 'magX','magY','magZ'
    Langmuir Probes (analog, digital): 'LPA', 'LPD'
    
    Example usage:

    from gir_fields_loader import GIRAFF_Fields_Loader

    v12 = GIRAFF_Fields_Loader('V12D')
    v12.plot_gainphase('381')
    v12dat = v12.load_data()

    """


    def __init__(self, pld, chn) -> None:

        self.pld = pld  #Payload specification (380 or 381)
        self.chn = chn  #Channel specification 
        self._gir_channelspecs()


        if (self.chn == 'V12D') or (self.chn == 'V34D'): 
            self.type = 'DC'
            filterpole = 4
        if (self.chn == 'V1SD') or (self.chn == 'V2SD') or (self.chn == 'V3SD') or (self.chn == 'V4SD'): 
            self.type = 'skins'
            filterpole = 1
        if (self.chn == 'VLF12D') or (self.chn == 'VLF34D') or (self.chn == 'VLF13D') or (self.chn == 'VLF32D') or (self.chn == 'VLF24D') or (self.chn == 'VLF41D'): 
            self.type = 'VLF'
            filterpole = 4
        if (self.chn == 'HF12') or (self.chn == 'HF34'):
            self.type = 'HF'
            filterpole = 2
        if (self.chn == 'V12A') or (self.chn == 'V34A') or (self.chn == 'VLF12A') or (self.chn == 'V1SA') or (self.chn == 'V2SA') or (self.chn == 'V3SA') or (self.chn == 'V4SA'):
            self.type = 'analog'
            filterpole = "nan"
        if (self.chn == 'magX') or (self.chn == 'magY') or (self.chn == 'magZ'):
            self.type = 'mag'
            filterpole = "nan"
        if (self.chn == 'LPA'):
            self.type = 'LPA'
            filterpole = "nan"
        if (self.chn == 'LPD'):
            self.type = 'LPD'
            filterpole = "nan"


        #**********NEED TO FINISH THIS ROUTINE*************
        #if self.type != 'mag':
        #    self._gir_load_gainphase()



        self.chnspecs["folder"] = "efield_" + self.type
        #low pass filter pole
        self.chnspecs["lpf_pole"] = filterpole

        return 
    

    def __str__(self):
        return "GIRAFF " + self.chn + " object"




    #---------------------------------------------------------------------------------
    #Load gain/phase calibrated files. These are produced by gir_cal_transfer_function.py.

    #NOTE: Unfortunately the time values are not sampled evenly. Over the course of the mission this causes any spectrogram
    #that is produced from the waveform data to be significantly off (e.g. 2 sec by end of mission). 
    #This is b/c the spectrograms use a single "fs" (sample freq) value. To get around this, interpolate the times to a regular grid
    #NOTE 2: For future rocket missions this step should be done prior to applying gain/phase calibration.

    #---------------------------------------------------------------------------------

    def load_data_gainphase_corrected(self):

        import pickle 
        import numpy as np
        path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/'


        if self.type == 'DC':
            if self.type == '380':
                goo = pickle.load(open(path + 'efield_DC/' + 'Giraff_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
            if self.type == '381':
                goo = pickle.load(open(path + 'efield_DC/' + 'Giraff_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
        if self.type == 'skins':
            if self.type == '380':
                goo = pickle.load(open(path + 'efield_skins/' + 'Giraff_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
            if self.type == '381':
                goo = pickle.load(open(path + 'efield_skins/' + 'Giraff_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))

        if self.type == 'VLF':
            if self.type == '380':
                goo = pickle.load(open(path + 'efield_VLF/' + 'Giraff_Analog 1_' + self.chn + '_6-30000-100_gainphase_corrected.pkl', 'rb'))
            if self.type == '381':
                goo = pickle.load(open(path + 'efield_skins/' + 'Giraff_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
        
        return goo['wf'], goo['tvals']



    #---------------------------------------------------------------------------------
    #Load the data sent to me by Steve Martin. These are gain corrected but not 
    #phase corrected (I use the unity gain to "gain/phase" correct these). 
    #plt = '380' or '381' for the two payloads
    #---------------------------------------------------------------------------------

    def load_data(self):
        import numpy as np
        from scipy.io import readsav 

        path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/'


        if self.type == 'analog':
            print('no analog data from Steve Martin available to load...RETURNING')
            return -1

        if self.type == 'skins':

            folder = "efield_skins"

            if self.pld == '380':
                fn = "36380_TM1_LFDSP_Skins_S9_volts.sav"
            if self.pld == '381':
                fn = "36381_TM1_LFDSP_Skins_S9_volts.sav"
            
            vals = readsav(path + folder + '/' + fn)

            pol = 1
            if self.chnspecs["polarity"] == 'Neg':
                pol = -1

            if self.chn == 'V1SD': 
                 wf = vals.dv1sd_volts * pol
            if self.chn == 'V2SD': 
                 wf = vals.dv2sd_volts * pol
            if self.chn == 'V3SD': 
                 wf = vals.dv3sd_volts * pol
            if self.chn == 'V4SD': 
                 wf = vals.dv4sd_volts * pol

            t = vals.tskins


            #Get rid of all times t<0
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good]


        if self.type  == "DC":

            folder = "efield_DC"

            if self.pld == '380':
                fn = "36380_TM2_V12D_V34D_mvm.sav"
            if self.pld == '381': 
                fn = "36381_TM2_V12D_V34D_mvm.sav"
            
            vals = readsav(path + folder + '/' + fn)

            t = vals.tv12v34
            pol = 1 
            if self.chnspecs['polarity'] == 'Neg':
                pol = -1

            if self.chn == 'V12D': wf = vals.dv12d_mvm * pol
            if self.chn == 'V34D': wf = vals.dv34d_mvm * pol

            #Get rid of all times t<0
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good]


        if self.type == 'VLF':

            folder = 'efield_VLF'
            
            if self.pld == '380':
                if self.chn == 'VLF12D': 
                    fn = '36380_GIRAFF_TM2_LFDSP_VLF12D_mvm.sav'
                if self.chn == 'VLF34D':
                    fn = '36380_GIRAFF_TM2_LFDSP_VLF34D_mvm.sav'
            if self.pld == '381':
                
                #****WHICH FILE SHOULD I BE LOADING???
                #if self.chn == 'VLF12D':
                #    fn = '36381_GIRAFF_TM2_LFDSP_IT_VLF12D_mvm_FFT131072.sav'
                if self.chn == 'VLF12D':
                    fn = '36381_GIRAFF_TM2_LFDSP_IT_VLF12D_mvm.sav'
                if self.chn == 'VLF34D':
                    fn = '36381_GIRAFF_TM2_LFDSP_IT_VLF34D_mvm.sav'


            vals = readsav(path + folder + '/' + fn)

            t = vals.tvlf12vlf34
            pol = 1
            if self.chnspecs["polarity"] == 'Neg':
                pol = -1

            if self.chn == 'VLF12D': wf = vals.dvlf12d_mvm * pol
            if self.chn == 'VLF34D': wf = vals.dvlf34d_mvm * pol
            if self.chn == 'VLF13D': wf = vals.dvlf13d_mvm * pol
            if self.chn == 'VLF32D': wf = vals.dvlf32d_mvm * pol
            if self.chn == 'VLF24D': wf = vals.dvlf24d_mvm * pol
            if self.chn == 'VLF41D': wf = vals.dvlf41d_mvm * pol


            #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good]


        if self.type == 'HF':

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
            
            #****WHICH FILE SHOULD I BE LOADING???
            if self.pld == '380':
                #fn = '36380_TM1_HF12HF34_S9_counts_FFT2120.sav'
                fn = '36380_TM1_HF12HF34_S9_counts.sav'
            if self.pld == '381':
                #fn = '36381_TM1_HF12HF34_S9_counts_FFT2120.sav'
                fn = '36381_TM1_HF12HF34_S9_counts.sav'


            vals = readsav(path + folder + '/' + fn)

            t = vals.atimesfft
            f = vals.afreq
            if self.chn == 'HF12': wf = vals.afftpow12
            if self.chn == 'HF34': wf = vals.afftpow34

            #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good], f[good]

        

        #Return mag data (nT)
        if self.type == 'mag':   
            from scipy.io import readsav 

            folder = 'mag'

            if self.pld == '380':
                fn = '36380_TM1_MagXYZ_S9_nT.sav'
            if self.pld == '381':
                fn = '36381_TM1_MagXYZ_S9_nT.sav'

            vals = readsav(path + folder + '/' + fn)

            #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
            good = np.squeeze(np.where(vals.times >= 0.))
            return vals.dmagxnt[good], vals.dmagynt[good], vals.dmagznt[good], vals.times[good]





    #-------------------------------------------------------------------------------------
    #load the gain/phase files for desired channel

    """
    NOTES: Select desired gain/phase files (from Paulo) from the calibration testing.
    From these files Paulo derives a fit (ax + b) that allows calibration from counts to volts.
    These values (a, b) are found in the Giraff channel list document as the yellow boxes in far right column

    See end_gainphase_test.py for comparisons to theoretical behavior
    """
    #-------------------------------------------------------------------------------------

    def _gir_load_gainphase(self):

        import numpy as np
        from math import remainder

        path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/gain_phase_files/'

        #Load gain/phase file
        if self.chn == "V1SD": fn = "Giraff_Analog 1_V1SD_10-10000-100.txt"
        elif self.chn == "V2SD": fn = "Giraff_Analog 1_V2SD_10-10000-100.txt"
        elif self.chn == "V3SD": fn = "Giraff_Analog 1_V3SD_10-10000-100.txt"
        elif self.chn == "V4SD": fn = "Giraff_Analog 1_V4SD_10-10000-100.txt"
        elif self.chn == "V1SA": fn = "Giraff_Analog 1_V1SA_10-10000-100.txt"
        elif self.chn == "V2SA": fn = "Giraff_Analog 1_V2SA_10-10000-100.txt"
        elif self.chn == "V3SA": fn = "Giraff_Analog 1_V3SA_10-10000-100.txt"
        elif self.chn == "V4SA": fn = "Giraff_Analog 1_V4SA_10-10000-100.txt"
        elif self.chn == "V12D": fn = "Giraff_Analog 1_V12D_10-10000-100.txt"
        elif self.chn == "V34D": fn = "Giraff_Analog 1_V34D_10-10000-100.txt"
        elif self.chn == "V13D": fn = "Giraff_Analog 1_V13D_10-10000-100.txt"
        elif self.chn == "V32D": fn = "Giraff_Analog 1_V32D_10-10000-100.txt"
        elif self.chn == "V24D": fn = "Giraff_Analog 1_V24D_10-10000-100.txt"
        elif self.chn == "V41D": fn = "Giraff_Analog 1_V41D_10-10000-100.txt"
        elif self.chn == "VLF12D": fn = "Giraff_Analog 1_VLF12D_6-30000-100.txt"
        elif self.chn == "VLF34D": fn = "Giraff_Analog 1_VLF34D_6-30000-100.txt"
        elif self.chn == "VLF13D": fn = "Giraff_Analog 1_VLF13D_6-30000-100.txt"
        elif self.chn == "VLF32D": fn = "Giraff_Analog 1_VLF32D_6-30000-100.txt"
        elif self.chn == "VLF24D": fn = "Giraff_Analog 1_VLF24D_6-30000-100.txt"
        elif self.chn == "VLF41D": fn = "Giraff_Analog 1_VLF41D_6-30000-100.txt"
        elif self.chn == "VLF12A": fn = "Giraff_Analog 1_VLF12A_6-100000-100.txt"
        elif self.chn == "V12A": fn = "Giraff_Analog 1_V12A_10-10000-100.txt"
        elif self.chn == "V34A": fn = "Giraff_Analog 1_V34A_10-10000-100.txt"
        elif self.chn == "HF12": fn = "Giraff_Analog 1_HF12_1000-20000000-100.txt"
        elif self.chn == "HF34": fn = "Giraff_Analog 1_HF34_1000-20000000-100.txt"
#        elif self.chn == "mag": fn = "placeholder.txt"

        self.chnspecs["gainphase_file"] = fn
        self.chnspecs["gainphase_path"] = path


        with open(path + fn) as f:
            lines = f.readlines()


        f = lines[0].split()  #freq in Hz
        p = lines[1].split()  #phase in deg
        g = lines[2].split()  #gain in dB
        f = [float(i) for i in f]
        p = [float(i) for i in p]
        g = [float(i) for i in g]

        #--change to radians
        prad = [np.deg2rad(i) for i in p]


        #--"Fix" the phase data so that it can only vary b/t -pi and pi. Some of the curves 
        #--extend a bit beyond this for some reason. To fix, unwrap and then rewrap angles.
        prad = np.unwrap(prad)
        prad = np.asarray([remainder(prad[i], 2*np.pi) for i in range(len(prad))])


        #--change gain from dB to linear scale for calculation of transfer function
        #--From Steve Martin email on Nov 7, 2022: 
        #--Gain=10^(0.05 * (opchan+gainoffset))
        #--NOTE: This gives the correct max gain values for all the channels EXCEPT for HF.
        offset = 0.
        Hmag = [10**(0.05*i + offset) for i in g]

 

        #--Remove bad data points that can occur at very low frequencies. This happens b/c 
        #--the signal/noise ratio can become high leading to artificial values at very low freqs
        #--during the gain/phase tests. 
        if self.chn == 'V13D':
            Hmag[0] = Hmag[3]
            Hmag[1] = Hmag[3]
            Hmag[2] = Hmag[3]
        if self.chn == 'V34D':
            Hmag[0] = Hmag[1]
        if self.chn == 'V32D':
            Hmag[0] = Hmag[1]
        if self.chn == 'V24D':
            Hmag[0] = Hmag[1]

        self.phase = prad 
        self.gain = Hmag 
        self.gaindB = g
        self.freq_gainphase = f     



    
    #-------------------------------------------------
    #Plot the gain phase curves for requested channel
    #-------------------------------------------------

    def plot_gainphase(self):

        import matplotlib.pyplot as plt
        import numpy as np


        index3dB = self._determine_3dB()


        #slope of idealized filter
        decades = 1
        pole = self.chnspecs["lpf_pole"]

        slope = -20*pole   #20 dB/decade of normalized freq


        x1 = self.freq_gainphase[index3dB]
        x2 = x1*10*decades

        y1 = self.gaindB[index3dB]
        y2 = slope + y1

        xv = [x1,x2]
        yv = [y1,y2]


        fig, axs = plt.subplots(3)
        axs[0].plot(self.freq_gainphase,self.gaindB)
        axs[1].plot(self.freq_gainphase,self.gain)
        axs[2].plot(self.freq_gainphase,self.phase)
        for i in range(3):
            axs[i].axvline(self.chnspecs["fs"]/2,color='r',linestyle='--')

        #3dB point (actual)
        axs[0].plot(self.freq_gainphase[index3dB], self.gaindB[index3dB], 'X')
        #3dB point (expected)
        axs[0].plot(self.chnspecs["lpf"], self.chnspecs["filter_3dB_dBvalue_measured"], 'o',markersize=4)


        #idealized filter falloff
        axs[0].plot(xv,yv,linestyle='--')

        axs[0].set_title('gain/phase; \n fn='+ self.chnspecs["gainphase_file"])
        for i in range(3):
            axs[i].set_xscale('log')
            axs[i].set_xlim(1, 50000)
        
        axs[0].set_yscale('linear')
        axs[0].set(ylabel='gain(dB)',xlabel='freq(kHz)')
        axs[1].set(ylabel='gain(linear)',xlabel='freq(kHz)')
        axs[2].set(ylabel='phase(deg)',xlabel='freq(kHz)')
        plt.show()



    #--------------------------------------------------------------
    #Returns actual 3dB value below peak response and its frequency 
    #--------------------------------------------------------------

    def _determine_3dB(self):

        import numpy as np 

        #determine 3dB point. 
        threedB = np.nanmax(self.gaindB) - 3
        c1 = self.gaindB <= threedB
        c2 = np.asarray(self.freq_gainphase) > self.chnspecs['fs']/6
        c3 = (c1 & c2)
        index3dB = [i for i in range(len(c3)) if c3[i] == True][0]

        #Expected 3dB point 
        #threedBfreq = self.chnspecs["lpf"]
        #goo = [i for i in range(len(self.freq_gainphase)) if self.freq_gainphase[i] >= threedBfreq]
        #index3dB_exp = goo[0]

        self.chnspecs["filter_3dB_freq_measured"] = self.freq_gainphase[index3dB]
        self.chnspecs["filter_3dB_dBvalue_measured"] = threedB


        return index3dB


    #-------------------------------------------------------------
    #Return channel specifications from the MEB for desired channel
    #-------------------------------------------------------------

    def _gir_channelspecs(self):

        import numpy as np
        import sys  
        sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
        import gir_cal_channel_data as gcd

        """
        #NOTE: For reference, use this code to "quickly" read in the MEB file in order to populate values
        
        import tabula
        path = "/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/Endurance Channel Specifications DC AC Cal 2-11-2022.pdf"
        df = tabula.read_pdf(path)
        df.rename(columns={})
        d1 = df[0]
        d2 = df[1]
        df1.keys()
            Index(['Unnamed: 0', 'Unnamed: 1', 'Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4',
            'Unnamed: 5', 'Unnamed: 6', 'Unnamed: 7', 'Unnamed: 8', 'Unnamed: 9',
            'Unnamed: 10', 'Desired\rGaindB', 'Measured\rGaindB', 'Unnamed: 11',
            'Unnamed: 12', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15',
            'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18'],
            dtype='object')
        nrows = 25 
        for i in range(nrows):
            print(d1['Measured\rGaindB'][i], end=",")

        """

        #Get basic data for every channel. This is identical b/t the two payloads, and so I'll just load "380" here
        chn_dat380, goo = gcd.gir_cal_channel_data()


        #extract channel names
        chn_names = [""]*len(chn_dat380)        
        for i in range((len(chn_dat380))):
            chn_names[i] = chn_dat380[i]['signal']
        chn_names = np.asarray(chn_names)


        #Select desired dictionary from list of dictionaries
        good = np.where(chn_names == 'V34D')[0]
        self.chnspecs = chn_dat380[good[0]]

        """
        self.chnspecs = {"chn":chns[i],
            "tm":tm[i],
            "fs":fs[i],
            "max_sig_pm_mV":max_sig_pm_mV[i],
            "max_sig_pm_mVm":max_sig_pm_mVm[i],
            "boom_length":boom_length[i],
            "gain_desired":gain_desired[i],
            "gain_measured":gain_measured[i],
            "gaindB_desired":gaindB_desired[i],
            "gaindB_measured":gaindB_measured[i],
            "measured_mVm":measured_mVm[i],
            "hpf":hpf[i],
            "lpf":lpf[i],
            "bits":bits[i],
            "coeffa":coeffa[i],
            "coeffb":coeffb[i],
            "polarity":polarity[i]}
        """


if __name__ == '__main__':
    
    #import sys 
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    #from gir_fields_loader import GIRAFF_Fields_Loader
    #import numpy as np

    v12 = GIRAFF_Fields_Loader('381','V12D')
    v1 = GIRAFF_Fields_Loader('381','V1SD')
    hf12 = GIRAFF_Fields_Loader('381','HF12')
    vlf12 = GIRAFF_Fields_Loader('381','VLF12D')
    v1a = GIRAFF_Fields_Loader('381','V1SA')

 
    v12.chnspecs


    v12.plot_gainphase()
    v12.freq_gainphase
    #x1 = self.freq_gainphase[index3dB]

    #Load Steve's save files
    v12dat = v12.load_data()
    v1dat = v1.load_data()
    vlf12dat = vlf12.load_data()
    v1adat = v1a.load_data()
    hf12dat = hf12.load_data()


    #Load Aaron's gain/phase corrected files
    v12datgp = v12.load_data_gainphase_corrected()
    v1datgp = v1.load_data_gainphase_corrected()
    vlf12datgp = vlf12.load_data_gainphase_corrected()







