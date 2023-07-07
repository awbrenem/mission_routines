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
#   - Loads Endurance data and other info for the Fields instrument
#   
#   T-0: May 11, 2022  01:31:00.0 U.T.
#   Data array original sample rate 30kHz
#   Column 1, time (sec) since T-0, format:f13.7
#   Column 2, VLF12 (mV/m), format:f10.3
#--------------------------------------------------------------------------


class Endurance_Fields_Loader:

    """
    Loadable channels:
    DC: 'V12D','V34D','V13D','V32D','V24D','V41D'
    skins: 'V1SD','V2SD','V3SD','V4SD'
    VLF: 'VLF12D','VLF34D','VLF13D','VLF32D','VLF24D','VLF41D'
    HF: 'HF12','HF34'
    analog: 'V12A','V34A','VLF12A','V1SA','V2SA','V3SA','V4SA'

    Example usage:

    from end_fields_loader import Endurance_Fields_Loader

    v12 = Endurance_Fields_Loader('V12D')
    v12.plot_gainphase()
    v12dat = v12.load_data()

    """


    def __init__(self, chn) -> None:

        self.chn = chn  #Channel specification 
        self._end_channelspecs()
        self._end_load_gainphase()


        if (self.chn == 'V12D') or (self.chn == 'V34D') or (self.chn == 'V13D') or (self.chn == 'V32D') or (self.chn == 'V24D') or (self.chn == 'V41D'): 
            self.type = 'DC'
            filterpole = 4
        if (self.chn == 'V1SD') or (self.chn == 'V2SD') or (self.chn == 'V3SD') or (self.chn == 'V4SD'): 
            self.type = 'skins'
            filterpole = 1
        if (self.chn == 'VLF12D') or (self.chn == 'VLF34D') or (self.chn == 'VLF13D') or (self.chn == 'VLF32D') or (self.chn == 'VLF24D') or (self.chn == 'VLF41D'): 
            self.type = 'VLF'
            filterpole = 2
        if (self.chn == 'HF12') or (self.chn == 'HF34'):
            self.type = 'HF'
            filterpole = 2
        if (self.chn == 'V12A') or (self.chn == 'V34A') or (self.chn == 'VLF12A') or (self.chn == 'V1SA') or (self.chn == 'V2SA') or (self.chn == 'V3SA') or (self.chn == 'V4SA'):
            self.type = 'analog'
            filterpole = "nan"

        self.chnspecs["folder"] = "efield_" + self.type
        #low pass filter pole
        self.chnspecs["lpf_pole"] = filterpole

        return 
    

    def __str__(self):
        return "Endurance " + self.chn + " object"

    #---------------------------------------------------------------------------------
    #Load gain/phase calibrated files. These are produced by end_transfer_function.py.
    #---------------------------------------------------------------------------------

    def load_data_gainphase_corrected(self):

        import pickle 
        path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'


        if self.type == 'DC':
            goo = pickle.load(open(path + 'efield_DC/' + 'Endurance_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
        if self.type == 'skins':
            goo = pickle.load(open(path + 'efield_skins/' + 'Endurance_Analog 1_' + self.chn + '_10-10000-100_gainphase_corrected.pkl', 'rb'))
        if self.type == 'VLF':
            goo = pickle.load(open(path + 'efield_VLF/' + 'Endurance_Analog 1_' + self.chn + '_6-30000-100_gainphase_corrected.pkl', 'rb'))
        
          
        return goo['wf'], goo['tvals']



    #---------------------------------------------------------------------------------
    #Load the data sent to me by Steve Martin. These are gain corrected but not 
    #phase corrected (I use the unity gain to "gain/phase" correct these). 
    #---------------------------------------------------------------------------------

    def load_data(self):
        import numpy as np
        from scipy.io import readsav 

        path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'


        if self.type == 'analog':
            print('no analog data from Steve Martin available to load...RETURNING')
            return -1

        if self.type == 'skins':

            folder = "efield_skins"
            fn = "47001_TM1_LFDSP_S5Skins_V1SD2SD3SD4SD_cal.sav"
            vals = readsav(path + folder + '/' + fn)

            if self.chn == 'V1SD': 
                 wf = vals.dv1s_volts * self.chnspecs["polarity"]
                 t = vals.tv1s
            if self.chn == 'V2SD': 
                 wf = vals.dv2s_volts * self.chnspecs["polarity"]
                 t = vals.tv2s
            if self.chn == 'V3SD': 
                 wf = vals.dv3s_volts * self.chnspecs["polarity"]
                 t = vals.tv3s
            if self.chn == 'V4SD': 
                 wf = vals.dv4s_volts * self.chnspecs["polarity"]
                 t = vals.tv4s

            #Get rid of all times t<0
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good]


        if self.type  == "DC":

            folder = "efield_DC"
            fn = "47001_TM1_LFDSP_S5DCE_DCES5_calibrated.sav"

            vals = readsav(path + folder + '/' + fn)

            t = vals.times
            if self.chn == 'V12D': wf = vals.dv12_mvm * self.chnspecs["polarity"]
            if self.chn == 'V34D': wf = vals.dv34_mvm * self.chnspecs["polarity"]
            if self.chn == 'V13D': wf = vals.dv13_mvm * self.chnspecs["polarity"]
            if self.chn == 'V32D': wf = vals.dv32_mvm * self.chnspecs["polarity"]
            if self.chn == 'V24D': wf = vals.dv24_mvm * self.chnspecs["polarity"]
            if self.chn == 'V41D': wf = vals.dv41_mvm * self.chnspecs["polarity"]

            #Get rid of all times t<0
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good]


        if self.type == 'VLF':

            folder = 'efield_VLF'
            fn = '47001_TM1_LFDSP_S5_VLF_mvm.sav'

            vals = readsav(path + folder + '/' + fn)

            t = vals.tvlf
            if self.chn == 'VLF12D': wf = vals.dvlf12_mvm * self.chnspecs["polarity"]
            if self.chn == 'VLF34D': wf = vals.dvlf34_mvm * self.chnspecs["polarity"]
            if self.chn == 'VLF13D': wf = vals.dvlf13_mvm * self.chnspecs["polarity"]
            if self.chn == 'VLF32D': wf = vals.dvlf32_mvm * self.chnspecs["polarity"]
            if self.chn == 'VLF24D': wf = vals.dvlf24_mvm * self.chnspecs["polarity"]
            #NOTE: VLF41D missing from Steve's file
            #if self.chn == 'VLF41D': wf = vals.dvlf41_mvm

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
            fn = '47001_TM2_0-LOS_S1_HFsnippets_HF1234_method2_mblkskips_4_FFT8400.sav'
            vals = readsav(path + folder + '/' + fn)
            #dict_keys(['afftpow12', 'afftpow34', 'afreq', 'atimesfft', 'fftsize', 'overlap', 'weight', 'samplerate', 'nfreq', 'nfftlines', 'in_file'])

            t = vals.atimesfft
            f = vals.afreq
            if self.chn == 'HF12': wf = vals.afftpow12
            if self.chn == 'HF34': wf = vals.afftpow34

            #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
            good = np.squeeze(np.where(t >= 0.))
            return wf[good], t[good], f[good]

        

        if self.type == 'mag':   #NOTE: NOT YET IMPLEMENTED

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

            #def mag_dc():

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

            """







    #-------------------------------------------------------------------------------------
    #load the gain/phase files for desired channel

    """
    NOTES: Select desired gain/phase files (from Paulo) from the calibration testing.
    From these files Paulo derives a fit (ax + b) that allows calibration from counts to volts.
    These values (a, b) are found in the Endurance channel list document as the yellow boxes in far right column

    See end_gainphase_test.py for comparisons to theoretical behavior
    """
    #-------------------------------------------------------------------------------------

    def _end_load_gainphase(self):

        import numpy as np
        from math import remainder

        path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'

        #Load gain/phase file
        if self.chn == "V1SD": fn = "Endurance_Analog 1_V1SD_10-10000-100.txt"
        elif self.chn == "V2SD": fn = "Endurance_Analog 1_V2SD_10-10000-100.txt"
        elif self.chn == "V3SD": fn = "Endurance_Analog 1_V3SD_10-10000-100.txt"
        elif self.chn == "V4SD": fn = "Endurance_Analog 1_V4SD_10-10000-100.txt"
        elif self.chn == "V1SA": fn = "Endurance_Analog 1_V1SA_10-10000-100.txt"
        elif self.chn == "V2SA": fn = "Endurance_Analog 1_V2SA_10-10000-100.txt"
        elif self.chn == "V3SA": fn = "Endurance_Analog 1_V3SA_10-10000-100.txt"
        elif self.chn == "V4SA": fn = "Endurance_Analog 1_V4SA_10-10000-100.txt"
        elif self.chn == "V12D": fn = "Endurance_Analog 1_V12D_10-10000-100.txt"
        elif self.chn == "V34D": fn = "Endurance_Analog 1_V34D_10-10000-100.txt"
        elif self.chn == "V13D": fn = "Endurance_Analog 1_V13D_10-10000-100.txt"
        elif self.chn == "V32D": fn = "Endurance_Analog 1_V32D_10-10000-100.txt"
        elif self.chn == "V24D": fn = "Endurance_Analog 1_V24D_10-10000-100.txt"
        elif self.chn == "V41D": fn = "Endurance_Analog 1_V41D_10-10000-100.txt"
        elif self.chn == "VLF12D": fn = "Endurance_Analog 1_VLF12D_6-30000-100.txt"
        elif self.chn == "VLF34D": fn = "Endurance_Analog 1_VLF34D_6-30000-100.txt"
        elif self.chn == "VLF13D": fn = "Endurance_Analog 1_VLF13D_6-30000-100.txt"
        elif self.chn == "VLF32D": fn = "Endurance_Analog 1_VLF32D_6-30000-100.txt"
        elif self.chn == "VLF24D": fn = "Endurance_Analog 1_VLF24D_6-30000-100.txt"
        elif self.chn == "VLF41D": fn = "Endurance_Analog 1_VLF41D_6-30000-100.txt"
        elif self.chn == "VLF12A": fn = "Endurance_Analog 1_VLF12A_6-100000-100.txt"
        elif self.chn == "V12A": fn = "Endurance_Analog 1_V12A_10-10000-100.txt"
        elif self.chn == "V34A": fn = "Endurance_Analog 1_V34A_10-10000-100.txt"
        elif self.chn == "HF12": fn = "Endurance_Analog 1_HF12_1000-20000000-100.txt"
        elif self.chn == "HF34": fn = "Endurance_Analog 1_HF34_1000-20000000-100.txt"

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

    def _end_channelspecs(self):

        import numpy as np

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

        #NOTE: the sample freqs in the MEB doc are wrong. Correct values are:
        fsDC = 10000
        fsVLF = 30000
        fsskins = 2000

        chns = ['V12D','V34D','V13D','V32D','V24D','V41D','VLF12D','VLF34D','VLF13D','VLF32D','VLF24D','VLF41D','V1SD','V2SD','V3SD','V4SD','V12A','V34A','VLF12A','V1SA','V2SA','V3SA','V4SA','HF12','HF34']
        tm = [8000,8000,8000,8000,8000,8000,32000,32000,32000,32000,32000,32000,2000,2000,2000,2000,2000,2000,64000,2000,2000,2000,2000,10e6,10e6]
        fs = [fsDC,fsDC,fsDC,fsDC,fsDC,fsDC,fsVLF,fsVLF,fsVLF,fsVLF,fsVLF,fsVLF,fsskins,fsskins,fsskins,fsskins,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),10e6,10e6]
        max_sig_pm_mV = [1251,1251,1251,1251,1251,1251,125,125,125,125,125,125,10000,10000,10000,10000,1251,1251,125,10000,10000,10000,10000,63,63]
        boom_length = [3.212,3.212,2.271227,2.271227,2.271227,2.271227,3.212,3.212,2.271227,2.271227,2.271227,2.271227,float("nan"),float("nan"),float("nan"),float("nan"),3.0,3.0,3.0,float("nan"),float("nan"),float("nan"),float("nan"),3.212,3.212]
        max_sig_pm_mVm = [417.0,417.0,417.0,417.0,417.0,417.0,41.7,41.7,41.7,41.7,41.7,41.7,float("nan"),float("nan"),float("nan"),float("nan"),417.0,417.0,41.7,float("nan"),float("nan"),float("nan"),float("nan"),21.0,21.0]
        gain_desired = [2.0,2.0,2.0,2.0,2.0,2.0,19.98,19.98,19.98,19.98,19.98,19.98,0.25,0.25,0.25,0.25,2.0,2.0,19.98,0.25,0.25,0.25,0.25,39.68,39.68]
        gaindB_desired = [6.0,6.0,6.0,6.0,6.0,6.0,26.0,26.0,26.0,26.0,26.0,26.0,-12.0,-12.0,-12.0,-12.0,6.0,6.0,26.0,-12.0,-12.0,-12.0,-12.0,32.0,32.0]
        gain_measured = [1.991,1.992,1.991,1.992,1.990,1.992,19.763,20.080,19.747,20.016,20.259,19.952,0.203,0.203,0.203,0.203,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]
        gaindB_measured = [5.982,5.987,5.981,5.985,5.979,5.986,25.917,26.055,25.910,26.028,26.132,26.000,-13.854,-13.853,-13.857,-13.849,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]
        measured_mVm = [419,418,419,418,419,418,42,42,42,42,41,42,12321,12320,12326,12314,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]
        hpf = [float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),16,16,16,16,16,16,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),16,float("nan"),float("nan"),float("nan"),float("nan"),3000,3000]
        lpf = [3500,3500,3500,3500,3500,3500,14000,14000,14000,14000,14000,14000,875,875,875,875,875,875,28000,875,875,875,875,4.375e6,4.375e6]
        bits = [18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,14.0,14.0]
        coeffa = [1.2555,1.2549,1.2557,1.2552,-1.256,-1.255,0.1265,0.1245,0.1266,0.1249,0.1234,0.1253,-12.321,-12.32,-12.326,-12.314,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]
        coeffb = [1.1e-03,-4.6288e-05,0.00034567,0.00066907,-0.00074203,-0.00057493,5.5604e-05,5.1386e-05,4.4484e-05,5.4683e-05,5.7427e-05,-0.0005279,0.1449,0.1517,0.1444,0.1192,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]
        polarity = [1,1,1,1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan"),float("nan")]


        #Select desired channel
        i = chns.index(self.chn)



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



if __name__ == '__main__':
    
    #import sys 
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    #from end_fields_loader import Endurance_Fields_Loader
    #import numpy as np

    v12 = Endurance_Fields_Loader('V12D')
    v1 = Endurance_Fields_Loader('V1SD')
    hf12 = Endurance_Fields_Loader('HF12')
    vlf12 = Endurance_Fields_Loader('VLF12D')
    v1a = Endurance_Fields_Loader('V1SA')


    v12.plot_gainphase()
    
    #Load Steve's save files
    v12dat = v12.load_data()
    v1dat = v1.load_data()
    hf12dat = hf12.load_data()
    vlf12dat = vlf12.load_data()
    v1adat = v1a.load_data

    #Load Aaron's gain/phase corrected files
    v12datgp = v12.load_data_gainphase_corrected()
    v1datgp = v1.load_data_gainphase_corrected()
    vlf12datgp = vlf12.load_data_gainphase_corrected()







