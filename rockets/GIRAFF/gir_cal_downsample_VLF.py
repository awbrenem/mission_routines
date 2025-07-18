
#Downsample the GIRAFF VLF data to be at a lower cadence for use with waves < 50 kHz (nearly all the interesting waves).
#Not doing this means the waveform files are HUGE and difficult to use. 
#NOTE: VLF12 and VLF34 are 800 kS/sec while the diagonals are 100 kS/sec


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency as bp
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from scipy import signal
import numpy as np 
import pickle
import matplotlib.pyplot as plt
import plot_spectrogram as ps


#----------SETUP---------------

testing = 0   #If set, will make plots for testing


#Load E-fields data
pld = '381'
chn = 'VLF34DF'   #full resolution data
vdat = GFL(pld,chn)
wf_orig, tdat_orig = vdat.load_data()

srt_old = 1/(tdat_orig[1]-tdat_orig[0])




#Determine number of points to resample to. 
nyq = 50000. #Max useable frequency 


#----------END SETUP---------------



#New sample rate - make 10% higher than Nyquist to avoid aliasing
srt_new = 2*nyq + 2*nyq*0.1 




#original and renamed files
if chn[0:3] == 'VLF':
    chngoo = chn[:-1]
else: 
    chngoo = chn

if pld == '381':
    fnorig = '36' + pld + '_GIRAFF_TM2_LFDSP_IT_'+chngoo+'_mvm'
else:
    fnorig = '36' + pld + '_GIRAFF_TM2_LFDSP_'+chngoo+'_mvm'
fnnew = fnorig + '_downsampled_50kHz'






#Resample in chunks so that we don't have to pass huge arrays (very slow).
chunksz = int(2**22)  #number of points for each chunk (power of 2 may lead to faster computation)
nchunks = int(len(wf_orig) / chunksz)  
nsamp = int(np.round(srt_new * (tdat_orig[int(chunksz)] - tdat_orig[0]))) #samples that will be in each chunk after resampling 


#Filter to new sample rate to avoid aliasing, which can affect the final signal amplitude. 
valsbp = bp.butter_lowpass_filter(wf_orig, nyq, srt_old, order=5)



#Downsample each chunk
wf = np.zeros((nsamp, nchunks))
tv = np.zeros((nsamp, nchunks))

for i in range(nchunks):

    y, x = signal.resample(valsbp[i*chunksz:(i+1)*chunksz], nsamp, t=tdat_orig[i*chunksz:(i+1)*chunksz], window='hann')
    wf[:,i] = y 
    tv[:,i] = x


    if testing:

        #Test plot to see that the phase is not altered. 
        #---> Result: phase is NOT altered. 
        #---> However, amplitude is decreased. 
        srt1 = int(1/(x[1]-x[0]))  #exact downsampled sample rate
        srt2 = int(1/(tdat_orig[1]-tdat_orig[0]))

        #dtplot = 0.001  #sec 
        dtplot = 0.001  #sec 
        npts1 = int(srt1*dtplot)
        npts2 = int(srt2*dtplot)


        #plt.plot(tdat_orig[i*chunksz:(i*chunksz)+npts2],wf_orig[i*chunksz:(i*chunksz)+npts2])
        plt.plot(tdat_orig[i*chunksz:(i*chunksz)+npts2],valsbp[i*chunksz:(i*chunksz)+npts2],'*')
        plt.plot(x[0:npts1],y[0:npts1],'*')
        plt.xlim(tdat_orig[i*chunksz],tdat_orig[i*chunksz + npts2])
        #plt.ylim(-5,5)
        plt.show()
        print('h')



        #FFT this chunk
        fs = 1/ (tdat_orig[1] - tdat_orig[0])
        wftst = wf_orig[i*chunksz:(i*chunksz)+100*npts2]
        wftst = valsbp[i*chunksz:(i*chunksz)+100*npts2]
        nps = 1024
        fspecx, tspecx, powercAx = signal.spectrogram(wftst, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
        ps.plot_spectrogram(tspecx,fspecx,np.abs(powercAx),yr=[0,400000], yscale='linear',vr=[-70,-20])

        print('h')








#Reshape the chunked arrays to an array of size [n]
wf = wf.reshape(nsamp*nchunks, order='F')
times = tv.reshape(nsamp*nchunks, order='F')



notes = {'VLF waveform data from ' + fnorig + ' downsampled to 50kHz from gir_cal_downsample_VLF.py'}


path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/efield_VLF/' + fnnew + '.pkl'
pickle.dump([times, wf, notes],open(path,'wb'))





"""
#Test the results
import plot_spectrogram as ps
import matplotlib.pyplot as plt

plt.plot(times)
plt.plot(times,wf)
#plt.plot(tdat,wf12)

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tv2, wf2, fs_thres=0.1, nfft=nfft, noverlap=2)

vr=[-80,-50]
xr = [100,550]
yr = [0,50000]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=[0,1400],xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='log',yr=[300,50000],xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")

"""


