#Test gainphase curves against theoretical filter response

#RESULTS:

#PROBLEM: V12D channel is falling off like a 1-pole filter. Paulo indicates this should be 4-pole



import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft
import end_load_gainphase as gainphase


#----Channels tested
#filename = "Endurance_Analog 1_V12D_10-10000-100.txt"
#filename = "Endurance_Analog 1_V12A_10-10000-100.txt"
#filename = "Endurance_Analog 1_VLF12D_6-30000-100.txt"
#filename = "Endurance_Analog 1_VLF12A_6-100000-100.txt"
#filename = "Endurance_Analog 1_V1SD_10-10000-100.txt"
filename = "Endurance_Analog 1_V1SA_10-10000-100.txt"


#filename = "Endurance_Analog 1_V3SD_10-10000-100.txt"
#filename = "Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt"



#----Channels not yet tested
#filename = "Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt"


prad, Hmag, f = gainphase.end_load_gainphase(filename)



#Normalize gain curve
maxv = np.nanmax(Hmag)
Hmag_unity = [i/maxv for i in Hmag]

#Revert to dB (with max value 0 dB)
HmagdB = [10*np.log10(i) for i in Hmag_unity]



#deterimine 3dB point (corner freq)
#---First interpolate freqs to a much higher resolution
interp = interp1d(f,HmagdB,kind='cubic', bounds_error=False)
fnew = range(0,int(np.max(f)),1)
HmagdB2 = interp(fnew)

bad = list(map(tuple, np.where(np.isnan(HmagdB2))))
HmagdB2[bad] = 0


#fig,axs = plt.subplots(2)
#axs[0].plot(f,HmagdB,'.')
#axs[1].plot(fnew,HmagdB2,'.')

#plt.plot(f,HmagdB)
#plt.xscale('log')
#plt.xlim(1e3,1e4)
tmp = np.asarray(HmagdB2)
tmp2 = np.asarray(list(fnew))
goo = np.squeeze(list(map(tuple, np.where((tmp <= -3) & (tmp2 > 1000)))))
fc = int(fnew[goo[0]])


#fc = 3300 #filename = "Endurance_Analog 1_V12D_10-10000-100.txt"
#fc = 3800
#fc = 1400 #filename = "Endurance_Analog 1_V12A_10-10000-100.txt"
#fc = 15000 #filename = "Endurance_Analog 1_VLF12D_6-30000-100.txt"
#fc = 7e6   #filename = "Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt"
#fc = 1500  #filename = "Endurance_Analog 1_V1SA_10-10000-100.txt"
#fc = 68000 #filename = "Endurance_Analog 1_VLF12A_6-100000-100.txt"
#fc = 1500  #filename = "Endurance_Analog 1_V1SD_10-10000-100.txt"



#Normalize freq to corner freq
fnorm = [i/fc for i in f]
#fig, axs = plt.subplots(2)
#axs[0].plot(fnorm,HmagdB)
#axs[1].plot(f,Hmag_unity)

#Define filter data points
#f2 = f 
f2 = range(fc,200000,2000)
f2norm = [i/fc for i in f2]
#plt.plot(f2norm)

os = 3  #offset in dB (3 dB point used to define reference frequency fc)
att_1pole = [(-20 * (i-1)/10)-os for i in f2norm]
att_2pole = [(-40 * (i-1)/10)-os for i in f2norm]
att_4pole = [(-80 * (i-1)/10)-os for i in f2norm]



#plt.plot(fnorm,HmagdB,f2norm,att_1pole,'.',f2norm,att_4pole,'.')
plt.plot(f,HmagdB,f2,att_1pole,'.',f2,att_4pole,'.',f2,att_2pole,'.')
plt.xlabel('freq (Hz) \n' + filename + '\n' + 'freq 3dB (measured) = ' + str(fc) + ' Hz')
plt.ylabel('dB from max response\n1pole and 4pole comparison (dots)')
#plt.plot(fnorm,HmagdB,f2norm,att_1pole,'.')
plt.xscale('log')
plt.xlim(1e2,1e5)
plt.ylim(-80,10)

print('Here')












fig, axs = plt.subplots(3)
axs[1].plot(f,Hmag)
#axs[2].plot(f,prad)
#axs[0].set_title('gain/phase; \n fn='+ fn)
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


