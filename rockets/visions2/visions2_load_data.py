#Load VISIONS 2 rocket data 



#VISIONS2EEAShiresGSFC.sav
#VISIONS2EIAShiresGSFC.sav
#35040_LFDSP_S1_VLF12_mvm_AaronB.sav


from os.path import dirname, join as pjoin
#import scipy.io as sio 
from scipy.io import readsav 
import matplotlib.pyplot as plt 
from matplotlib import interactive 
interactive(True)
import numpy as np

from scipy import signal
from scipy.fft import fftshift



path = "/Users/abrenema/Desktop/Research/Rocket_missions/VISIONS2/VISIONS2_data/"
eea_data = readsav(path + 'VISIONS2EEAShiresGSFC.sav')
eia_data = readsav(path + 'VISIONS2EIAShiresGSFC.sav')
vlf12_data = readsav(path + '35040_LFDSP_S1_VLF12_mvm_AaronB.sav') #Low flyer (35.040)


#np.ndarray.shape(eea_data['eeas_hires_eflux_gsfc'])

eea_data.keys()

#Mag data (Bart = Bartington mag; Bon = Bonalsky (Fluxgate))
mag_bart_nwu = readsav(path + '35040_Bart_mag_NWU.sav') 
mag_bart_xyz = readsav(path + '35040_Bart_mag_XYZ.sav')
mag_bon_nwu_lmc = readsav(path + '35040_BON_mag_model_NWU_LMC.sav')
mag_bon_xyz = readsav(path + '35040_Bonalsky_mag_XYZ.sav')

bo_bon = np.sqrt(mag_bon_xyz['xmag']**2 + mag_bon_xyz['ymag']**2 + mag_bon_xyz['zmag']**2)
fce_bon = [28.*bo_bon[i] for i in range(len(bo_bon))]
fci_bon = [fce_bon[i]/1836. for i in range(len(fce_bon))]
flh_bon = [np.sqrt(fce_bon[i]*fci_bon[i]) for i in range(len(fce_bon))]
times_bon = mag_bon_xyz['t']


print("Here")
len(eea_data["eeas_hires_eflux_gsfc"])
len(eea_data["eeas_hires_ftime_gsfc"])
len(eea_data["eeas_hires_energy_gsfc"])
len(eea_data["eeas_hires_pitchnom_gsfc"])

np.shape(eea_data["eeas_hires_eflux_gsfc"]) #(12480, 49, 20)
np.shape(eea_data["eeas_hires_ftime_gsfc"]) #12480
np.shape(eea_data["eeas_hires_energy_gsfc"]) #49 -> 3 to 30000 eV 
np.shape(eea_data["eeas_hires_pitchnom_gsfc"]) #20 -> -180 to 180 alternating



print(vlf12_data["calibrations"])

vlf12 = vlf12_data["dvlf12"]
times = vlf12_data["tvlf12"]
#plt.plot(times, vlf12)
#plt.show()
#plt.close()




#sampling freq
sr = [1/(times[i+1]-times[i]) for i in range(times.size-1)]
fs = np.mean(sr)


#npfft = 256
#specfreqs, spectimes, Sxx = signal.spectrogram(vlf12, fs, nfft=npfft)
specfreqs, spectimes, Sxx = signal.spectrogram(vlf12, fs, nperseg=1024)


#Sxx2, specfreqs2, spectimes2, im = plt.specgram(vlf12, Fs=fs, NFFT=npfft)


#(257, 41452)



#interpolate mag times to spec times 
fci_bon_interp = np.interp(spectimes, times_bon, fci_bon)
flh_bon_interp = np.interp(spectimes, times_bon, flh_bon)

#Turn power into dB
SxxdB = 10*np.log10(Sxx)



#Recreate Rob's poster plot of EEAs electrons from 3-30 eV for 0 deg pitch angle bin
# for times of 326 - 344 sec from launch

tmpt = eea_data["eeas_hires_ftime_gsfc"]
tmpe = eea_data["eeas_hires_energy_gsfc"]
tmpv = eea_data["eeas_hires_eflux_gsfc"][:,:,0]

#Turn flux into dB 
tmpvdB = 20.*np.log10(tmpv)

goodt = np.squeeze(np.where((tmpt > 326) & (tmpt < 344)))

"""
plt.pcolormesh(tmpt[goodt], tmpe, np.transpose(tmpv[goodt,:]), shading="gouraud")
plt.yscale('log')
plt.ylim(1, 30000)
plt.title("EAAs electrons")
plt.ylabel("Energy [eV]")
plt.xlabel("MET (sec)")
"""





fig, axs = plt.subplots(2)

fig.suptitle('VISIONS-2')

#pcm = axs[0].pcolormesh(spectimes, specfreqs, SxxdB, shading='gouraud', cmap='turbo', vmin=np.min(SxxdB)/2., vmax=np.max(SxxdB)/10)
pcm = axs[0].pcolormesh(spectimes, specfreqs, SxxdB, shading='gouraud', cmap='turbo', vmin=-70, vmax=-20)
axs[0].plot(spectimes, fci_bon_interp)
axs[0].plot(spectimes, flh_bon_interp)
axs[0].set_ylabel('Frequency [Hz]')
axs[0].set_xlabel('Time [sec]')
axs[0].set_yscale('log')
axs[0].set_ylim(np.min(specfreqs), np.max(specfreqs))
fig.colorbar(pcm, label="VLF12 [mV/m]^2/Hz dB", ax=axs[0])

pcm1 = axs[1].pcolormesh(tmpt, tmpe, np.transpose(tmpvdB), shading='gouraud', cmap='turbo')
axs[1].set_yscale('log')
axs[1].set_ylim(1, 30000)
axs[1].set_title("EAAs electrons")
axs[1].set_ylabel("Energy [eV]")
axs[1].set_xlabel("MET (sec)")
fig.colorbar(pcm1, label="Electronflux", ax=axs[1])


"""
pcm2 = axs[1].pcolormesh(tmpt, tmpe, np.transpose(tmpv[:,:]), shading="gouraud")
axs[1].set_yscale('log')
axs[1].set_ylim(1, 30000)
axs[1].set_title("EAAs electrons")
axs[1].set_ylabel("Energy [eV]")
axs[1].set_xlabel("MET (sec)")
fig.colorbar(pcm2, label="test2", ax=axs[1])
"""


#for i in axs: i.set_xlim(0, np.max(spectimes))
#for i in axs: i.set_xlim(328, 342)
for i in axs: i.set_xlim(200, 700)



"""
plt.pcolormesh(spectimes, specfreqs, SxxdB, shading='gouraud', vmin=-16, vmax=-3)
plt.plot(spectimes, fci_bon_interp)
plt.plot(spectimes, flh_bon_interp)
plt.colorbar(label="test")
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()
plt.close()
"""




print("Here")





#plt.show()
#plt.close()







"""
 fig, ax = plt.subplots(1)
    ax[0].set_title(str(len(df_surviving_list_tmp))+" events: L="+str(lrange[0])+"-"+str(lrange[1])+"  MLT="+str(mltrange[0])+"-"+str(mltrange[1])+"  dL="+str(dlrange[0])+"-"+str(dlrange[1])+"  dMLT="+str(dmltrange[0])+"-"+str(dmltrange[0]))
    ax[0].scatter(df_surviving_list_tmp[xkey], df_surviving_list_tmp[ykey], color='lightgray')
    ax[0].set_xlabel(xkey)
    ax[0].set_ylabel(ykey)
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    ax[0].set_xlim([xmin, xmax])
    ax[0].set_ylim([ymin, ymax])
""" 

