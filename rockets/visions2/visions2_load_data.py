#Load VISIONS 2 rocket data 



#VISIONS2EEAShiresGSFC.sav
#VISIONS2EIAShiresGSFC.sav
#35040_LFDSP_S1_VLF12_mvm_AaronB.sav


from os.path import dirname, join as pjoin
#import scipy.io as sio 
from scipy.io import readsav 
import matplotlib.pyplot as plt 



path = "/Users/abrenema/Desktop/Research/Rocket_missions/VISIONS2/VISIONS2_data/"
eea_data = readsav(path + 'VISIONS2EEAShiresGSFC.sav')
eia_data = readsav(path + 'VISIONS2EIAShiresGSFC.sav')
vlf12_data = readsav(path + '35040_LFDSP_S1_VLF12_mvm_AaronB.sav')


print("Here")
len(eea_data["eeas_hires_eflux_gsfc"])
len(eea_data["eeas_hires_ftime_gsfc"])
len(eea_data["eeas_hires_energy_gsfc"])
len(eea_data["eeas_hires_pitchnom_gsfc"])

print(vlf12_data["calibrations"])

#data = vlf12_data[""]
#plt.plot()

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

