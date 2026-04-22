"""
Compare the density values from all the different PFISR beams at various time cadences

I've been told by the PFISR team that 20 min integrations are the best 
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/PFISR/')
import numpy as np 
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import gir_load_data
import pfisr_load_data as pfld
from datetime import datetime, timezone
from matplotlib.ticker import MultipleLocator



pld = '380'




#------------------
#load altitude data
ephem  = gir_load_data.load_ephemeris(pld)

alts = ephem['altitude']
time_ephem = ephem['time']


#Load PFISR data for comparison
#pftype = '1min'
#pftype = '5min'
pftype = '20min'

if pld == '380':
    fn = '20250209.002_lp_'+pftype+'-fitcal.h5'
    tslice = datetime(2025,2,9,8,35,00, tzinfo=timezone.utc)  #GIRAFF 36.380 launch time
if pld == '381':
    fn = '20250202.002_lp_'+pftype+'-fitcal.h5'
    tslice = datetime(2025,2,2,7,7,00, tzinfo=timezone.utc)  #GIRAFF 36.381 launch time



filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/pfisr/36.'+pld+'/' + fn

type = 'Ne_Mod'  #Density in m-3
pfisrTimes, pfisrAlts, pfisrNeMod = pfld.pfisr_load_data(filename, type, plot_test=0)

pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 1, plot_test=0)
pvstBeam1 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 2, plot_test=0)
pvstBeam2 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 3, plot_test=0)
pvstBeam3 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 4, plot_test=0)
pvstBeam4 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 5, plot_test=0)
pvstBeam5 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 6, plot_test=0)
pvstBeam6 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 7, plot_test=0)
pvstBeam7 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 8, plot_test=0)
pvstBeam8 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 9, plot_test=0)
pvstBeam9 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 10, plot_test=0)
pvstBeam10 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 11, plot_test=0)
pvstBeam11 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 12, plot_test=0)
pvstBeam12 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, 13, plot_test=0)
pvstBeam13 = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)






#----------------------------------------------------
#Marilia's HF spectral data (E||) for comparison
#----------------------------------------------------

times, freqs, power = gir_load_data.load_marilia_HF(pld)


#-----------------------------------------------------
#Final plot
#-----------------------------------------------------


title = '36.'+pld+ ' (gir_cal_pfisr_compare_beams.py)\n'+str(tslice)+'\n'+pftype+' beam integration'

colors = ['blue','orange','green','red','purple','maroon','magenta','navy','gold','cyan','darkblue','orangered','teal']



fig = plt.figure(figsize=(13,8))
fig.suptitle(title,fontsize=8)
yr = [1e6, 5e6]
xr = [100, 520]
vr = [-95,-90]

ax1 = plt.subplot(2,2,1)
ax1.plot(pvstBeam2['times'],8980*np.sqrt(pvstBeam2['pfisrdat_interp']/1e6), label='Beam 2 (VERTICAL)', color=colors[1], linewidth=4)
ax1.plot(pvstBeam1['times'],8980*np.sqrt(pvstBeam1['pfisrdat_interp']/1e6), label='Beam 1', color=colors[0], linestyle='-.',linewidth=2)
ax1.set_ylim(yr)
ax1.set_xlim(xr)
ax1.set_ylabel('fpe (Hz)')
ax1.legend()
ax1.set_title('Near/Over Poker')

ax2 = plt.subplot(2,2,2)
ax2.plot(pvstBeam5['times'],8980*np.sqrt(pvstBeam5['pfisrdat_interp']/1e6), label='Beam 5', color=colors[4], linestyle='-.',linewidth=2)
ax2.plot(pvstBeam8['times'],8980*np.sqrt(pvstBeam8['pfisrdat_interp']/1e6), label='Beam 8', color=colors[7], linestyle='-.',linewidth=2)
ax2.plot(pvstBeam11['times'],8980*np.sqrt(pvstBeam11['pfisrdat_interp']/1e6), label='Beam 11', color=colors[10], linestyle='-.',linewidth=2)
ax2.plot(pvstBeam3['times'],8980*np.sqrt(pvstBeam3['pfisrdat_interp']/1e6), label='Beam 3', color=colors[2], linestyle='-.',linewidth=2)
ax2.plot(pvstBeam4['times'],8980*np.sqrt(pvstBeam4['pfisrdat_interp']/1e6), label='Beam 4', color=colors[3], linestyle='-.',linewidth=2)
ax2.set_ylim(yr)
ax2.set_xlim(xr)
ax2.set_ylabel('fpe (Hz)')
ax2.set_title('North of Poker')
ax2.legend()

ax3 = plt.subplot(2,2,3)
ax3.plot(pvstBeam7['times'],8980*np.sqrt(pvstBeam7['pfisrdat_interp']/1e6), label='Beam 7', color=colors[6], linestyle='-.',linewidth=2)
ax3.plot(pvstBeam10['times'],8980*np.sqrt(pvstBeam10['pfisrdat_interp']/1e6), label='Beam 10', color=colors[9], linestyle='-.',linewidth=2)
ax3.plot(pvstBeam13['times'],8980*np.sqrt(pvstBeam13['pfisrdat_interp']/1e6), label='Beam 13', color=colors[12], linestyle='-.',linewidth=2)
ps.plot_spectrogram(times,freqs,power,vr=vr,yscale='linear',xr=xr,yr=yr,ylabel="HF E|| (Marilia)\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=ax3,title=title,colorbar=0)
ax3.set_ylim(yr)
ax3.set_xlim(xr)
ax3.set_ylabel('fpe (Hz)')
ax3.set_title('West/North of Poker')
ax3.legend()
ax3.set_xlabel('time (sec)')    

ax4 = plt.subplot(2,2,4)
ax4.plot(pvstBeam6['times'],8980*np.sqrt(pvstBeam6['pfisrdat_interp']/1e6), label='Beam 6', color=colors[5], linestyle='-.',linewidth=2)
ax4.plot(pvstBeam9['times'],8980*np.sqrt(pvstBeam9['pfisrdat_interp']/1e6), label='Beam 9', color=colors[8], linestyle='-.',linewidth=2)
ax4.plot(pvstBeam12['times'],8980*np.sqrt(pvstBeam12['pfisrdat_interp']/1e6), label='Beam 12', color=colors[11], linestyle='-.',linewidth=2)
ax4.set_ylim(yr)
ax4.set_xlim(xr)
ax4.set_title('East/North of Poker')
ax4.set_xlabel('time (sec)')    
ax4.legend()

# Apply grid to all subplots
for ax in [ax1, ax2, ax3, ax4]:
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(.5e6))
    ax.grid(True, alpha=0.3)

    ax2 = ax.twinx()
    ax2.plot(time_ephem, alts, color='blue', linewidth=2)
    ax2.set_xlim(xr)
    ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
    ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis='y', labelleft=True, labelright=False)
    ax2.spines['left'].set_position(('outward', 60))



fig.tight_layout(pad=1)
fig.subplots_adjust(left=0.2)
#plt.show()




