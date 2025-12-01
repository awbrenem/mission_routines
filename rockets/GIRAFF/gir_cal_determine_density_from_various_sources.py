# gir_cal_determine_density_from_various_sources.py
# Makes a plot comparing plasma frequency determined from:
# --lower hybrid freq 
# --upper hybrid freq
# --IRI model
# --PFISR data

# Compare this to the cyclotron freq to deterimine over/under dense conditions

# NOTE: currently uses data from the EFW HF. Need to use Marilia's HF data, which will be much more accurate

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/PFISR/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pfisr_load_data as pfld
import iri_read_output as iri
import plasma_params_get_density_from_flhr_freq as dflh
import plot_spectrogram as ps
import iri_vals_to_mission_times as ivt
import numpy as np 
import matplotlib.pyplot as plt
import pickle
from gir_load_fields import GIRAFF_Fields_Loader as GFL
from gir_load_data import load_ephemeris
from datetime import datetime, timezone



#pld = '380'
pld = '381'


#------------------
#load altitude data
ephem  = load_ephemeris(pld)
#plt.plot(ephem['time'], ephem['altitude'])



#------------------
#IRI model (run by Aaron using automated inputs)
filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/IRI/iri_output_36'+pld+'.txt'
iridat = iri.read_iri_output(filename)
fpeiri = 8980.*np.sqrt(iridat['Ne_cm-3'])  #in Hz

iri_timevals = ivt.iri_vals_to_mission_times(iridat, ephem['time'], ephem['altitude'], test=0)



#------------------
#PFISR data

#pftype = '1min'
pftype = '5min'
#pftype = '20min'

if pld == '380':
    fn = '20250209.002_lp_'+pftype+'-fitcal.h5'
if pld == '381':
    fn = '20250202.002_lp_'+pftype+'-fitcal.h5'


#36.380 files
#fn = '20250209.002_lp_1min-fitcal.h5'
#fn = '20250209.002_lp_5min-fitcal.h5'
#fn = '20250209.002_lp_20min-fitcal.h5'

#36.381 files
#fn = '20250202.002_lp_1min-fitcal.h5'
#fn = '20250202.002_lp_5min-fitcal.h5'
#fn = '20250202.002_lp_20min-fitcal.h5' 


filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/pfisr/36.'+pld+'/' + fn

type = 'Ne_Mod'  #Density in m-3
pfisrTimes, pfisrAlts, pfisrNeMod = pfld.pfisr_load_data(filename, type, plot_test=0)

tslice = datetime(2025,2,2,7,7,00, tzinfo=timezone.utc)  #GIRAFF 36.381 launch time
beamnumber = 12
pfisrAltsS, pfisrNeModS = pfld.pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, beamnumber, plot_test=0)

#load altitude data for GIRAFF 36.380
pvst = pfld.pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], plot_test=0)




#------------------
#GIRAFF magnetic field 
magv = GFL(pld,'mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]
fce = 28*Bot



#------------------
#Load lower hybrid vertices as determined from gir_cal_determine_accurate_lower_hybrid_freq.py
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))
#pickle.dump([vertices], open(flhfile,'wb'))

times = vertices[0][:,0]
flh   = vertices[0][:,1]


#------------------
#Load upper hybrid vertices as determined from gir_cal_determine_accurate_upper_hybrid_freq.py

    #Lower cutoff of HF emission band
#fuhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/upper_hybrid_id/GIRAFF_'+pld+'_HF_lower_boundary_vals_byeye.pkl'

    #Middle freq of HF emission band (Kurth method - probably more accurate)
fuhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/upper_hybrid_id/GIRAFF_'+pld+'_HF_middle_vals_byeye.pkl'
vertices = pickle.load(open(fuhfile,'rb'))
timesHF = vertices[0][:,0]
fuhHF  = vertices[0][:,1]



#Interpolate Bo to times of fuh
BoHF = np.interp(timesHF, mag[3], Bo)
fceHF = 28*BoHF
fpeHF = np.sqrt(fuhHF**2 - fceHF**2)



#Interpolate Bo to times of flh
BoI = np.interp(times, mag[3], Bo)
fceI = 28*BoI

nH_ne = [0.03]*len(flh)
nO_ne = [0.97]*len(flh) 



#Interpret the IRI fractional composition values the times of the lower hybrid determination
nH_ne = np.interp(times, ephem['time'], iri_timevals['H+'])
nO_ne = np.interp(times, ephem['time'], iri_timevals['O+'])
nN_ne = np.interp(times, ephem['time'], iri_timevals['N+'])
nHe_ne = np.interp(times, ephem['time'], iri_timevals['He+'])
nO2_ne = np.interp(times, ephem['time'], iri_timevals['O2+'])
nNO_ne = np.interp(times, ephem['time'], iri_timevals['NO+'])

ne1 = dflh.dens_IonMassFractions(flh, fceI, nH_ne, nO_ne, nN_ne, nHe_ne, nO2_ne, nNO_ne)
ne2 = dflh.dens_singleion(flh, BoI, 'H+')
ne3 = dflh.dens_singleion(flh, BoI, 'O+')

fpe1 = 8980.*np.sqrt(ne1)
fpe2 = 8980.*np.sqrt(ne2)
fpe3 = 8980.*np.sqrt(ne3)


fig, axs = plt.subplots(3, 1, figsize=(10, 8))
axs[0].plot(times, fceI, label='fce', color='black')
axs[0].plot(timesHF, fpeHF, label='fpe (upper hybrid)', color='red', linestyle=':')
axs[0].plot(times, fpe1, label='fpe (lower hybrid w/ IRI composition)', linestyle='--')
#axs[0].plot(times, fpe2, label='fpe from H+ only', linestyle='--')
#axs[0].plot(times, fpe3, label='fpe from O+ only', linestyle='--') 
axs[0].plot(iri_timevals['times'], 8980*np.sqrt(iri_timevals['Ne_cm-3']), label='fpe (IRI)', color='blue', linestyle='-.')
#fce (Hz) vs time
axs[0].plot(pvst['times'],8980*np.sqrt(pvst['pfisrdat_interp']/1e6), label='fpe (PFISR LP '+pftype+')', color='purple', linestyle='-.')
axs[0].set_ylim(0,5e6)
axs[0].set_ylabel('Frequency (Hz)')
axs[0].set_title('GIRAFF 36.'+pld+' Comparison b/t plasma and cyclotron freq\ngir_cal_determine_density_from_various_sources.py', fontsize=8)
#axs[0].set_xticklabels([])

# top x-axis showing rocket altitude corresponding to time
ax_top = axs[0].twiny()
def _update_top_xlim(ax):
    # choose 10 tick positions across current x-limits and label with interpolated altitude
    xmin, xmax = ax.get_xlim()
    ticks = np.linspace(xmin, xmax, 10)
    alt_ticks = np.interp(ticks, ephem['time'], ephem['altitude'])
    ax_top.set_xticks(ticks)
    ax_top.set_xticklabels([f"{alt:.0f}" for alt in alt_ticks])
    ax_top.set_xlim(xmin, xmax)
    ax_top.set_xlabel('Altitude')
# initialize and connect the update so labels follow any later xlim changes
_update_top_xlim(axs[0])
axs[0].callbacks.connect('xlim_changed', _update_top_xlim)


axs[1].plot(ephem['time'], iri_timevals['O+'], label='O+ IRI', color='blue', linestyle='--')
axs[1].plot(ephem['time'], iri_timevals['N+'], label='N+ IRI', color='magenta', linestyle='--')
axs[1].plot(ephem['time'], iri_timevals['H+'], label='H+ IRI', color='red', linestyle='--')
axs[1].plot(ephem['time'], iri_timevals['He+'], label='He+ IRI', color='green', linestyle='--')
axs[1].plot(ephem['time'], iri_timevals['O2+'], label='O2+ IRI', color='orange', linestyle='--')
axs[1].plot(ephem['time'], iri_timevals['NO+'], label='NO+ IRI', color='brown', linestyle='--')
#axs[1].set_xticklabels([])

axs[2].plot(ephem['time'], iri_timevals['Tn_K'], label='Tn IRI', color='black', linestyle='--')
axs[2].plot(ephem['time'], iri_timevals['Ti_K'], label='Ti IRI', color='blue', linestyle='-')
axs[2].plot(ephem['time'], iri_timevals['Te_K'], label='Te IRI', color='red', linestyle='--')
axs[2].set_ylabel('Tn, Ti, Te (K)')
axs[2].set_xlabel('Time (s)')

for i in range(3):
    axs[i].set_xlim(100,500)
    axs[i].grid()
    axs[i].legend(fontsize=6, loc='upper right')
plt.show()
