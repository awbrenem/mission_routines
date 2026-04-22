#gir_analysis_DC_Efield.py


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import numpy as np 
import matplotlib.pyplot as plt
from scipy.io import readsav
from matplotlib.ticker import MultipleLocator
import gir_load_data as gld



#Get altitude data
ephem = gld.load_ephemeris('380')
alts380 = ephem['altitude']
time_ephem380 = ephem['time']

ephem = gld.load_ephemeris('381')
alts381 = ephem['altitude']
time_ephem381 = ephem['time']



E380 = gld.load_EDCdespun_Henry('380')
E381 = gld.load_EDCdespun_Henry('381')

#-------------------------------------------
#Version 1 has y-axes adjusted for each plot
#-------------------------------------------

xlim = [100,520]

fig, axs = plt.subplots(2, 2, figsize=(10, 6))
plt.title('GIRAFF DC E-field Data')
axs[0,0].plot(E380['time'], E380['Ezon'], label='Ezon')
axs[0,0].set_xlabel('Time (sec)', fontsize=6)
axs[0,0].set_ylabel('Ezon (mV/m)', fontsize=6)
axs[0,0].set_title('DC E-field Data\n380', fontsize=6)
axs[0,0].tick_params(labelsize=6)
axs[0,0].set_ylim(-15, 20)
axs[0,0].legend(fontsize=6)
axs[1,0].plot(E380['time'], E380['Emerid'], label='Emerid', color='orange')
axs[1,0].set_xlabel('Time (sec)', fontsize=6)
axs[1,0].set_ylabel('Emerid (mV/m)', fontsize=6)
axs[1,0].tick_params(labelsize=6)
axs[1,0].set_ylim(-100, 50)
axs[1,0].legend(fontsize=6)

axs[0,1].plot(E381['time'], E381['Ezon'], label='Ezon')
axs[0,1].set_xlabel('Time (sec)', fontsize=6)
axs[0,1].set_ylabel('Ezon (mV/m)', fontsize=6)
axs[0,1].set_title('DC E-field Data\n381', fontsize=6)
axs[0,1].set_ylim(-50, -10)
axs[0,1].tick_params(labelsize=6)
axs[0,1].legend(fontsize=6)
axs[1,1].plot(E381['time'], E381['Emerid'], label='Emerid', color='orange')
axs[1,1].set_xlabel('Time (sec)', fontsize=6)
axs[1,1].set_ylabel('Emerid (mV/m)', fontsize=6)
axs[1,1].set_ylim(-40, 30)
axs[1,1].tick_params(labelsize=6)
axs[1,1].legend(fontsize=6)

for i in range(2):
    for j in range(2):
        axs[i,j].set_xlim(xlim)
        axs[i,j].xaxis.set_major_locator(MultipleLocator(25))
        axs[i,j].yaxis.set_major_locator(MultipleLocator(10))
        axs[i,j].grid(True, alpha=0.3)

for i in range(2):
    ax2 = axs[i,0].twinx()
    ax2.plot(time_ephem380, alts380, color='blue', linewidth=2)
    ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
    ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis='y', labelleft=True, labelright=False)
    ax2.spines['left'].set_position(('outward', 60))

    ax2 = axs[i,1].twinx()
    ax2.plot(time_ephem381, alts381, color='blue', linewidth=2)
    ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
    ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis='y', labelleft=True, labelright=False)
    ax2.spines['left'].set_position(('outward', 60))

fig.tight_layout(pad=1)
fig.subplots_adjust(left=0.2)
plt.show()



#-------------------------------------------
#Version 2 has y-axes adjusted for comparison b/t 380 and 381
#-------------------------------------------

xlim = [100,520]

fig, axs = plt.subplots(2, 2, figsize=(10, 6))
plt.title('GIRAFF DC E-field Data')
axs[0,0].plot(E380['time'], E380['Ezon'], label='Ezon')
axs[0,0].set_xlabel('Time (sec)', fontsize=6)
axs[0,0].set_ylabel('Ezon (mV/m)', fontsize=6)
axs[0,0].set_title('DC E-field Data\n380', fontsize=6)
axs[0,0].tick_params(labelsize=6)
axs[0,0].set_ylim(-50, 20)
axs[0,0].legend(fontsize=6)
axs[1,0].plot(E380['time'], E380['Emerid'], label='Emerid', color='orange')
axs[1,0].set_xlabel('Time (sec)', fontsize=6)
axs[1,0].set_ylabel('Emerid (mV/m)', fontsize=6)
axs[1,0].tick_params(labelsize=6)
axs[1,0].set_ylim(-100, 50)
axs[1,0].legend(fontsize=6)

axs[0,1].plot(E381['time'], E381['Ezon'], label='Ezon')
axs[0,1].set_xlabel('Time (sec)', fontsize=6)
axs[0,1].set_ylabel('Ezon (mV/m)', fontsize=6)
axs[0,1].set_title('DC E-field Data\n381', fontsize=6)
axs[0,1].set_ylim(-50, 20)
axs[0,1].tick_params(labelsize=6)
axs[0,1].legend(fontsize=6)
axs[1,1].plot(E381['time'], E381['Emerid'], label='Emerid', color='orange')
axs[1,1].set_xlabel('Time (sec)', fontsize=6)
axs[1,1].set_ylabel('Emerid (mV/m)', fontsize=6)
axs[1,1].set_ylim(-100, 50)
axs[1,1].tick_params(labelsize=6)
axs[1,1].legend(fontsize=6)

for i in range(2):
    for j in range(2):
        axs[i,j].set_xlim(xlim)
        axs[i,j].xaxis.set_major_locator(MultipleLocator(25))
        axs[i,j].yaxis.set_major_locator(MultipleLocator(10))
        axs[i,j].grid(True, alpha=0.3)

for i in range(2):
    for j in range(2):
        ax2 = axs[i,j].twinx()
        ax2.plot(time_ephem, alts, color='blue', linewidth=2)
        ax2.set_ylabel('Altitude (km)', color='blue', fontsize=6)
        ax2.tick_params(axis='y', labelcolor='blue', labelsize=6)
        ax2.yaxis.set_label_position("left")
        ax2.tick_params(axis='y', labelleft=True, labelright=False)
        ax2.spines['left'].set_position(('outward', 60))


axs[0,0].set_ylim(-50, 20)

fig.tight_layout(pad=1)
fig.subplots_adjust(left=0.2)
plt.tight_layout()
plt.show()







#    v12 = GIRAFF_Fields_Loader('V12D')
#    v12.plot_gainphase('381')
#    v12dat = v12.load_data()



"""
import xarray as xr

folder = "efield_DC"
fn = "Giraff_380_efield_DC_not_gp_calibrated.netcdf"
path = "/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/"

# Load the netcdf file
ds = xr.open_dataset(path + folder + '/' + fn)

time, E12, E34 

pld = ''

plt.plot(ds['time'], ds['E12'], label='E12')


ef = GFL(pld).load_fields()
ef.load_efield()

"""


#Load the Zesta ENU data from Alex Hoffman
#zmag = gld.load_zesta_mag(pld)

print('h')
