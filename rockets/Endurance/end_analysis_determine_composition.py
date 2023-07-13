#Endurance - identify harmonics of the Bernstein waves

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
import end_load_data
import end_functions as end
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
import plasma_params_get_density_from_flhr_freq as dflh
import plasma_params_get_flhr_freq as dflh2
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps
from scipy import signal
from scipy.interpolate import interp1d
import numpy as np 
import matplotlib.pyplot as plt
from astropy import units as u  
import pickle


#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]

"""Enable auto module reloading"""
#%load_ext autoreload
#%autoreload 2


"""Load E-field VLF data"""
evlf = end_load_data.efield_vlf()
print(evlf.keys())

"""Load mag data"""
mag = end_load_data.mag_dc()





#Get cyclotron freqs
Bo = signal.decimate(mag.Bmag, 10)
tvals = signal.decimate(mag.tsec,10)
fce = [28*i for i in Bo]
fcH = [i/1836.15 for i in fce]
fcO = [i/(1836.15*15) for i in fce]
fcN = [i/(1836.15*14) for i in fce]
fcHe = [i/(1836.15*4) for i in fce]



#Use Langmuir probe density values to determine fpe and fuh and then O+ to H+ fractions


#----------NEED BETTER VALUES (these are extremely rough)-----------
# Langmuir probe density estimate from plot
tLP = [135,200,250,300,350,400,450,500,550,600,650,700,750,800,840,870]
dLP = [2.75,1.5,0.75,0.4,0.3,0.25,0.2,0.2,0.2,0.2,0.4,0.6,1.2,2,2.4,1]
dLP = [i*1e5 for i in dLP]


#----------NEED BETTER VALUES-----------
#flh identification 
tlh = [117,139,167,206,264,335,420,479,574,666,759,817,850,868,879]
flh = [7538,8241,7882,7338,6956,6613,6476,6438,6632,6997,7426,7858,7515,5733,4711]
fpeLP = [8980.*np.sqrt(i) for i in dLP]

#Interpret cyclotron freq values to times of flh determination
interp = interp1d(tvals,fce,kind='cubic', bounds_error=False)
fce2 = interp(tlh)
interp = interp1d(tvals,fcH,kind='cubic', bounds_error=False)
fcH2 = interp(tlh)
interp = interp1d(tvals,fcO,kind='cubic', bounds_error=False)
fcO2 = interp(tlh)
interp = interp1d(tvals,Bo,kind='cubic', bounds_error=False)
Bo2 = interp(tlh)


#Interpolate Langmuir probe values to flh identification times 
interpd = interp1d(tLP,dLP,kind='cubic', bounds_error=False)
dLP2 = interpd(tlh)

#Add units
fce_u = fce2 * u.Hz
Bo_u = Bo2 * u.nT
flh_u = flh * u.Hz

#Select %H+ for each time of flh identification.
#Adjust so that the determined density matches the Langmuir probe values for ne. 
nH_ne = [0.0001,0.0073,0.0064,0.0035,0.011,0.04,0.05,0.05,0.05,0.04,0.015,0.011,0.00004,0.000,0.5] * u.dimensionless_unscaled 
nO_ne = [1-i for i in nH_ne] * u.dimensionless_unscaled 

ne1 = dflh.dens_IonMassFractions(flh_u, fce_u, nH_ne, nO_ne)
ne2 = dflh.dens_singleion(flh_u, Bo_u, 'H+')
ne3 = dflh.dens_singleion(flh_u, Bo_u, 'O+')

ne1 = [i.value if not np.isnan(i) else 0 for i in ne1]
ne2 = [i.value if not np.isnan(i) else 0 for i in ne2]
ne3 = [i.value if not np.isnan(i) else 0 for i in ne3]



plt.plot(tlh,dLP2,'.',tlh,ne1,'x')
plt.show()
print("here")


info = ["From end_determine_composition.py." +
        "Use Langmuir probe derived density values along with id of the flhr line to determine plasma composition." + 
        "ne_flhID_HpOp are the density values using a mixed plasma of H+ and O+ with the fractional ratios relative to ne of nH_ne and nO_ne." + 
        "ne_flhID_Hp and ne_flhID_Op are the same but for 100% H+ and O+ plasmas."]


dict_fin = {'times':tlh,'ne_langmuirprobe':dLP,
            'ne_flhID_HpOp':ne1,
            'nH_ne':nH_ne,'nO_ne':nO_ne,
            'ne_flhID_Hp':ne2,
            'ne_flhID_Op':ne3,
            'info':info}

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/plasma_composition.pkl'
pickle.dump(dict_fin, open(path,'wb'))



#----------------------------------------------
#Rough altitude comparison with IRI model O+ % 


import pandas as pd 

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'
folder = 'ephemeris'
fn = 'Endurance_GPS_velocity_position_altitude.csv'

header = ['time','xvel','yvel','zvel','lat','long','alt']
ephem =  pd.read_csv(path + folder + '/' + fn, skiprows=1, names=header)


#IRI O+ % 
alt_iri = [100.00,150.00,200.00,250.00,300.00,350.00,400.00,450.00,500.00,550.00,600.00,650.00,700.00,750.00,800.00,850.00,900.00]
Op_iri = [0.0,1.2,24.7,76.1,91.3,96.8,96.2,95.5,94.7,93.8,93.1,92.5,91.7,90.8,89.7,88.4,87.0]
Hp_iri = [0.0,0.0,0.0,0.0,0.0,0.1,0.2,0.3,0.4,0.7,1.0,1.3,1.7,2.2,2.9,3.6,4.4]
Op_iri = [i/100 for i in Op_iri]
Hp_iri = [i/100 for i in Hp_iri]

#Associate rocket times with altitudes so I can compare to IRI
altv = np.asarray(ephem.alt)
alts = np.empty(len(tlh))
for i in range(len(tlh)):
    good = np.where(ephem.time >= tlh[i])
    good2 = good[0][0]
    alts[i] = altv[good2]


plt.plot(alt_iri, Op_iri, alts, nO_ne)
plt.show()


