#Endurance - determine ion composition by requiring that the SLP density matches 
#the density determined from the by-eye lower hybrid line identification (from end_cal_determine_accurate_lower_hybrid_freq.py) 

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
#import end_load_data
import end_functions as end
import plasma_params_get_density_from_flhr_freq as dflh
import plasma_params_get_flhr_freq as dflh2
import plot_spectrogram as ps
from scipy import signal
from scipy.interpolate import interp1d
import numpy as np 
import matplotlib.pyplot as plt
from astropy import units as u  
import pickle
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
import pyIGRF
#import iri2016


#iri2016.IRI(’datetime’, [alt.min, alt.max, alt.step], latitude, longitude): Compute IRI altitude profile at a set time and location (latitude/longitude)
#vals = iri2016.IRI('2003-11-21T12', [100,1000,10], -76.77, -11.95)
#altprof = ion.IRI('2015-12-28T12', [60,450,10], 45.5017, -73.5673)
#vals = iri2016.altprofile('2003-11-21T12', -11.95, -76.77)

#"/Users/abrenema/opt/anaconda3/lib/python3.8/site-packages/iri2016/altitude.py"
#def main(time: str, alt_km: T.Sequence[float], glat: float, glon: float):



#Get timeline of data to separate out science collection times 
#tl, gsS, gsE = end_data_loader.load_timeline()
tl = end_data_loader.load_timeline()

slp = end_data_loader.load_slp()
ephem = end_data_loader.load_ephemeris()
#"Official" IRI run on Box is for start of mission
iri = end_data_loader.load_iri()



plot_times = ephem[0]['Time']
plot_alt = ephem[0]['Altitude']
plot_times = plot_times[0::20]
plot_alt = plot_alt[0::20]


#%load_ext nb_black
plt.rcParams['figure.figsize'] = [10, 4]


"""
Load mag data
NOTE: use the IGRF data because the onboard data become inaccurate after the rocket flip.
"""
#mag = EFL('mag')
#magDCx, magDCy, magDCz, tmagDC = mag.load_data()
#Bmag = np.sqrt(magDCx**2 + magDCy**2 + magDCz**2)

#Location of Ny-Alesund
glat = 78 + (55/60)
glon = 11 + (55/60)
Bmag = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])



#Get cyclotron freqs
#Bo = signal.decimate(Bmag, 5)
#tvals = list(signal.decimate(tmagDC,5))
tvals = plot_times
fce = [28*i for i in Bmag]
fcH = [i/1836.15 for i in fce]
fcO = [i/(1836.15*15) for i in fce]
fcN = [i/(1836.15*14) for i in fce]
fcHe = [i/(1836.15*4) for i in fce]



#Load SLP density
tLP = slp['ToF [s]']
dLP = slp['SLP Ni [/m3]']/(100**3)
goodv = np.where(~np.isnan(dLP))
tLP = tLP[goodv[0]]
dLP = dLP[goodv[0]]


#-------------------------------------------------------------------------------
#by-eye lower hybrid values from end_cal_determine_accurate_lower_hybrid_freq.py
#from VLF12-VLF34 phase

#load the by-eye lower hybrid values
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
verticesLOW, verticesHIG = pickle.load(open(flhfile,'rb'))



#Use the upper bound values. These seem to be more accurate and I'm able to identify them 
#over more of the mission.
tlh = verticesHIG[:,0]
flh = verticesHIG[:,1]


#fpeLP = [8980.*np.sqrt(i) for i in dLP]

#Interpret cyclotron freq values to times of flh determination
interp = interp1d(tvals,fce,kind='cubic', bounds_error=False)
fce2 = interp(tlh)
interp = interp1d(tvals,fcH,kind='cubic', bounds_error=False)
fcH2 = interp(tlh)
interp = interp1d(tvals,fcO,kind='cubic', bounds_error=False)
fcO2 = interp(tlh)
interp = interp1d(tvals,Bmag,kind='cubic', bounds_error=False)
Bo2 = interp(tlh)


#Interpolate Langmuir probe values to flh identification times 
interpd = interp1d(tLP,dLP,kind='cubic', bounds_error=False)
dLP2 = interpd(tlh)


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#FOR EACH VALUE OF NE1 FIND THE BEST RATIO OF O+ TO H+ THAT MINIMIZES 
#THE DIFFERENCE B/T NE1 AND DLP2
#nOv = np.arange(0.95,1,0.01)[np.newaxis]
#nHv = np.arange(0.,0.06,0.01)[np.newaxis]
#nf = nOv.T * nHv
nOv = np.arange(0.90,1,0.005)
nHv = np.arange(0.,0.1,0.005)

diffv = np.zeros([len(nHv),len(nOv)])
nHfin = np.zeros(len(flh))
#nOfin = np.zeros(len(flh))

for t in range(len(flh)):
    for h in range(len(nHv)):
        for o in range(len(nOv)):
            netst = dflh.dens_IonMassFractions([flh[t]], [fce2[t]], [nHv[h]], [nOv[o]])
            if not np.isnan(netst[0]):
                diffv[h,o] = np.abs(dLP2[t] - netst[0].value)
            else:
                diffv[h,o] = 1e31

    minv = np.min(diffv)
    goo = np.where(diffv == np.min(diffv))
    if not np.isnan(minv):
        nHfin[t] = nHv[goo[0]][0]

nOfin = 1 - nHfin
tots = nHfin+nOfin


#recalculate density from identified fractions and see how it aligns with SLP density
ne_redo = dflh.dens_IonMassFractions([flh], [fce2], [nHfin], [nOfin])

#ne1 = dflh.dens_IonMassFractions(flh, fce2, nH_ne, nO_ne)
#ne2 = dflh.dens_singleion(flh, Bmag, 'H+')
#ne3 = dflh.dens_singleion(flh, Bmag, 'O+')

#ne1 = [i.value if not np.isnan(i) else 0 for i in ne1]
#ne2 = [i.value if not np.isnan(i) else 0 for i in ne2]
#ne3 = [i.value if not np.isnan(i) else 0 for i in ne3]






fig, axs = plt.subplots(3)
axs[0].title.set_text('end_analysis_determine_composition.py')
axs[0].plot(tlh,nOfin,'.')
axs[1].plot(tlh,nHfin,'.')
axs[2].plot(tlh,ne_redo[0],'.',tLP,dLP)
axs[0].set_ylim(0.9,1)
axs[1].set_ylim(0,0.1)
axs[2].set_ylim(0,300000)
axs[0].set_ylabel('O+ fraction')
axs[1].set_ylabel('H+ fraction')
axs[2].set_ylabel('SLP dens(orange) vs\nreconstructed dens\nfrom by-eye LH determination')
axs[0].set_xlim(100,900)
axs[1].set_xlim(100,900)
axs[2].set_xlim(100,900)
plt.xlabel('time (sec)')
#axs[0].plot(iri['times_upleg'],iri['O_ions']*0.01,'.')
#axs[1].plot(iri['times_upleg'],iri['H_ions']*0.01,'.')
#axs[0].plot(iri['times_downleg'],iri['O_ions']*0.01,'.')
#axs[1].plot(iri['times_downleg'],iri['H_ions']*0.01,'.')



plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)



info = ["From end_determine_composition.py." +
        "Use Langmuir probe derived density values along with id of the flhr line to determine plasma composition." + 
        "ne_flhID_HpOp are the density values using a mixed plasma of H+ and O+ with the fractional ratios relative to ne of nH_ne and nO_ne." + 
        "ne_flhID_Hp and ne_flhID_Op are the same but for 100% H+ and O+ plasmas."]


dict_fin = {'times':tlh,'ne_langmuirprobe':dLP2,
            'ne_flhID_HpOp':ne1,
            'nH_ne':nH_ne,'nO_ne':nO_ne,
            'ne_flhID_Hp':ne2,
            'ne_flhID_Op':ne3,
            'flh':flh,
            'info':info}

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/plasma_composition/plasma_composition.pkl'
pickle.dump(dict_fin, open(path,'wb'))



#----------------------------------------------
#Rough altitude comparison with IRI model O+ % 

import pandas as pd 

#path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/'
#folder = 'ephemeris'
#fn = 'Endurance_GPS_velocity_position_altitude.csv'
#
#header = ['time','xvel','yvel','zvel','lat','long','alt']
#ephem =  pd.read_csv(path + folder + '/' + fn, skiprows=1, names=header)


#IRI O+ % 
#alt_iri = [100.00,150.00,200.00,250.00,300.00,350.00,400.00,450.00,500.00,550.00,600.00,650.00,700.00,750.00,800.00,850.00,900.00]
#Op_iri = [0.0,1.2,24.7,76.1,91.3,96.8,96.2,95.5,94.7,93.8,93.1,92.5,91.7,90.8,89.7,88.4,87.0]
#Hp_iri = [0.0,0.0,0.0,0.0,0.0,0.1,0.2,0.3,0.4,0.7,1.0,1.3,1.7,2.2,2.9,3.6,4.4]

alt_iri = iri['Height(km)']
Op_iri = iri['O_ions']
Hp_iri = iri['H_ions']
Op_iri = [i/100 for i in Op_iri]
Hp_iri = [i/100 for i in Hp_iri]


#Associate rocket times with altitudes so I can compare to IRI
#altv = np.asarray(ephem.alt)
altv = ephem[" Altitde (km)"]
alts = np.empty(len(tlh))

for i in range(len(tlh)):
    good = np.where(ephem['Flight Time'] >= tlh[i])
    good2 = good[0][0]
    alts[i] = altv[good2]


fig, axs = plt.subplots(2)
axs[0].plot(Op_iri, alt_iri, nO_ne, alts, 'x')
axs[1].plot(Hp_iri, alt_iri, nH_ne, alts, 'x')
#axs[1].set_xlim(0,0.1)
plt.show()


