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
vertices = [[ 113.97920064, 7542.68022999],
       [ 114.65698007, 7688.43231072],
       [ 117.36809777, 8005.06614129],
       [ 118.04587719, 8221.18129548],
       [ 123.4681126 , 8400.3223835 ],
       [ 128.27659056, 8320.57029883],
       [ 136.07702933, 8254.57144857],
       [ 142.57739498, 8134.59904222],
       [ 148.42772406, 7899.65307979],
       [ 153.62801658, 7779.68067345],
       [ 162.07849192, 7659.7082671 ],
       [ 169.87893069, 7539.73586075],
       [ 179.62947916, 7454.75540626],
       [ 186.77988137, 7369.77495176],
       [ 191.98017388, 7289.79334753],
       [ 203.03079548, 7224.80829409],
       [ 214.08141708, 7129.83013907],
       [ 223.83196555, 7044.84968457],
       [ 235.53262371, 6959.86923008],
       [ 245.28317218, 6889.88532637],
       [ 260.88404973, 6809.90372214],
       [ 277.78500041, 6759.9152195 ],
       [ 291.43576826, 6744.91866871],
       [ 303.13642642, 6704.92786659],
       [ 319.38734054, 6669.93591474],
       [ 340.1885106 , 6614.94856183],
       [ 353.83927846, 6589.95431051],
       [ 373.99041196, 6584.95546024],
       [ 392.84147233, 6589.95431051],
       [ 412.99260583, 6559.96120892],
       [ 429.24351995, 6504.97385601],
       [ 442.8942878 , 6484.97845496],
       [ 464.995531  , 6469.98190416],
       [ 479.94637198, 6479.97960469],
       [ 496.1972861 , 6524.96925707],
       [ 516.3484196 , 6534.9669576 ],
       [ 533.24937028, 6514.97155654],
       [ 552.75046721, 6509.97270628],
       [ 570.95149102, 6519.97040681],
       [ 585.25229544, 6534.9669576 ],
       [ 598.9030633 , 6639.94281315],
       [ 610.94006746, 6588.7105021 ],
       [ 623.39358298, 6624.50148277],
       [ 639.30640837, 6660.29246344],
       [ 651.7599239 , 6757.43941098],
       [ 660.75412955, 6823.90837509],
       [ 677.35881692, 6915.94232539],
       [ 697.42281415, 6926.16831987],
       [ 710.56819165, 6967.07229778],
       [ 727.17287901, 7053.99325084]]


vertices2 = [[ 118.55265075, 7964.38255622],
       [ 126.59693574, 8254.90835695],
       [ 134.64122073, 8530.14332606],
       [ 145.90321972, 8514.85249444],
       [ 156.36079021, 8438.39833636],
       [ 162.7962182 , 8316.07168342],
       [ 170.8405032 , 8102.00004078],
       [ 179.68921669, 8010.25505107],
       [ 186.92907318, 7903.21922975],
       [ 199.79992917, 7796.18340843],
       [ 208.64864266, 7643.27509225],
       [ 223.12835565, 7551.53010255],
       [ 233.58592614, 7475.07594446],
       [ 243.23906813, 7199.84097535],
       [ 255.30549562, 7459.78511285],
       [ 266.56749461, 7398.62178638],
       [ 270.58963711, 7169.25931212],
       [ 281.8516361 , 7199.84097535],
       [ 293.11363508, 7062.22349079],
       [ 302.76677708, 7352.74929152],
       [ 311.61549057, 7291.58596506],
       [ 320.46420406, 7092.80515403],
       [ 330.92177455, 6955.18766947],
       [ 338.16163104, 7016.35099594],
       [ 348.61920153, 6939.89683786],
       [ 358.27234352, 6955.18766947],
       [ 367.12105702, 6832.86101653],
       [ 377.57862751, 6909.31517462],
       [ 388.036198  , 6894.024343  ],
       [ 399.29819699, 7138.67764888],
       [ 404.92919648, 6848.15184815],
       [ 409.75576748, 7092.80515403],
       [ 425.03990896, 6955.18766947],
       [ 434.69305095, 6848.15184815],
       [ 453.99933494, 6832.86101653],
       [ 468.47904792, 6817.57018492],
       [ 474.91447592, 6786.98852168],
       [ 499.0473309 , 6802.2793533 ],
       [ 511.91818688, 6863.44267977],
       [ 520.76690038, 6924.60600624],
       [ 532.02889936, 6955.18766947],
       [ 542.48646986, 6955.18766947],
       [ 565.81489633, 6924.60600624],
       [ 577.88132382, 7092.80515403],
       [ 588.33889431, 6939.89683786],
       [ 626.95146228, 7184.55014373],
       [ 641.43117526, 7230.42263859],
       [ 651.88874575, 7337.45845991],
       [ 664.75960174, 7429.20344961],
       [ 680.04374323, 7429.20344961],
       [ 693.71902772, 7520.94843932],
       [ 706.5898837 , 7444.49428123],
       [ 725.09173919, 7475.07594446],
       [ 736.35373818, 7475.07594446],
       [ 746.81130867, 7612.69342902],
       [ 763.70430715, 7673.85675549],
       [ 775.77073464, 7811.47424005],
       [ 786.22830513, 7903.21922975],
       [ 797.49030412, 8086.70920916],
       [ 807.14344611, 8224.32669371],
       [ 819.2098736 , 8407.81667312],
       [ 832.88515809, 8499.56166283],
       [ 841.73387158, 8468.97999959],
       [ 852.99587057, 8331.36251504],
       [ 857.01801307, 8147.87253563],
       [ 859.43129856, 7887.92839813],
       [ 863.45344106, 7536.23927093],
       [ 869.08444056, 7306.87679667],
       [ 871.49772605, 7016.35099594],
       [ 877.12872555, 6710.5343636 ],
       [ 877.93315405, 6312.97274154],
       [ 881.95529654, 5961.28361434]]


#Use the upper bound values. These seem to be more accurate and I'm able to identify them 
#over more of the mission.
vertices2 = np.asarray(vertices2)
tlh = vertices2[:,0]
flh = vertices2[:,1]


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


