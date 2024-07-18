#Make nice looking spectral plots for paper

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import filter_wave_frequency as filt
import plasma_params_get_flhr_freq as dflh
#import igrf
import pyIGRF

glat = 78 + (55/60)
glon = 11 + (55/60)
#alts = np.arange(100,1000,10)


#mag = igrf.igrf('2022-05-11', glat=glat, glon=glon, alt_km=100)


v34 = EFL('VLF34D')
fs = v34.chnspecs['fs']
wf34, tdat = v34.load_data_gainphase_corrected()

#nps = 2*4096
nps = 4096
fspec, tspec, powerc = signal.spectrogram(wf34, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')


ephem = end_data_loader.load_ephemeris()
alt = ephem[0]['Altitude'] 

iri = end_data_loader.load_iri()
#igrf = end_data_loader.load_igrf()
ig = end_data_loader.load_ig()
slp = end_data_loader.load_slp()
tl = end_data_loader.load_timeline()

slp_times = slp['ToF [s]']
slp_alt = slp['Alt [km]']


#Use this as the plotted time cadence
plot_times = ephem[0]['Time']
plot_alt = ephem[0]['Altitude']
plot_times = plot_times[0::20]
plot_alt = plot_alt[0::20]


BoIGRF = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])



#science collection times
sstart = tl[1]
send = tl[2]
bsstart = tl[3]
bsend = tl[4]

#magnetic field 
#...NOTE: the onboard magnetometer gets wonky after the flip maneuver, likely due to its proximity to sc body. 
#Better to use IGRF data





magv = EFL('mag')
mag = magv.load_data()
Bo = np.sqrt(mag[0]**2 + mag[1]**2 + mag[2]**2)
Bot = mag[3]

#BoIGRF = np.asarray([52617.5,51557.4,50522.1,49511.1,48524.0,47560.2,46619.1,45700.2,44803.1,43927.1,43071.9,42236.8,41421.4,40625.3,39847.8,39088.6,38347.3,37623.2,36916.1])
#altIGRF = np.asarray([100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.])
#Bo1IGRF = igrf['Bo_upleg']
#Bo2IGRF = igrf['Bo_downleg']
#times1IGRF = igrf['times_upleg']
#times2IGRF = igrf['times_downleg']

#plt.plot(times1IGRF,Bo1IGRF)
#plt.plot(times2IGRF,Bo2IGRF)

#Boz = Bo[np.where(Bot > t)]
#Boz = Boz[0]  #nT




#-------------------
#composition
#98% O+ at 270 km (comparison of flh to SLP densities)
#-------------------
fOp2 = 0.98
fHp2 = 0.02

fOp = iri['O_ions']/100
fHp = iri['H_ions']/100

#Interpolate IRI alt to SLP alt
irialt = iri['Height(km)']
fOp = np.interp(plot_alt, irialt, fOp)
fHp = np.interp(plot_alt, irialt, fHp)

plt.plot(plot_times,fOp)


#------------------
#Densities
#------------------

nneut = ig['Dens_neutral']
ni = slp['SLP Ni [/m3]'] / 100**3
ni = np.asarray(ni)
ne = ni #cm-3

ne = np.interp(plot_times, slp_times, ne)
ni = np.interp(plot_times, slp_times, ni)


nOp = fOp * ne 
nHp = fHp * ne


#-------------------------------------------------
#Derived quantities from the NRL formulary online
#-------------------------------------------------

##Gyrofrequencies (Hz)
#fce = 28 * Bo   
##Interpolate high-cadence fce to low-cadence fpe
#fce = np.interp(slp_times, Bot, fce)
#fcH = fce / 1836
#fcO = fcH / 16


#BoIGRF = np.concatenate((Bo1IGRF, Bo2IGRF))
fceIGRF = 28 * BoIGRF
fcHIGRF = fceIGRF / 1836
fcOIGRF = fcHIGRF / 16
#timesIGRF = np.concatenate((times1IGRF, times2IGRF))



#plasma freqs (Hz)
fpe = 8980 * np.sqrt(ne) 
#fpH = 104219
#fpO = 26054


#beta 
#beta_e = fpe / fce


x = np.linspace(0, 2*np.pi, 10)
y = np.sin(x)
xvals = np.linspace(0, 2*np.pi, 50)
yinterp = np.interp(xvals, x, y)


#Lower hybrid freq calculated based on fractional composition (f >> fH)
flh = dflh.flhr_IonMassFractions(ne, fceIGRF, fHp, fOp)
flh2 = dflh.flhr_IonMassFractions(ne, fceIGRF, fHp2, fOp2)

#Lower hybrid freq identified by eye


#Remove data during non science collection times 
for i in range(len(bsstart)):
    goo = np.where((slp_times >= bsstart[i]) & (plot_times < bsend[i]))
    flh[goo] = "nan"
    flh2[goo] = "nan"
    fcHIGRF[goo] = "nan"
    fcOIGRF[goo] = "nan"


"""
#spectral plots 
vr = [-40,-25]
yr = [4000,10000]
ys = 'linear'
xr = [100,900]

fig, axs = plt.subplots()
axs2 = axs.twiny()
axs2.set_ylabel('altitude')

locs = [1000, 2400, 5000, 6000]
axs2.xaxis.set_major_locator(mticker.FixedLocator(locs))
labels = [4,5,6,7]
axs2.set_xticklabels(labels)

ps.plot_spectrogram(tspec,fspec,np.abs(powerc),ax=axs,vr=vr,xr=xr,yr=yr,yscale=ys,xlabel='time (sec)',ylabel="freq (Hz)")#,title="VLF34 power spec (dB)")
"""

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


vertices = np.asarray(vertices)
vertices2 = np.asarray(vertices2)



#spectral plots 
#vr = [-40,-25]
vr = [-86,-60]
yr = [4000,10000]
ys = 'linear'
#xr = [100,900]
xr = [850,900]

#axs.plot(slp_times, flh, color='white')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,xr=xr,yr=yr,yscale=ys,xlabel='time (sec)',ylabel="freq (Hz)")#,title="VLF34 power spec (dB)")
plt.plot(slp_times, flh, color='white')
plt.plot(slp_times, flh2, color='white')
plt.plot(slp_times, fcH, color='white')
plt.plot(slp_times, fcO, color='white')
plt.plot(slp_times, 2*fcH, color='white',linestyle='--')
plt.plot(slp_times, 3*fcH, color='white',linestyle='--')
plt.plot(slp_times, 4*fcH, color='white',linestyle='--')
plt.plot(slp_times, 5*fcH, color='white',linestyle='--')
plt.plot(slp_times, 6*fcH, color='white',linestyle='--')
plt.plot(slp_times, 7*fcH, color='white',linestyle='--')
plt.plot(slp_times, 8*fcH, color='white',linestyle='--')
plt.plot(slp_times, 9*fcH, color='white',linestyle='--')
plt.show()


vr = [-40,-25]
yr = [20,15000]
xr = [100,900]
ys = 'linear'
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yr=yr,xr=xr,yscale=ys,xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
#plt.plot(plot_times, flh,color='white')
plt.plot(vertices[:,0],vertices[:,1],color='magenta',linestyle='dotted')
plt.plot(vertices2[:,0],vertices2[:,1],color='magenta',linestyle='dotted')
plt.plot(plot_times, fcHIGRF, color='white',linestyle='dotted')
plt.plot(plot_times, fcOIGRF, color='white',linestyle='dotted')
#plt.plot(plot_times, 2*fcHIGRF, color='white',linestyle='--')
#plt.plot(plot_times, 3*fcHIGRF, color='white',linestyle='--')
#plt.plot(plot_times, 4*fcHIGRF, color='white',linestyle='--')
#plt.plot(plot_times, 5*fcHIGRF, color='white',linestyle='--')
plt.plot(plot_times, 6*fcHIGRF, color='white',linestyle='dotted')
plt.plot(plot_times, 7*fcHIGRF, color='white',linestyle='dotted')
plt.plot(plot_times, 8*fcHIGRF, color='white',linestyle='dotted')
plt.plot(plot_times, 9*fcHIGRF, color='white',linestyle='dotted')
plt.plot(plot_times, 10*fcHIGRF, color='white',linestyle='dotted')


#plt.plot(slp_times, 2*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 3*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 4*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 5*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 6*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 7*fcH, color='white',linestyle='--')
#plt.plot(slp_times, 8*fcH, color='white',linestyle='--')





#This version shows the fcH+ fundamental
vr = [-40,-37]
#yr = [0,15000]
yr = [500,800]
ys = 'linear'
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yr=yr,yscale=ys,xlabel='time(sec)',ylabel="power spectrum x'\nfreq(Hz)")
plt.plot(slp_times, flh, color='white')
plt.plot(slp_times, fcH, color='white')
plt.plot(slp_times, fcO, color='white')
plt.plot(slp_times, 2*fcH, color='white',linestyle='--')
plt.plot(slp_times, 3*fcH, color='white',linestyle='--')
plt.plot(slp_times, 4*fcH, color='white',linestyle='--')
plt.plot(slp_times, 5*fcH, color='white',linestyle='--')
plt.plot(slp_times, 6*fcH, color='white',linestyle='--')
plt.plot(slp_times, 7*fcH, color='white',linestyle='--')
plt.plot(slp_times, 8*fcH, color='white',linestyle='--')





