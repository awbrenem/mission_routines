"""

Calibrate the GIRAFF Langmuir probe data from current to density. 

Note that GIRAFF had cylindrical probes (identical to 45.007 Dissipation)
However, I think the below numbers should still scale as the collecting area

NOTE: Three part comparison. 
1) determine density values from HF freq, radar, or LH freq (assuming some ion composition)
2) compare these to calibrated LP values. 

NOTE 2: the density from ion composition is VERY sensitive to the exact O+ to H+ ratio, so this method 
is somewhat less than ideal. 


---------------------
GIRAFF cylindrical probes (tip only):
---------------------

radius = 0.25 cm 
length = 5.7 cm 
area = 2*pi*r*l = 8.95 cm^2


36.380: Biased at +5 V for electrons
36.381: Biased at -2.2 V for ions

---------------------------
REFERENCE PAYLOAD
---------------------------
From Rob's notes (reference from a different mission): 
LP1 (electrons): +5V; 1.25" diameter sphere = 3.175 cm. Area = 31.67 cm2
LP2 (ions):      -5V; 1" diameter sphere = 2.54 cm. Area = 20.26 cm2

Density to current relationships:
LP1: 1.3e5/cm3 = 32 microAmps
LP2: 1.3e6/cm3 = 3.9 microAmps

---------------------------

Using these references, adjust linearly for the different area and voltages on GIRAFF. 
(e.g. doubling area/voltage means doubling the amount of current corresponding to the same density)



Ratio of areas: 

area1_ratio = 31.67 cm2 / 8.95 cm2 = 3.54
area2_ratio = 20.26 cm2 / 8.95 cm2 = 2.26

Ratio of voltages:

volt1_ratio = +5 / +5 = 1 
volt2_ratio = -5 / -2.2 = 2.27 

---------------------------

---------------------------
GIRAFF CURRENT TO DENSITY SCALED RESPONSE
---------------------------

For Electrons:
REFERENCE: +5V: 1.3e5/cm3 = 32 microAmps
GIRAFF: 1.3e5/cm3 = 32 / (area1_ratio * volt1_ratio) = 32 / (3.54*1) = 9.04 microamps

For Ions:
REFERENCE: -5V: 1.3e6/cm3 = 3.9 microAmps
GIRAFF: 1.3e6/cm3 = 3.9 / (area2_ratio * volt2_ratio) = 3.9 / (2.26*2.27) = 0.76 microamps



------------------------------------------------------------------------------
ACTUAL DENSITY VALUES FOR CALIBRATION
------------------------------------------------------------------------------

******
PFISR values (MADRIGAL database) for beam 12 [these are from spectral data and very rough]
36.380 (08:35 UT, Feb 9): density near 300 km ranges from about 10^11.3 to 10^11.6 / m^3
36.381 (07:07 UT, Feb 2): density near 300 km ranges from about 10^11.2 to 10^11.6 / m^3

******
Upper hybrid values from E-FIELDS HF
36.380
36/381 UH near 300 km is about 2.1 MHz 

******
Lower hybrid values from E-FIELDS VLF
36.380 
36.381 LH near 300 km is about 7 kHz. Assuming 100% O+,

------------------------------------------
Useful references:
Steigies, C. (2005): Surface Property Effects on Langmuir Probes Launched on Sounding Rockets

Brace, Larry (1998) in Measurement techniques in Space Plasmas. 
Larry Brace paper on LP theory
"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GLF
import gir_load_data
from scipy import signal
import numpy as np 
import pickle
import matplotlib.pyplot as plt
import fft_spectrum_piecewise as fftspec
import plot_spectrogram as ps
import plasma_params_get_density_from_flhr_freq as flhr


pld = '381'
#pld = '380'



#Test values for 36.381 at 300 sec (~300 km alt from PFISR)
if pld == '380':
    densPFISR = np.asarray([10**11.2, 10**11.6])  #m-3
    densPFISR = densPFISR / (100**3)
if pld == '381':
    densPFISR = np.asarray([10**11.2, 10**11.6])  #m-3
    densPFISR = densPFISR / (100**3)




times, lp = gir_load_data.load_slp(pld)
plt.plot(times,lp)



#flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_381_lower_hybrid_freqs_byeye.pkl'
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
vertices = pickle.load(open(flhfile,'rb'))[0]



#Load E-fields data
mag = GLF(pld,'mag')
magx,magy,magz,tmag = mag.load_data()
Bo = np.sqrt(magx**2 + magy**2 + magz**2)


v12 = GLF(pld,'VLF12D')
wf12, tdat = v12.load_data()



#---------------------------------
#Density from UH freq (**NOTE: very sensitive to the exact UH freq)
#---------------------------------

#36.381
if pld == '381':
    fuh = 2.2e6  #MHz at 300 km alt (~300 sec)
if pld == '380':
    fuh = 0


Bo300 = np.interp(300,tmag,Bo)
fce300 = 28*Bo300

fpe = np.sqrt(fuh**2 - fce300**2)
ne_uh = (fpe/8980)**2



#---------------------------------
#Overall spectral plot (Fig 1)
#---------------------------------

nfft=16384
fspec, tspec, powerc, fs = fftspec.fft_spectrum_piecewise(tdat, wf12, fs_thres=0.1, nfft=nfft, noverlap=8)

xr = [100,550]
yr = [0,50000]
vr = [-90,-50]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")
plt.plot(vertices[:,0],vertices[:,1],'*',color='magenta')



#---------------------------------------------------------------
#Interpolate the Bo values to the times of the identified lower hybrid values (vertices)
#---------------------------------------------------------------

Boz = np.zeros(len(vertices[:,0]))
for i in range(len(vertices)):
    Boz[i] = np.interp(vertices[i,0], tmag, Bo)


#----------------------------------------------------
#Load and calibrate Langmuir probe values. 
#----------------------------------------------------

lpt, lpv = gir_load_data.load_slp(pld)


if pld == '380':
    lpv_cal = lpv * (1.3e5 / 9.04) * 1e6
if pld == '381':
    lpv_cal = lpv * (1.3e6 / 0.76) * 1e6



#----------------------------------------------------
#Turn lower hybrid freqs into density values 
#----------------------------------------------------

neO = flhr.dens_singleion(vertices[:,1], Boz, 'O+')
#neH = flhr.dens_singleion(vertices[:,1], Boz, 'H+')

nH_ne = [0.01]*len(vertices[:,0])
##nN_ne = [0.5]*len(vertices[:,0])
nO_ne = [0.99]*len(vertices[:,0])
ne2 = flhr.dens_IonMassFractions(vertices[:,1], 28*Boz, nH_ne=nH_ne, nO_ne=nO_ne)

#nH_ne = [0.01]*len(vertices[:,0])
#nN_ne = [0.5]*len(vertices[:,0])
#nO_ne = [0.99]*len(vertices[:,0])
#ne2 = flhr.dens_IonMassFractions([7000,7000], 28*[Boz[0],Boz[0]], nH_ne=[0.01,0.01], nO_ne=[0.99,0.99])



#----------------------------------------------------
#compare to "official" density values
#----------------------------------------------------


plt.plot(lpt,np.abs(lpv_cal))
if pld == '380':
    plt.title('36.380')
    plt.ylabel('Langmuir Probe density\n(cm-3)')
if pld == '381':
    plt.title('36.381')
    plt.ylabel('Langmuir Probe density\n(cm-3)')
plt.xlabel('time since launch (sec)')
plt.plot(vertices[:,0],neO,'.',color='red')  #100% O+
plt.plot([300,300],densPFISR)
plt.plot(vertices[:,0],ne2, '.',color='blue')
plt.plot(305,ne_uh,'*',color='green')


plt.yscale('log')
plt.ylim(1e3,1e6)
#plt.ylim(1e2,1e6)
plt.xlim(100,500)
#plt.xlim(50,520)

#Value at 300 km (~300 sec)









#quick look at harmonics at 450 sec 

Bo450 = np.interp(450,tmag,Bo)
fcH450 = 28*Bo450/1836


xr = [400,500]
yr = [1000,6000]
vr = [-90,-50]
ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz")
plt.plot(vertices[:,0],vertices[:,1],'*',color='magenta')
plt.plot([450,450,450,450,450],[fcH450,2*fcH450,3*fcH450,4*fcH450,5*fcH450],'*')
