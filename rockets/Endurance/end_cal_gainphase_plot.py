"""
Load and plot Endurance gain/phase curves from testing. 

(this is also done within end_fields_loader.py. This routine is for convience in plotting.)


Valid channels are:

'V12D','V34D','V13D','V32D','V24D','V41D'
'VLF12D','VLF34D','VLF13D','VLF32D','VLF24D','VLF41D'
'V1SD','V2SD','V3SD','V4SD'
'V12A','V34A'
'VLF12A'
'V1SA','V2SA','V3SA','V4SA'
'HF12','HF34'
'mag'

"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
import numpy as np 
import matplotlib.pyplot as plt
from math import remainder


chn = 'VLF12D'

path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'


#if chn == 'HF12': fn = 'Endurance_Analog 1_HF12_1000-20000000-100.txt'
#if chn == 'HF12': fn = 'Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt'
if chn == 'HF12': fn = 'Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt'    
if chn == 'HF34': fn = 'Endurance_Analog 1_HF34_1000-20000000-100.txt'
if chn == 'V1SA': fn = 'Endurance_Analog 1_V1SA_10-10000-100.txt'
if chn == 'V1SD': fn = 'Endurance_Analog 1_V1SD_10-10000-100.txt'
if chn == 'V2SA': fn = 'Endurance_Analog 1_V2SA_10-10000-100.txt'
if chn == 'V2SD': fn = 'Endurance_Analog 1_V2SD_10-10000-100.txt'
if chn == 'V3SA': fn = 'Endurance_Analog 1_V3SA_10-10000-100.txt'
if chn == 'V3SD': fn = 'Endurance_Analog 1_V3SD_10-10000-100.txt'
if chn == 'V4SA': fn = 'Endurance_Analog 1_V4SA_10-10000-100.txt'
if chn == 'V4SD': fn = 'Endurance_Analog 1_V4SD_10-10000-100.txt'
if chn == 'V12A': fn = 'Endurance_Analog 1_V12A_10-10000-100.txt'
if chn == 'V12D': fn = 'Endurance_Analog 1_V12D_10-10000-100.txt'
if chn == 'V13D': fn = 'Endurance_Analog 1_V13D_10-10000-100.txt'
if chn == 'V24D': fn = 'Endurance_Analog 1_V24D_10-10000-100.txt'
if chn == 'V32D': fn = 'Endurance_Analog 1_V32D_10-10000-100.txt'
if chn == 'V34A': fn = 'Endurance_Analog 1_V34A_10-10000-100.txt'
if chn == 'V34D': fn = 'Endurance_Analog 1_V34D_10-10000-100.txt'
if chn == 'V41D': fn = 'Endurance_Analog 1_V41D_10-10000-100.txt'
if chn == 'V42D': fn = 'Endurance_Analog 1_V42D_10-10000-100.txt'
if chn == 'VLF12A': fn = 'Endurance_Analog 1_VLF12A_6-100000-100.txt'
if chn == 'VLF12D': fn = 'Endurance_Analog 1_VLF12D_6-30000-100.txt'
if chn == 'VLF13D': fn = 'Endurance_Analog 1_VLF13D_6-30000-100.txt'
if chn == 'VLF24D': fn = 'Endurance_Analog 1_VLF24D_6-30000-100.txt'
if chn == 'VLF32D': fn = 'Endurance_Analog 1_VLF32D_6-30000-100.txt'
if chn == 'VLF34D': fn = 'Endurance_Analog 1_VLF34D_6-30000-100.txt'
if chn == 'VLF41D': fn = 'Endurance_Analog 1_VLF41D_6-30000-100.txt'
if chn == 'VLF42D': fn = 'Endurance_Analog 1_VLF42D_6-30000-100.txt'



with open(path + fn) as f:
    lines = f.readlines()


f = lines[0].split()  #freq in Hz
p = lines[1].split()  #phase in deg
g = lines[2].split()  #gain in dB
f = [float(i) for i in f]
p = [float(i) for i in p]
g = [float(i) for i in g]

#--change to radians
prad = [np.deg2rad(i) for i in p]


#--"Fix" the phase data so that it can only vary b/t -pi and pi. Some of the curves 
#--extend a bit beyond this for some reason. To fix, unwrap and then rewrap angles.
prad = np.unwrap(prad)
prad = np.asarray([remainder(prad[i], 2*np.pi) for i in range(len(prad))])


#--change gain from dB to linear scale for calculation of transfer function
#--From Steve Martin email on Nov 7, 2022: 
#--Gain=10^(0.05 * (opchan+gainoffset))
#--NOTE: This gives the correct max gain values for all the channels EXCEPT for HF.
offset = 0.
Hmag = [10**(0.05*i + offset) for i in g]



#--Remove bad data points that can occur at very low frequencies. This happens b/c 
#--the signal/noise ratio can become high leading to artificial values at very low freqs
#--during the gain/phase tests. 

fig, axs = plt.subplots(3)
axs[0].set_title(chn)
axs[0].plot(f,g)
axs[1].plot(f,Hmag)
axs[2].plot(f,prad)
axs[0].set_ylabel('gain (dB)')
axs[1].set_ylabel('gain (linear)')
axs[2].set_ylabel('phase (rad)')
axs[2].set_xlabel('frequency')
for i in range(3):
    axs[i].set_xlim(0,8000)
print('h')




"""
#compare phase of VLF12 and VLF34
prad12 = prad
f12 = f
prad34 = prad 
f34 = f

diff = prad12 - prad
diff = diff * (180/3.14)
plt.plot(f,diff)
plt.xlim(0,10000)
plt.ylim(-5,5)

plt.plot(f,prad,f,prad12)
plt.xlim(0,10000)
"""