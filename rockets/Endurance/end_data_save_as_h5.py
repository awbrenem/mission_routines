"""
Save Endurance waveform data into an .h5 file for use in IDL
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
from end_fields_loader import Endurance_Fields_Loader as EFL
#import end_data_loader
import h5py
import numpy as np



#Creating an .h5 file - this works fine on Python end but I can't load it with IDL
#https://www.christopherlovell.co.uk/blog/2016/04/27/h5py-intro.html
v12 = EFL('VLF12D')
v34 = EFL('VLF34D')

fs = v12.chnspecs['fs']
wf12, tvals = v12.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()


hf = h5py.File('datat.h5', 'w')
hf.create_dataset('times', data=tvals)
hf.create_dataset('wf12', data=wf12)
hf.create_dataset('wf34', data=wf34)
hf.close()







