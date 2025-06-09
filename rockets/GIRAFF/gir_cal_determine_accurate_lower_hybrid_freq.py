"""
Use phase and coherence spectra to determine accurate lower hybrid frequency 
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import numpy as np 
import plot_spectrogram as ps
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker
import pickle
import correlation_analysis as ca


pld = '381'

#Load data for two channels of interest
c1s = 'VLF12D'
c2s = 'VLF34D'
#c1s = 'V12D'
#c2s = 'V34D'


c1 = GFL(pld,c1s)
fs = c1.chnspecs['fs']
wf1, tdat = c1.load_data_gainphase_corrected()
c2 = GFL(pld,c2s)
wf2, tdat = c2.load_data_gainphase_corrected()




nfft = 16384
freqs_fin, tcenter_fin, csd_fin, coh_fin, phase_fin, fs_fin, spec_fin1, spec_fin2 = ca.csd_spectrum_piecewise(tdat[tdat > 100], wf2[tdat > 100], wf1[tdat > 100], nfft=nfft, noverlap=8, fs_thres=0.2)
phase_fin = phase_fin*(180/3.14)


phase_finz = np.copy(phase_fin)
coh_finz = np.copy(coh_fin)

#Filter phase and coherence values by spectral power
pow = np.abs(spec_fin1)
minpow = 5e-8
for i in range(len(tcenter_fin)):
      goo = np.where(pow[:,i] < minpow)[0]
      phase_finz[goo,i] = np.nan
      #coh_finz[goo,i] = np.nan


#cohmin = 0.01  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result
#cohmin = 0.3 

yr = [6000,8000]
yscale='linear'
xr = [100,550]
vr = [-80,-60]

fig,axs = plt.subplots(3)
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tcenter_fin,freqs_fin,np.abs(spec_fin1),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_finz,vr=[-60,60], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='phase',ax=axs[1])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin**2,vr=[0.7,1], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='coherence^2',ax=axs[2])
axs[0].get_xaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)






yr = [4000,8000]
yscale='linear'
xr = [400,500]
vr = [-80,-60]

fig,axs = plt.subplots(2)
ps.plot_spectrogram(tcenter_fin,freqs_fin,np.abs(spec_fin1),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_finz,vr=[-60,60], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='phase',ax=axs[1])
klicker = clicker(axs[0],["event"], markers=["x"])
vertices1 = (klicker.get_positions())['event']



vertices = [[ 109.09934342, 6071.26805778],
[ 112.09677419, 5782.3434992 ],
[ 116.33123989, 6007.06260032],
[ 119.18593586, 5840.12841091],
[ 121.99305357, 5987.80096308],
[ 124.70501475, 6141.894061  ],
[ 127.60728899, 6064.84751204],
[ 131.36597202, 6141.894061  ],
[ 134.17308973, 6109.79133226],
[ 136.17137692, 6103.37078652],
[ 137.97935103, 6103.37078652],
[ 142.26139499, 6340.93097913],
[ 144.73546484, 6597.75280899],
[ 147.11437815, 6899.51845907],
[ 150.15938719, 6777.52808989],
[ 153.01408317, 6713.32263242],
[ 155.53573128, 6706.90208668],
[ 159.00894471, 6777.52808989],
[ 161.95879722, 6848.1540931 ],
[ 166.24084118, 6674.79935795],
[ 169.66647635, 6796.78972713],
[ 173.47273765, 6745.42536116],
[ 177.18384242, 6764.68699839],
[ 180.41916453, 6873.83627608],
[ 184.89152155, 6950.88282504],
[ 187.69863926, 6681.21990369],
[ 190.26786564, 6674.79935795],
[ 192.40888762, 6732.58426966],
[ 195.45389666, 6661.95826645],
[ 198.45132743, 6816.05136437],
[ 200.32764827, 6842.72634791],
[ 203.60413099, 6790.40837088],
[ 205.7884528,  6900.85743351],
[ 209.32705414, 6807.84769656],
[ 212.69090973, 6860.16567359],
[ 216.01107888, 6958.98851911],
[ 219.15650229, 7011.30649615],
[ 221.73400203, 6970.61473623],
[ 226.01527278, 6970.61473623],
[ 228.68014539, 6819.47391368],
[ 231.43239087, 6819.47391368],
[ 235.01467864, 6691.58552536],
[ 237.81061056, 6802.034588  ],
[ 241.21815259, 6906.67054207],
[ 244.66938105, 7069.43758175],
[ 247.07213504, 6970.61473623],
[ 254.06196484, 6964.80162767],
[ 257.73162548, 6953.17541055],
[ 260.44018453, 6825.28702224],
[ 263.49823507, 6871.79189071],
[ 267.12420928, 6906.67054207],
[ 267.997938,   7354.27990118],
[ 270.61912417, 7319.40124982],
[ 273.89560689, 6831.10013079],
[ 277.21577605, 6906.67054207],
[ 281.89022472, 7017.11960471],
[ 286.60835984, 6958.98851911],
[ 290.62751197, 6813.66080512],
[ 294.82140985, 6848.53945647],
[ 299.40848565, 6796.22147944],
[ 300.74506699, 6573.44701583],
[ 303.32500044, 6824.01252828],
[ 307.10593741, 6803.13206891],
[ 310.84239275, 6663.92900644],
[ 313.68921588, 6643.04854707],
[ 317.24774478, 6698.72977205],
[ 321.38453463, 6935.37497825],
[ 324.54272904, 6914.49451888],
[ 328.50159244, 6796.17191578],
[ 332.10460296, 7060.65773447],
[ 335.26279736, 7067.61788759],
[ 339.57751366, 7025.85696885],
[ 343.40293223, 7213.78110318],
[ 346.4276818 , 7074.57804072],
[ 350.43102681, 7039.7772751 ],
[ 353.63370283, 6984.09605011],
[ 357.14775012, 6942.33513137],
[ 361.06213191, 6984.09605011],
[ 365.51029304, 6956.25543762],
[ 370.13638062, 7088.49834696],
[ 378.58788677, 7004.97650948],
[ 382.19089728, 6949.2952845 ],
[ 387.12835614, 6914.49451888],
[ 390.10862409, 6977.13589699],
[ 394.91263812, 7276.42248129],
[ 400.82725107, 7005.55914674],
[ 405.2234996 , 6889.20491273],
[ 409.14703324, 6785.77892696],
[ 413.07056688, 6869.8125404 ],
[ 417.27772946, 7096.05688429],
[ 421.34307756, 6966.77440207],
[ 425.03025375, 6869.8125404 ]]



vertices = np.asarray(vertices)

#Save vertices as a pickle file
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_381_lower_hybrid_freqs_byeye.pkl'
#vertices = pickle.load(open(flhfile,'rb'))
pickle.dump([vertices], open(flhfile,'wb'))




#----------------------------------------------------
#Test plot the vertices
#----------------------------------------------------

yr = [4000,8000]
yscale='linear'
xr = [100,550]
vr = [-80,-60]

fig,axs = plt.subplots(3)
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)
ps.plot_spectrogram(tcenter_fin,freqs_fin,np.abs(spec_fin1),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_finz,vr=[-60,60], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='phase',ax=axs[1])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin**2,vr=[0.7,1], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='coherence^2',ax=axs[2])
axs[0].get_xaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[0].plot(vertices[:,0],vertices[:,1],'*',color='magenta')
axs[1].plot(vertices[:,0],vertices[:,1],'*',color='magenta')
axs[2].plot(vertices[:,0],vertices[:,1],'*',color='magenta')















plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)

