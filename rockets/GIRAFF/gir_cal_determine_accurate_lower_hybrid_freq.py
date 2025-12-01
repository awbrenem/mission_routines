"""
Use phase and coherence spectra to determine accurate lower hybrid frequency 

NOTE: Much harder to identify lower hybrid freq on GIRAFF than on Endurance. 
36.381 isn't too bad, but on 36.380 the LH freq is only lit up when there's lots of DC irregularities, it seems.

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


pld = '380'

#Load data for two channels of interest
#c1s = 'VLF12D'
#c2s = 'VLF34D'
c1s = 'V12D'
c2s = 'V34D'


c1 = GFL(pld,c1s)
fs = c1.chnspecs['fs']
#wf1, tdat = c1.load_data_gainphase_corrected()
wf1, tdat = c1.load_data()
c2 = GFL(pld,c2s)
#wf2, tdat = c2.load_data_gainphase_corrected()
wf2, tdat = c2.load_data()




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

yr = [0,50000]
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






yr = [1000,12000]
yscale='linear'
xr = [100,500]
vr = [-80,-60]

fig,axs = plt.subplots(3)
ps.plot_spectrogram(tcenter_fin,freqs_fin,np.abs(spec_fin1),vr=vr,yscale='linear',yr=yr,xr=xr,ylabel="power spectrum VLF12\nfreq(Hz)\ndB of (mV/m)^2/Hz",ax=axs[0])
ps.plot_spectrogram(tcenter_fin,freqs_fin,phase_finz,vr=[-60,60], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='phase',ax=axs[1])
ps.plot_spectrogram(tcenter_fin,freqs_fin,coh_fin**2,vr=[0.7,1], zscale='linear',xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',ylabel='coh^2',ax=axs[2])
klicker = clicker(axs[1],["event"], markers=["x"])
vertices1 = (klicker.get_positions())['event']


#36.380

vertices = ([[ 113.66086668, 5197.21011334],
       [ 121.51970235, 5451.2890771 ],
       [ 129.1329494 , 5705.36804085],
       [ 134.2903103 , 6721.68389588],
       [ 141.41238012, 6975.76285963],
       [ 147.55209548, 7102.80234151],
       [ 155.16534253, 7060.45584755],
       [ 160.81388067, 6806.3768838 ],
       [ 167.44477326, 6552.29792004],
       [ 175.30360892, 6764.03038984],
       [ 183.4080332 , 6679.33740192],
       [ 189.79333718, 6594.644414  ],
       [ 197.16099562, 6382.9119442 ],
       [ 203.79188821, 6425.25843816],
       [ 211.65072387, 6425.25843816],
       [ 219.50955954, 6636.99090796],
       [ 226.63162936, 6933.41636567],
       [ 233.50811056, 7060.45584755],
       [ 239.40223731, 7187.49532943],
       [ 246.76989575, 7187.49532943],
       [ 251.43607942, 7018.10935359],
       [ 257.57579479, 7272.18831735],
       [ 263.71551015, 7441.57429319],
       [ 272.80228889, 7441.57429319],
       [ 280.41553594, 7356.88130527],
       [ 287.04642853, 7229.84182339],
       [ 295.39644142, 6975.76285963],
       [ 304.48322016, 6975.76285963],
       [ 311.60528998, 7102.80234151],
       [ 319.46412564, 7229.84182339],
       [ 328.55090438, 7060.45584755],
       [ 337.14650589, 7145.14883547],
       [ 347.70681631, 7229.84182339],
       [ 356.05682921, 7314.53481131],
       [ 362.19654457, 7272.18831735],
       [ 370.30096885, 7060.45584755],
       [ 378.89657035, 7018.10935359],
       [ 389.21129216, 7018.10935359],
       [ 398.78924813, 6679.33740192],
       [ 407.38484964, 6721.68389588],
       [ 414.99809669, 6467.60493212],
       [ 426.29517296, 6171.17947441],
       [ 433.66283139, 5705.36804085],
       [ 440.78490121, 4773.74517375]])

#36.381
vertices = ([[ 109.66355237, 5987.06199461],
       [ 116.4237056 , 5968.5309973 ],
       [ 123.53053336, 6190.90296496],
       [ 128.3839767 , 6005.59299191],
       [ 134.62411815, 5912.93800539],
       [ 138.26420066, 7006.26684636],
       [ 140.3442478 , 6431.80592992],
       [ 145.37102841, 6505.92991914],
       [ 149.53112271, 6654.17789757],
       [ 160.10469571, 6580.05390836],
       [ 167.0381862 , 6765.3638814 ],
       [ 174.66502574, 6542.99191375],
       [ 180.73182992, 6617.11590296],
       [ 185.23859874, 6802.42587601],
       [ 190.9587284 , 6524.46091644],
       [ 196.85219532, 6431.80592992],
       [ 203.43901128, 6524.46091644],
       [ 209.85248999, 6709.77088949],
       [ 216.61264322, 6932.14285714],
       [ 223.89280824, 7098.92183288],
       [ 228.74625158, 6802.42587601],
       [ 236.0264166 , 6450.33692722],
       [ 241.74654626, 6691.23989218],
       [ 247.81335044, 6783.89487871],
       [ 254.05349188, 6820.95687332],
       [ 262.02700595, 6672.70889488],
       [ 269.30717096, 7228.63881402],
       [ 272.77391621, 6876.54986523],
       [ 280.22741849, 6895.08086253],
       [ 288.02759529, 6709.77088949],
       [ 297.04113293, 6728.30188679],
       [ 304.66797247, 6728.30188679],
       [ 311.60146297, 6635.64690027],
       [ 319.22830251, 6709.77088949],
       [ 325.98845574, 6932.14285714],
       [ 332.22859718, 6932.14285714],
       [ 342.10882113, 7024.79784367],
       [ 349.21564889, 6876.54986523],
       [ 354.93577854, 6876.54986523],
       [ 366.02936333, 6895.08086253],
       [ 381.97639146, 6876.54986523],
       [ 391.16326637, 6895.08086253],
       [ 408.15031807, 6987.73584906],
       [ 419.41724012, 6839.48787062]])


vertices = np.asarray(vertices)

#Save vertices as a pickle file
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/lower_hybrid_id/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
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

