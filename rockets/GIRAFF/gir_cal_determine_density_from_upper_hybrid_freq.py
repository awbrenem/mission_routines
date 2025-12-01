"""
Determine the density from the identification of the upper hybrid freq
See Kurth+15 RBSP paper.
He uses the center freq of the HF power band as fUH. 

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
from mpl_point_clicker import clicker






pld = '380'
chn = 'HF12'


c1 = GFL(pld,chn)
fs = c1.chnspecs['fs']
#wf1, tdat = c1.load_data_gainphase_corrected()
powerHF, tHF, fHF = c1.load_data()


plt.figure(figsize=(10,6))
axs = plt.subplot(111)
ps.plot_spectrogram(tHF,fHF,powerHF,vr=[-30,-2],yr=[1.5e6,4.0e6], yscale='linear',ax=axs)
klicker = clicker(axs,["event"], markers=["x"])
vertices1 = (klicker.get_positions())['event']


#36.380 (upper hybrid middle values - Kurth method)
vertices1 = ([[9.33831795e+01, 2.14606742e+06],
       [9.81698688e+01, 2.29980613e+06],
       [1.06717528e+02, 2.46396781e+06],
       [1.18684252e+02, 2.63334098e+06],
       [1.26890005e+02, 2.80271414e+06],
       [1.34753851e+02, 2.94342415e+06],
       [1.43643417e+02, 3.12843176e+06],
       [1.51507264e+02, 3.30562215e+06],
       [1.62106362e+02, 3.46196660e+06],
       [1.69970208e+02, 3.60007088e+06],
       [1.76466429e+02, 3.73296367e+06],
       [1.90484591e+02, 3.65479144e+06],
       [2.03818940e+02, 3.65479144e+06],
       [2.19204727e+02, 3.54274457e+06],
       [2.33222888e+02, 3.48541827e+06],
       [2.41770548e+02, 3.35773697e+06],
       [2.50660113e+02, 3.20660399e+06],
       [2.59207773e+02, 3.10237435e+06],
       [2.74593560e+02, 3.01899064e+06],
       [2.86218377e+02, 2.95124137e+06],
       [2.97843194e+02, 2.84440600e+06],
       [2.99552726e+02, 2.71672469e+06],
       [3.10151823e+02, 2.60728357e+06],
       [3.33401457e+02, 2.52389986e+06],
       [3.53232027e+02, 2.48220800e+06],
       [3.69643533e+02, 2.43791040e+06],
       [3.87422665e+02, 2.34149799e+06],
       [4.01098920e+02, 2.26332576e+06],
       [4.13065643e+02, 2.18515353e+06],
       [4.19219958e+02, 2.06789519e+06],
       [4.29135243e+02, 1.91415647e+06]])


#36.381 (middle values - Kurth method for fUH)
vertices1 = ([[9.56076957e+01, 2.13072183e+06],
       [1.00315388e+02, 2.10576767e+06],
       [1.02253850e+02, 2.05368942e+06],
       [1.08900004e+02, 2.04175482e+06],
       [1.12500004e+02, 2.10576767e+06],
       [1.17761543e+02, 2.14591132e+06],
       [1.25515389e+02, 2.15567599e+06],
       [1.32992312e+02, 2.21534898e+06],
       [1.42684620e+02, 2.20124445e+06],
       [1.50438467e+02, 2.17737526e+06],
       [1.57084621e+02, 2.15025117e+06],
       [1.66776929e+02, 2.11227745e+06],
       [1.73146160e+02, 2.07972854e+06],
       [1.78407699e+02, 2.04500971e+06],
       [1.88100007e+02, 2.04609467e+06],
       [1.98346161e+02, 2.05585935e+06],
       [2.03884623e+02, 2.07647365e+06],
       [2.11915392e+02, 2.06453905e+06],
       [2.20776931e+02, 2.05260445e+06],
       [2.30469239e+02, 2.02439540e+06],
       [2.35176931e+02, 1.93542840e+06],
       [2.40992316e+02, 1.99835628e+06],
       [2.46807701e+02, 2.05151949e+06],
       [2.50407701e+02, 2.09600299e+06],
       [2.53730778e+02, 2.13723161e+06],
       [2.57330778e+02, 2.16110080e+06],
       [2.64807702e+02, 2.12204212e+06],
       [2.72284625e+02, 2.07430372e+06],
       [2.80592318e+02, 2.04934956e+06],
       [2.91115395e+02, 2.04609467e+06],
       [2.99423088e+02, 2.05151949e+06],
       [3.09392319e+02, 2.04826460e+06],
       [3.17976934e+02, 2.04175482e+06],
       [3.28223089e+02, 2.03958489e+06],
       [3.35700012e+02, 2.04717964e+06],
       [3.44284628e+02, 2.06236913e+06],
       [3.50653859e+02, 2.08840825e+06],
       [3.58407705e+02, 2.11444737e+06],
       [3.63946167e+02, 2.14591132e+06],
       [3.71976936e+02, 2.17303540e+06],
       [3.82776937e+02, 2.21534898e+06]])


vertices = np.asarray(vertices1)

#Save vertices as a pickle file
fuhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/upper_hybrid_id/GIRAFF_'+pld+'_HF_middle_vals_byeye.pkl'
#vertices = pickle.load(open(flhfile,'rb'))
pickle.dump([vertices], open(fuhfile,'wb'))



#36.381 (lower boundary values)
vertices1 = ([[1.09772710e+02, 1.83999334e+06],
       [1.14823879e+02, 1.89092280e+06],
       [1.16267070e+02, 1.93107873e+06],
       [1.26730204e+02, 1.94283168e+06],
       [1.27812597e+02, 2.02412294e+06],
       [1.29616586e+02, 2.07897005e+06],
       [1.34667754e+02, 2.11227008e+06],
       [1.38275731e+02, 2.06917592e+06],
       [1.39718922e+02, 2.02216411e+06],
       [1.40440518e+02, 1.96144052e+06],
       [1.48738866e+02, 1.98396701e+06],
       [1.55594023e+02, 1.95556404e+06],
       [1.59923595e+02, 2.01334940e+06],
       [1.64974764e+02, 2.09072300e+06],
       [1.70386730e+02, 2.02510235e+06],
       [1.73994707e+02, 1.94968757e+06],
       [1.80849864e+02, 1.89386104e+06],
       [1.90591403e+02, 1.91051106e+06],
       [1.94560178e+02, 1.95752287e+06],
       [1.99972144e+02, 1.99082290e+06],
       [2.05023312e+02, 1.94087285e+06],
       [2.08992087e+02, 1.89777869e+06],
       [2.16568840e+02, 1.87231396e+06],
       [2.19816019e+02, 1.90757282e+06],
       [2.20898412e+02, 1.96241993e+06],
       [2.25588783e+02, 2.02118470e+06],
       [2.30279153e+02, 1.95556404e+06],
       [2.31722344e+02, 1.88308750e+06],
       [2.36773513e+02, 1.84195217e+06],
       [2.43628670e+02, 1.90561399e+06],
       [2.53009411e+02, 1.93303755e+06],
       [2.59503770e+02, 1.97417289e+06],
       [2.66358927e+02, 2.03979354e+06],
       [2.72492488e+02, 1.97613171e+06],
       [2.75739668e+02, 1.93695520e+06],
       [2.82955623e+02, 1.90071693e+06],
       [2.95222746e+02, 1.90365517e+06],
       [2.99191521e+02, 1.94185227e+06],
       [3.01717105e+02, 1.97906995e+06],
       [3.06407475e+02, 1.99865820e+06],
       [3.10737048e+02, 1.95948169e+06],
       [3.13984228e+02, 1.92716107e+06],
       [3.26972946e+02, 1.92716107e+06],
       [3.36353687e+02, 1.93107873e+06],
       [3.42848046e+02, 1.97711112e+06],
       [3.49342406e+02, 2.00453468e+06],
       [3.55475967e+02, 1.96731700e+06],
       [3.59444742e+02, 1.92422284e+06],
       [3.68103888e+02, 1.96829641e+06],
       [3.70990270e+02, 2.01530822e+06],
       [3.76041438e+02, 2.05252591e+06],
       [3.84339786e+02, 2.09170242e+06],
       [3.93359729e+02, 2.11227008e+06],
       [4.02018875e+02, 2.06623768e+06],
       [4.09234829e+02, 2.02902000e+06],
       [4.16089986e+02, 2.07701122e+06],
       [4.19697964e+02, 2.15046718e+06]])









vertices = np.asarray(vertices1)

#Save vertices as a pickle file
fuhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/upper_hybrid_id/GIRAFF_'+pld+'_HF_lower_boundary_vals_byeye.pkl'
#vertices = pickle.load(open(flhfile,'rb'))
pickle.dump([vertices], open(fuhfile,'wb'))




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

