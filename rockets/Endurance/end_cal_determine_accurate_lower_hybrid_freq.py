"""
Use phase and coherence spectra to determine accurate lower hybrid frequency 
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
from end_fields_loader import Endurance_Fields_Loader as EFL
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker
import pickle




#Load data for two channels of interest
#c1s = 'VLF24D'
#c2s = 'VLF13D'
c1s = 'VLF12D'
c2s = 'VLF34D'
#c1s = 'VLF24D'
#c2s = 'VLF32D'
#c1s = 'VLF13D'
#c2s = 'VLF32D'
#c1s = 'VLF24D'
#c2s = 'VLF41D'


c1 = EFL(c1s)
fs = c1.chnspecs['fs']
wf1, tdat = c1.load_data_gainphase_corrected()
c2 = EFL(c2s)
wf2, tdat = c2.load_data_gainphase_corrected()

#------see if this needs to be flipped
#wf2 = -1 * wf2
#cs2 = 'VLF31D'
#-------------------------------------

#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 4096
fspec, tspec, powerc1 = signal.spectrogram(wf1, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc2 = signal.spectrogram(wf2, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')

cohmin = 0.01  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result
#cohmin = 0.3 



##NOTE: + sense of phase defined as pointing towards center of potential of "powerc1"
##Nval = 3
#Nval = 3
#gx,cohx,phasex = correlation_analysis.interferometric_coherence_2D(powerc1,powerc2,Nval,coh_min=cohmin)
#phasex = np.degrees(phasex)


#tchunk = 0.5  #delta-time (sec) for each time chunk to divide up the spectra into (and average over)
tchunk = 2  #delta-time (sec) for each time chunk to divide up the spectra into (and average over)
nchunks = int(np.ceil((wf1.size/fs)/tchunk)) #number of chunks in ENTIRE timerange
cohx2, phasex2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wf1,wf2,tdat,fs,tchunk,coh_min=cohmin,nperseg=512)


phasex2 = np.degrees(phasex2)

ptmp_diff = np.abs(np.abs(powerc1) - np.abs(powerc2))
ptmp_sum = np.abs(powerc1) + np.abs(powerc2)
ptmp_fracdiff = ptmp_diff/ptmp_sum




mask = np.where(cohx2**2 > 0.6)




#yr = [0,15000]
yr = [4000,9000]
yscale='linear'
xr = [100,900]
#xr = [100,300]
vr = [-45,-30]


fig,axs = plt.subplots(3)
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.1,wspace=0.4)

ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=vr,xr=xr,yr=yr,yscale=yscale,ax=axs[0],
                    xlabel='time(sec)',ylabel=c1s + "\nfreq(Hz)")
axs[0].get_xaxis().set_visible(False)
ps.plot_spectrogram(tchunks2,freqs2,phasex2,vr=[-180,180], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,ax=axs[1],xlabel='time(sec)',
                    ylabel='phase')
axs[1].get_xaxis().set_visible(False)
ps.plot_spectrogram(tchunks2,freqs2,cohx2**2,vr=[-0.1,0.6], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,ax=axs[2],xlabel='time(sec)',
                    ylabel='coh**2\n(coh >'+str(cohmin)+')\nfreq(Hz)')



fig,axs = plt.subplots(1)
ps.plot_spectrogram(tchunks2,freqs2,cohx2**2,vr=[-0.1,0.6], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,ax=axs,xlabel='time(sec)',
                    ylabel='coh**2\n(coh >'+str(cohmin)+')\nfreq(Hz)')
klicker = clicker(axs, ["event"], markers=["x"])
vertices1 = (klicker.get_positions())['event']



fig,axs = plt.subplots(1)
ps.plot_spectrogram(tchunks2,freqs2,phasex2,vr=[-180,180], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',
                    ylabel='phase',ax=axs)

klicker = clicker(axs, ["event"], markers=["x"])
vertices1 = (klicker.get_positions())['event']


#from VLF12-VLF34 phase (May 12, 2025)
verticesLOW = [[110.1798756 , 6051.89455241],
       [ 112.4213161 , 6677.83915364],
       [ 114.28918318, 7243.59677398],
       [ 115.78347685, 7656.8807094 ],
       [ 124.00209201, 8078.18957561],
       [ 135.20929451, 8166.46381425],
       [ 144.92220334, 8050.10231787],
       [ 148.28436408, 7893.61616756],
       [ 154.26153875, 7717.06769029],
       [ 173.31378299, 7528.48181684],
       [ 188.25671965, 7367.98320114],
       [ 206.56181706, 7211.49705084],
       [ 222.62547397, 7115.19788142],
       [ 245.41345238, 6986.79898886],
       [ 262.59782954, 6886.48735405],
       [ 288.37439528, 6818.27544237],
       [ 306.30591927, 6734.01366913],
       [ 331.33533818, 6653.76436128],
       [ 353.37616975, 6589.564915  ],
       [ 381.02060257, 6549.44026108],
       [ 409.03860881, 6545.42779569],
       [ 435.56232139, 6513.32807255],
       [ 463.20675421, 6473.20341862],
       [ 488.23617311, 6477.21588401],
       [ 514.38631227, 6509.31560715],
       [ 540.53645143, 6497.27821098],
       [ 562.577283  , 6509.31560715],
       [ 581.62952724, 6537.4028649 ],
       [ 603.2967854 , 6577.52751883],
       [ 625.33761698, 6633.70203432],
       [ 646.2577283 , 6734.01366913],
       [ 666.05711938, 6830.31283855],
       [ 685.10936362, 6926.61200797],
       [ 705.65590153, 6982.78652347],
       [ 724.70814577, 7046.98596975],
       [ 749.36399126, 7151.31006995],
       [ 766.174795  , 7267.67156633],
       [ 782.23845191, 7460.26990517],
       [ 794.56637466, 7660.8931748 ],
       [ 804.6528569 , 7825.40425589],
       [ 812.12432523, 7965.84054463],
       [ 829.68227581, 8074.17711022],
       [ 839.39518464, 7949.79068306],
       [ 849.10809347, 7737.13001725],
       [ 854.3381213 , 7560.58153998],
       [ 858.44742888, 7339.8959434 ],
       [ 860.31529596, 6938.64940415],
       [ 864.79817696, 6569.50258804],
       [ 867.03961746, 6212.39316811],
       [ 870.77535163, 5831.20895582],
       [ 872.26964529, 5562.37377453],
       [ 875.63180604, 5261.43887009]]

"""
#from VLF12-VLF34 phase
verticesLOW = [[ 113.97920064, 7542.68022999],
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
"""

verticesHIG = [[ 118.55265075, 7964.38255622],
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



verticesLOW = np.asarray(verticesLOW)
verticesHIG = np.asarray(verticesHIG)

#Save vertices as a pickle file
flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/lower_hybrid_id/lower_hybrid_freqs_byeye.pkl'
#vertices = pickle.load(open(flhfile,'rb'))
pickle.dump([verticesLOW, verticesHIG], open(flhfile,'wb'))



fig, axs = plt.subplots(3)
yr=[4000,10000]
xr = [100,900]
vr = [-40,-35]

ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=vr,xr=xr,yr=yr,yscale=yscale,ax=axs[0],
                    xlabel='time(sec)',ylabel=c1s + "\nfreq(Hz)",title='end_cal_determine_accurate_lower_hybrid_freq.py')
ps.plot_spectrogram(tchunks2,freqs2,phasex2,vr=[-120,120], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',
                    ylabel='freq (Hz)',ax=axs[1])
ps.plot_spectrogram(tchunks2,freqs2,cohx2**2,vr=[0,0.6], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,ax=axs[2],xlabel='time(sec)',
                    ylabel='coh**2\n(coh >'+str(cohmin)+')\nfreq(Hz)')

#klicker = clicker(axs, ["event"], markers=["x"])
#verticestmp = (klicker.get_positions())['event']
#axs.plot(verticesLOW[:,0],verticesLOW[:,1],'*',color='magenta')
axs[0].plot(verticesHIG[:,0],verticesHIG[:,1],'*',color='magenta')
axs[1].plot(verticesHIG[:,0],verticesHIG[:,1],'*',color='magenta')
axs[2].plot(verticesHIG[:,0],verticesHIG[:,1],'*',color='magenta')
#axs.plot(verticesHIG[:,0],verticesHIG[:,1],color='magenta')

plt.savefig("/Users/abrenema/Desktop/tst.pdf", dpi=350)















fig,axs = plt.subplots(2)
ps.plot_spectrogram(tchunks2,freqs2,phasex2,vr=[-180,180], zscale='linear',
                    xr=xr,yr=yr,yscale=yscale,xlabel='time(sec)',
                    ylabel='phase',ax=axs[0],title='end_cal_determine_accurate_lower_hybrid_freq.py')
axs[0].plot(vertices[:,0],vertices[:,1],'*',color='magenta')
axs[0].plot(vertices2[:,0],vertices2[:,1],'*',color='black')
ps.plot_spectrogram(tspec,fspec,np.abs(powerc1),vr=vr,xr=xr,yr=yr,yscale=yscale,ax=axs[1],
                    xlabel='time(sec)',ylabel=c1s + "\nfreq(Hz)")
axs[1].plot(vertices[:,0],vertices[:,1],'*',color='magenta')
axs[1].plot(vertices2[:,0],vertices2[:,1],'*',color='black')

print('h')



