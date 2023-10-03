"""
Interferometry analysis on Endurance


Todo:
--Test different Bernstein bands and compare
--Test cross-correlation b/t one band and another on single VLF channel. 
--Test aliased Bernstein waves with Skins. This has advantages for doing interferometry
--Plot difference in PSD for long and short booms (see Pfaff+96). Any systematic difference
    can indicate that the shorter booms are not responding the same.
"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
from scipy.interpolate import interp1d
import plot_spectrogram as ps
import filter_wave_frequency as filt
import pickle
import correlation_analysis
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader



#--------------------------------------------------------------
#Get timeline of data to separate out science collection times 
#--------------------------------------------------------------

tl, gsS, gsE = end_data_loader.load_timeline()


#ephem = end_data_loader.load_ephemeris()

#---------------------------------------------
#Load gain/phase corrected data
#---------------------------------------------

v12DC = EFL('V12D')
v13DC = EFL('V13D')
v41DC = EFL('V41D')
v34DC = EFL('V34D')
v32DC = EFL('V32D')
v24DC = EFL('V24D')


wf12DC, tdatDC = v12DC.load_data_gainphase_corrected()
wf13DC, tdatDC = v13DC.load_data_gainphase_corrected()
wf34DC, tdatDC = v34DC.load_data_gainphase_corrected()
wf32DC, tdatDC = v32DC.load_data_gainphase_corrected()
wf24DC, tdatDC = v24DC.load_data_gainphase_corrected()
wf41DC, tdatDC = v41DC.load_data_gainphase_corrected()

fsDC = v12DC.chnspecs['fs']

v12 = EFL('VLF12D')
v13 = EFL('VLF13D')
v34 = EFL('VLF34D')
v24 = EFL('VLF24D')
v32 = EFL('VLF32D')
v41 = EFL('VLF41D')

wf12, tdat = v12.load_data_gainphase_corrected()
fs = v12.chnspecs['fs']
wf13, tgoo = v13.load_data_gainphase_corrected()
wf34, tgoo = v34.load_data_gainphase_corrected()
wf24, tgoo = v24.load_data_gainphase_corrected()
wf32, tgoo = v32.load_data_gainphase_corrected()
wf41, tgoo = v41.load_data_gainphase_corrected()


v1 = EFL('V1SD')
v2 = EFL('V2SD')
v3 = EFL('V3SD')
v4 = EFL('V4SD')

wf1, tdats = v1.load_data_gainphase_corrected()
wf2, tdats = v2.load_data_gainphase_corrected()
wf3, tdats = v3.load_data_gainphase_corrected()
wf4, tdats = v4.load_data_gainphase_corrected()


fss = v1.chnspecs['fs']



#----------------------------------------------------------------------
#Get spectral data for finding waves (use VLF12)
#----------------------------------------------------------------------

fspec, tspec, powerc12 = signal.spectrogram(wf12, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc13 = signal.spectrogram(wf13, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc34 = signal.spectrogram(wf34, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc32 = signal.spectrogram(wf32, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc24 = signal.spectrogram(wf24, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerc41 = signal.spectrogram(wf41, fs, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')

powerc = powerc12


fspecs, tspecs, powercs = signal.spectrogram(wf1, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs2, tspecs2, powercs2 = signal.spectrogram(wf1-wf2, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')

fspecDC, tspecDC, powerc12DC = signal.spectrogram(wf12DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powerc34DC = signal.spectrogram(wf34DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powerc32DC = signal.spectrogram(wf32DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecDC, tspecDC, powerc24DC = signal.spectrogram(wf24DC, fsDC, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')

fspecs, tspecs, powerc1s = signal.spectrogram(wf1, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs, tspecs, powerc2s = signal.spectrogram(wf2, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs, tspecs, powerc3s = signal.spectrogram(wf3, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')
fspecs, tspecs, powerc4s = signal.spectrogram(wf4, fss, nperseg=1024,noverlap=512,window='hann',return_onesided=True,mode='complex')


#--------------------------------------------------
#Test comparison to particle data 
#--------------------------------------------------


dfin, dfout = end_data_loader.load_PES()


#----------------------------------------------------------------------
#Compare all VLF waveforms for desired timerange
#----------------------------------------------------------------------

def plot_wf(tr, diagonals=1):

    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wf12bp[goot]), list(wf34bp[goot]), list(wf24bp[goot]), list(wf32bp[goot])])


    if diagonals:

        fig,axs = plt.subplots(6,figsize=(8,8))
        axs[0].plot(tdat[goot],wf13bp[goot],color='b')
        axs[0].plot(tdat[goot],wf24bp[goot],color='g')
        axs[0].title.set_text('VLF (13 vs 24)')
        axs[1].plot(tdat[goot],wf41bp[goot],color='b')
        axs[1].plot(tdat[goot],wf32bp[goot],color='r')
        axs[1].title.set_text('VLF (41 vs 32)')
        axs[2].plot(tdat[goot],wf13bp[goot],color='b')
        axs[2].plot(tdat[goot],wf32bp[goot],color='c')
        axs[2].title.set_text('VLF (13 vs 32)')
        axs[3].plot(tdat[goot],wf24bp[goot],color='g')
        axs[3].plot(tdat[goot],wf41bp[goot],color='r')
        axs[3].title.set_text('VLF (24 vs 41)')
        axs[4].plot(tdat[goot],wf24bp[goot],color='g')
        axs[4].plot(tdat[goot],wf32bp[goot],color='c')
        axs[4].title.set_text('VLF (24 vs 32)')
        axs[5].plot(tdat[goot],wf13bp[goot],color='r')
        axs[5].plot(tdat[goot],wf41bp[goot],color='c')
        axs[5].title.set_text('VLF (13 vs 41)')
        for i in range(6):
            axs[i].set_ylim(-maxv, maxv)

    else:
        fig,axs = plt.subplots(6,figsize=(8,8))
        axs[0].plot(tdat[goot],wf12bp[goot],color='b')
        axs[0].plot(tdat[goot],wf34bp[goot],color='g')
        axs[0].title.set_text('VLF (12 vs 34)')
        axs[1].plot(tdat[goot],wf12bp[goot],color='b')
        axs[1].plot(tdat[goot],wf24bp[goot],color='r')
        axs[1].title.set_text('VLF (12 vs 24)')
        axs[2].plot(tdat[goot],wf12bp[goot],color='b')
        axs[2].plot(tdat[goot],wf32bp[goot],color='c')
        axs[2].title.set_text('VLF (12 vs 32)')
        axs[3].plot(tdat[goot],wf34bp[goot],color='g')
        axs[3].plot(tdat[goot],wf24bp[goot],color='r')
        axs[3].title.set_text('VLF (34 vs 24)')
        axs[4].plot(tdat[goot],wf34bp[goot],color='g')
        axs[4].plot(tdat[goot],wf32bp[goot],color='c')
        axs[4].title.set_text('VLF (34 vs 32)')
        axs[5].plot(tdat[goot],wf24bp[goot],color='r')
        axs[5].plot(tdat[goot],wf32bp[goot],color='c')
        axs[5].title.set_text('VLF (24 vs 32)')
        for i in range(6):
            axs[i].set_ylim(-maxv, maxv)



#----------------------------------------------------------------------
#Compare all hodograms for desired timerange
#----------------------------------------------------------------------

def plot_hod(tr):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    maxv = np.max([list(wf12bp[goot]), list(wf34bp[goot]), list(wf24bp[goot]), list(wf32bp[goot])])

    fig,axs = plt.subplots(3,2,figsize=(8,8))
    axs[0,0].plot(wf12bp[goot],wf34bp[goot])    
    axs[0,1].plot(wf12bp[goot],wf24bp[goot])    
    axs[1,0].plot(wf12bp[goot],wf32bp[goot])    
    axs[1,1].plot(wf34bp[goot],wf24bp[goot])    
    axs[2,0].plot(wf34bp[goot],wf32bp[goot])    
    axs[2,1].plot(wf24bp[goot],wf32bp[goot])    
    axs[0,0].title.set_text('12 vs 34')
    axs[0,1].title.set_text('12 vs 24')
    axs[1,0].title.set_text('12 vs 32')
    axs[1,1].title.set_text('34 vs 24')
    axs[2,0].title.set_text('34 vs 32')
    axs[2,1].title.set_text('24 vs 32')
    for i in range(3):
        for j in range(2):
            axs[i,j].set(aspect='equal')
            axs[i,j].set_xlim(-maxv, maxv)
            axs[i,j].set_ylim(-maxv, maxv)



#Cross-spectral density for skins 
def plot_csd_skins(tr, xlm=[0,1000],ylm=[-180,180],nps=512,mincoh=0.5,ylmpsd=[0,0]):

    goot = np.where((tdats >= tr[0]) & (tdats <= tr[1]))
    window = np.hanning(len(goot[0]))
    wf1z = wf1[goot]*window
    wf2z = wf2[goot]*window
    wf3z = wf3[goot]*window
    wf4z = wf4[goot]*window
    
    csd12, angle12, f = correlation_analysis.cross_spectral_density(wf1z,wf2z,fss,nperseg=nps)
    #csd34, angle34, f = correlation_analysis.cross_spectral_density(wf3z,wf4z,fss,nperseg=nps)
 
    goodv12 = np.where(csd12 >= mincoh)
    #goodv34 = np.where(csd34 >= mincoh)


    psdS1, psdf1 = correlation_analysis.psd(wf1z, tdats[goot], fss, tr)
    psdS2, psdf2 = correlation_analysis.psd(wf2z, tdats[goot], fss, tr)


    fig, axs = plt.subplots(3,figsize=(12,8))
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.4,
                        hspace=0.4)

    axs[0].plot(f, csd12)
    axs[0].plot(f[goodv12], csd12[goodv12],'.')
    axs[0].set_ylabel('coherency\nV1,V2')
    axs[0].set_ylim(0,1)
    axs[0].set_xlabel('freq(Hz)')
    axs[1].plot(f, angle12, color='r')
    axs[1].plot(f[goodv12], angle12[goodv12],'.')
    axs[1].set_ylim(ylm)
    axs[1].set_ylabel('Phase\nV1-V2')
    axs[2].plot(psdf1, psdS1, psdf2, psdS2)
    axs[2].set_ylabel('PSD')
    axs[2].set_yscale('log')

    for i in range(3):
        axs[i].set_xlim(xlm)
    
    if ylmpsd != [0,0]:
        axs[2].set_ylim(ylmpsd)



#Cross-spectral density for VLF
def plot_csd(tr, xlm=[0,12000],ylm=[-180,180],nps=512,mincoh=0.5,diagonals=1):
    goot = np.where((tdat >= tr[0]) & (tdat <= tr[1]))
    window = np.hanning(len(goot[0]))
    wf12z = wf12[goot]*window
    wf13z = wf13[goot]*window
    wf24z = wf24[goot]*window
    wf34z = wf34[goot]*window
    wf32z = wf32[goot]*window
    wf41z = wf41[goot]*window
    csd1234, angle1234, f = correlation_analysis.cross_spectral_density(wf12z,wf34z,fs,nperseg=nps)
    csd1213, angle1213, f = correlation_analysis.cross_spectral_density(wf12z,wf13z,fs,nperseg=nps)
    csd1224, angle1224, f = correlation_analysis.cross_spectral_density(wf12z,wf24z,fs,nperseg=nps)
    csd1232, angle1232, f = correlation_analysis.cross_spectral_density(wf12z,wf32z,fs,nperseg=nps)
    csd1241, angle1241, f = correlation_analysis.cross_spectral_density(wf12z,wf41z,fs,nperseg=nps)

    csd3413, angle3413, f = correlation_analysis.cross_spectral_density(wf34z,wf13z,fs,nperseg=nps)
    csd3424, angle3424, f = correlation_analysis.cross_spectral_density(wf34z,wf24z,fs,nperseg=nps)
    csd3441, angle3441, f = correlation_analysis.cross_spectral_density(wf34z,wf41z,fs,nperseg=nps)
    csd3432, angle3432, f = correlation_analysis.cross_spectral_density(wf34z,wf32z,fs,nperseg=nps)

    csd1324, angle1324, f = correlation_analysis.cross_spectral_density(wf13z,wf24z,fs,nperseg=nps)
    csd1341, angle1341, f = correlation_analysis.cross_spectral_density(wf13z,wf41z,fs,nperseg=nps)
    csd1332, angle1332, f = correlation_analysis.cross_spectral_density(wf13z,wf32z,fs,nperseg=nps)

    csd2441, angle2441, f = correlation_analysis.cross_spectral_density(wf24z,wf41z,fs,nperseg=nps)
    csd2432, angle2432, f = correlation_analysis.cross_spectral_density(wf24z,wf32z,fs,nperseg=nps)

    csd4132, angle4132, f = correlation_analysis.cross_spectral_density(wf41z,wf32z,fs,nperseg=nps)



    goodv1234 = np.where(csd1234 >= mincoh)
    goodv1213 = np.where(csd1213 >= mincoh)
    goodv1224 = np.where(csd1224 >= mincoh)
    goodv1232 = np.where(csd1232 >= mincoh)
    goodv1241 = np.where(csd1241 >= mincoh)
    goodv3413 = np.where(csd3413 >= mincoh)
    goodv3424 = np.where(csd3424 >= mincoh)
    goodv3441 = np.where(csd3441 >= mincoh)
    goodv3432 = np.where(csd3432 >= mincoh)
    goodv1324 = np.where(csd1324 >= mincoh)
    goodv1341 = np.where(csd1341 >= mincoh)
    goodv1332 = np.where(csd1332 >= mincoh)
    goodv2441 = np.where(csd2441 >= mincoh)
    goodv2432 = np.where(csd2432 >= mincoh)
    goodv4132 = np.where(csd4132 >= mincoh)




    #Plot the diagonal pairs only
    if diagonals:

        fig, axs = plt.subplots(4,3,figsize=(12,8))
        plt.subplots_adjust(left=0.1,
                            bottom=0.1,
                            right=0.9,
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)

        axs[0,0].plot(f, csd1324)
        axs[0,0].plot(f[goodv1324], csd1324[goodv1324],'.')
        axs[0,0].set_ylabel('coherency\nVLF13-VLF24')
        axs[0,0].set_xlabel('freq(Hz)')
        axs[1,0].plot(f, angle1324, color='r')
        axs[1,0].plot(f[goodv1324], angle1324[goodv1324],'.')
        axs[1,0].set_ylim(ylm)
        axs[1,0].set_ylabel('Phase\nVLF13-VLF24')

        axs[0,1].plot(f, csd4132)
        axs[0,1].plot(f[goodv4132], csd4132[goodv4132],'.')
        axs[0,1].set_ylabel('coherency\nVLF41-VLF32')
        axs[0,1].set_xlabel('freq(Hz)')
        axs[1,1].plot(f, angle4132, color='r')
        axs[1,1].plot(f[goodv4132], angle4132[goodv4132],'.')
        axs[1,1].set_ylim(ylm)
        axs[1,1].set_ylabel('Phase\nVLF41-VLF32')

        axs[0,2].plot(f, csd1332)
        axs[0,2].plot(f[goodv1332], csd1332[goodv1332],'.')
        axs[0,2].set_ylabel('coherency\nVLF13-VLF32')
        axs[0,2].set_xlabel('freq(Hz)')
        axs[1,2].plot(f, angle1332, color='r')
        axs[1,2].plot(f[goodv1332], angle1332[goodv1332],'.')
        axs[1,2].set_ylim(ylm)
        axs[1,2].set_ylabel('Phase\nVLF13-VLF32')

        axs[2,0].plot(f, csd2441)
        axs[2,0].plot(f[goodv2441], csd2441[goodv2441],'.')
        axs[2,0].set_ylabel('coherency\nVLF24-VLF41')
        axs[2,0].set_xlabel('freq(Hz)')
        axs[3,0].plot(f, angle2441, color='r')
        axs[3,0].plot(f[goodv2441], angle2441[goodv2441],'.')
        axs[3,0].set_ylim(ylm)
        axs[3,0].set_ylabel('Phase\nVLF24-VLF41')

        axs[2,1].plot(f, csd2432)
        axs[2,1].plot(f[goodv2432], csd2432[goodv2432],'.')
        axs[2,1].set_ylabel('coherency\nVLF24-VLF32')
        axs[2,1].set_xlabel('freq(Hz)')
        axs[3,1].plot(f, angle2432, color='r')
        axs[3,1].plot(f[goodv2432], angle2432[goodv2432],'.')
        axs[3,1].set_ylim(ylm)
        axs[3,1].set_ylabel('Phase\nVLF24-VLF32')

        axs[2,2].plot(f, csd1341)
        axs[2,2].plot(f[goodv1341], csd1341[goodv1341],'.')
        axs[2,2].set_ylabel('coherency\nVLF13-VLF41')
        axs[2,2].set_xlabel('freq(Hz)')
        axs[3,2].plot(f, angle1341, color='r')
        axs[3,2].plot(f[goodv1341], angle1341[goodv1341],'.')
        axs[3,2].set_ylim(ylm)
        axs[3,2].set_ylabel('Phase\nVLF13-VLF41')


        for i in range(4):
            for j in range(3):
                axs[i,j].set_xlim(xlm)
                #axs[i,j].set_ylim(ylm)



    #Plot original channels (w/o 13 and 41)
    else: 

        fig, axs = plt.subplots(4,3,figsize=(12,8))
        plt.subplots_adjust(left=0.1,
                            bottom=0.1,
                            right=0.9,
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)

        axs[0,0].plot(f, csd1234)
        axs[0,0].plot(f[goodv1234], csd1234[goodv1234],'.')
        axs[0,0].set_ylabel('coherency\nVLF12-VLF34')
        axs[0,0].set_xlabel('freq(Hz)')
        axs[1,0].plot(f, angle1234, color='r')
        axs[1,0].plot(f[goodv1234], angle1234[goodv1234],'.')
        axs[1,0].set_ylim(ylm)
        axs[1,0].set_ylabel('Phase\nVLF12-VLF34')

        axs[0,1].plot(f, csd1224)
        axs[0,1].plot(f[goodv1224], csd1224[goodv1224],'.')
        axs[0,1].set_ylabel('coherency\nVLF12-VLF24')
        axs[0,1].set_xlabel('freq(Hz)')
        axs[1,1].plot(f, angle1224, color='r')
        axs[1,1].plot(f[goodv1224], angle1224[goodv1224],'.')
        axs[1,1].set_ylim(ylm)
        axs[1,1].set_ylabel('Phase\nVLF12-VLF24')

        axs[0,2].plot(f, csd1232)
        axs[0,2].plot(f[goodv1232], csd1232[goodv1232],'.')
        axs[0,2].set_ylabel('coherency\nVLF12-VLF32')
        axs[0,2].set_xlabel('freq(Hz)')
        axs[1,2].plot(f, angle1232, color='r')
        axs[1,2].plot(f[goodv1232], angle1232[goodv1232],'.')
        axs[1,2].set_ylim(ylm)
        axs[1,2].set_ylabel('Phase\nVLF12-VLF32')

        axs[2,0].plot(f, csd3424)
        axs[2,0].plot(f[goodv3424], csd3424[goodv3424],'.')
        axs[2,0].set_ylabel('coherency\nVLF34-VLF24')
        axs[2,0].set_xlabel('freq(Hz)')
        axs[3,0].plot(f, angle3424, color='r')
        axs[3,0].plot(f[goodv3424], angle3424[goodv3424],'.')
        axs[3,0].set_ylim(ylm)
        axs[3,0].set_ylabel('Phase\nVLF34-VLF24')

        axs[2,1].plot(f, csd3432)
        axs[2,1].plot(f[goodv3432], csd3432[goodv3432],'.')
        axs[2,1].set_ylabel('coherency\nVLF34-VLF32')
        axs[2,1].set_xlabel('freq(Hz)')
        axs[3,1].plot(f, angle3432, color='r')
        axs[3,1].plot(f[goodv3432], angle3432[goodv3432],'.')
        axs[3,1].set_ylim(ylm)
        axs[3,1].set_ylabel('Phase\nVLF34-VLF32')

        axs[2,2].plot(f, csd2432)
        axs[2,2].plot(f[goodv2432], csd2432[goodv2432],'.')
        axs[2,2].set_ylabel('coherency\nVLF24-VLF32')
        axs[2,2].set_xlabel('freq(Hz)')
        axs[3,2].plot(f, angle2432, color='r')
        axs[3,2].plot(f[goodv2432], angle2432[goodv2432],'.')
        axs[3,2].set_ylim(ylm)
        axs[3,2].set_ylabel('Phase\nVLF24-VLF32')


        for i in range(4):
            for j in range(3):
                axs[i,j].set_xlim(xlm)
                #axs[i,j].set_ylim(ylm)



#Makes nice plot of coh/phase spectrograms and resulting slice in coherence/phase (including average)
def plot_phase_coh(parr,carr,pharr,powa,coha,pha,xrspec=[100,200],yrspec=[5500,7000],
                   vr1=[-40,-25],vr2=[0.9,1],vr3=[0,140],
                   xrline=[5500,7000],
                   ylm1=[0,1],ylm2=[0,1],ylm3=[-180,180],
                   tzp1=0,tzp2=0):

    try:
        figgoo = plt.gcf()
    except NameError:
        print("No fig yet defined")
    else:
        plt.close(figgoo)


    fig,axs = plt.subplots(6, gridspec_kw={'height_ratios':[1,1,1,1,1,4]})
    ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=vr1,xr=xrspec,yr=yrspec,yscale='linear',ax=axs[0])
    ps.plot_spectrogram(tchunks,freqs,coh,vr=vr2, zscale='linear',xr=xrspec,yr=yrspec,yscale='linear',ax=axs[1])
    #ps.plot_spectrogram(tchunks,freqs,np.abs(phase),vr=vr3, zscale='linear',xr=xrspec,yr=yrspec,yscale='linear',ax=axs[2])
    ps.plot_spectrogram(tchunks,freqs,phase,vr=vr3, zscale='linear',xr=xrspec,yr=yrspec,yscale='linear',ax=axs[2])
    for i in range(3): axs[i].axvline(tzp1)
    for i in range(3): axs[i].axvline(tzp2)
    
    axs[3].plot(fspec,parr)
    axs[3].plot(fspec,powa,'.',color='black')
    axs[3].set_xlim(xrline)
    axs[3].set_ylim(ylm1)
    
    axs[4].plot(freqs,carr)
    axs[4].plot(freqs,coha,'.',color='black')
    axs[4].set_xlim(xrline)
    axs[4].set_ylim(ylm2)
    
    axs[5].plot(freqs,pharr)
    axs[5].plot(freqs,pha,'.',color='black')
    axs[5].set_xlim(xrline)
    axs[5].set_ylim(ylm3)




#----------------------------------------------------------------------------
#Bernstein at mission start (phase and coherence test)
#----------------------------------------------------------------------------

tchunk = 0.2  #sec
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,-wf24,tdat,fs,tchunk,coh_min=0.4,nperseg=2048)
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf32,-wf41,tdat,fs,tchunk,coh_min=0.,nperseg=2048)
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.,nperseg=512)
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(-wf41,wf32,tdat,fs,tchunk,coh_min=0.,nperseg=512)


fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-40,-25],yr=[5000,8000],xr=[110,220], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=[-40,-25],yr=[5000,8000],xr=[110,220], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc24),vr=[-40,-25],yr=[5000,8000],xr=[110,220], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc32),vr=[-40,-25],yr=[5000,8000],xr=[110,220], yscale='linear',ax=axs[3])

tslice = 140
nsec = 20
p12avg, p12vals, tv = ps.slice_spectrogram(tslice,tspec,np.abs(powerc12),nsec)
p34avg, p34vals, tv = ps.slice_spectrogram(tslice,tspec,np.abs(powerc34),nsec)
p24avg, p24vals, tv = ps.slice_spectrogram(tslice,tspec,np.abs(powerc24),nsec)
p32avg, p32vals, tv = ps.slice_spectrogram(tslice,tspec,np.abs(powerc32),nsec)


fig,axs = plt.subplots(1)
axs.plot(fspec,np.sqrt(p12avg),color='black')
axs.plot(fspec,np.sqrt(p34avg),color='blue')
axs.plot(fspec,np.sqrt(p24avg),color='red')
axs.plot(fspec,np.sqrt(p32avg),color='orange')
axs.set_xlim(5500,7500)
axs.set_xlim(0,7500)



goob = list(range(0,45))
goob = goob[1:]
thetarad = [i * (np.pi/180) for i in goob]
l45 = 45 * (np.pi/180)
rat = [np.cos(l45 + i)/np.cos(l45 - i) for i in thetarad]
plt.plot(goob, rat)
plt.yscale('linear')
plt.ylim(0,1)



#Slice the phase vs freq for particular time(s)
#tz = 140
#tz = 129.7
#tz = 183
#nsec = 8
#tz = 160
#nsec = 8
tz = 140
nsec = 2
pavg, phasearr, tarr_phase = ps.slice_spectrogram(tz,tchunks,phase,nsec)
cavg, coharr, tarr_coh = ps.slice_spectrogram(tz,tchunks,coh,nsec)
powavg, powarr, tarr_pow = ps.slice_spectrogram(tz,tspec,np.abs(powerc32),nsec)

tr = [120,220]
#yr = [5000,7000]
yr = [5000,9000]
plot_phase_coh(powarr,coharr,phasearr,powavg,cavg,pavg,xrspec=tr,yrspec=yr,
               vr1=[-60,-20],vr2=[0.9,1],vr3=[0,140],
               xrline=yr,
               ylm1=[0,0.002],ylm2=[0,1],ylm3=[-180,180],
               tzp1=tz,tzp2=tz+nsec)


#-------------------------------------------------------
#-------------------------------------------------------
#UNDER CONSTRUCTION 
#-------------------------------------------------------
#-------------------------------------------------------

#Make plot of E-field power (z-axis) vs f (y-axis) and k-value (x-axis) [Lalti+23 MMS paper; https://doi.org/10.1029/2022JA031150]
#Do this with nested for loop. For each freq (all times of interest) bin the freq vs k-values.  
#k-values come from the phase as k = phase (radians) / d 




nkbins = 100
klim = [-2,2]
kstep = (klim[1]-klim[0])/nkbins
kvals = np.arange(klim[0],klim[1],kstep)   #k-values to bin\
d = 2.26 #meters - effective length of spaced receiver

nfreqs = np.shape(phasearr)[0]

tdelta = tarr_phase[1] - tarr_phase[0]


#final values will be [nfreqs,nkbins]
pow_finv = np.empty((nfreqs,nkbins))
pow_finv2 = np.empty((nfreqs,nkbins))


for f in range(0,nfreqs-1,1):
    pslice = phasearr[f,:] #phase values for a particular freq and all times
    kslice = [(3.14/360)*i/d for i in pslice] #Change delta-phases into k-values

    if np.sum(kslice != 0):
        #for current freq slice bin the k-values for all times
        for k in range(0,nkbins-1,1):

            #If, at current frequency slice there are k-values within the current range of k, extract the corresponding power values.
            gootimeIDX = np.where((kslice >= kvals[k]) & (kslice < kvals[k+1])) #time indices satisfying condition
            #--Note that the power array is a different size than the phase array, so need to look over a range of values
            powavg_goo = np.empty(len(gootimeIDX[0])) #Can be multiple times for current freq that have k-values in current range (b/c we're considering all times)

            if len(gootimeIDX[0]) != 0:

                #relevant freq range of (larger) power array to average over 
                frange_goo = [freqs[f],freqs[f+1]]
                frangeIDX = np.where((fspec >= frange_goo[0]) & (fspec <= frange_goo[1]))[0]
                frangeIDX = [np.min(frangeIDX), np.max(frangeIDX)]


                for t in range(0,len(gootimeIDX[0]),1):
                    #relevant time range of (larger) power array to average over 
                    if t < len(gootimeIDX[0])-1:
                        trange_goo = [tarr_phase[t],tarr_phase[t+1]]
                    else:
                        trange_goo = [tarr_phase[t],tarr_phase[t]+tdelta]

                    trangeIDX = np.where((tarr_pow >= trange_goo[0]) & (tarr_pow <= trange_goo[1]))[0]
                    trangeIDX = [np.min(trangeIDX), np.max(trangeIDX)]

                    powarrtmp = powarr[frangeIDX[0]:frangeIDX[1],trangeIDX[0]:trangeIDX[1]]
                    powavg_goo[t] = np.mean(powarrtmp) #average over all the overlapping freq and time range of larger power array

                #if np.sum(powavg_goo != 0):

                if np.sum(powavg_goo > 1e-100):
                    pow_finv[f,k] = np.mean(powavg_goo)            
                    pow_finv2[f,k] = np.nanmax(powavg_goo)            
                print('here')

#change zero values to different number (conflicts with dB scale)

print("end of nested loop")


#Plot the results (freq vs k)
vr = [-40,-25]
yr = [4000,8000]
kr = [-1,1]

fig,axs = plt.subplots(2)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=vr,yr=yr,xr=[110,220], yscale='linear',ax=axs[0],xlabel='time(s)',ylabel='f(Hz)')
ps.plot_spectrogram(kvals,freqs,pow_finv2,vr=vr,xr=kr,yr=yr,yscale='linear',ax=axs[1],minzval=-120,xlabel='k(1/m)',ylabel='f(Hz)')

#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------







#wf12-wf32 results for 183s
df = (6200 - 6000)   #Hz
dp = (72 - 44.4) * (3.14/180)    #rad
#L = 0.446  m
#wf12-wf42 results for 183s
df = (6200 - 6000)   #Hz
dp = (75.8 - 53.4) * (3.14/180)    #rad
#L = 0.55 m


#wf34-wf32 results for 183s
df = (6200 - 6000)   #Hz
dp = (-23.09 + 26.3) * (3.14/180)    #rad
#L = 3 m  (not well determined)



#wf12-wf32 results
df = (6607 - 6064)   #Hz
dp = (72.9 - 27.3) * (3.14/180)    #rad



#phase speed in m/s
slope = dp/df  #rad * sec
d = 1.06 #separation vector in meters
vp = 2*np.pi*d/slope
ftst = 6200 #test freq in Hz
#wave vector in rad/m 
k = (2*np.pi*ftst) / vp
#wavelength in meters
wavelength = 2*np.pi/k  #L(6200 Hz) = 0.73 m




#----------------------------------------------------------------------------
#Bernstein at mission end (phase and coherence test)
#----------------------------------------------------------------------------

ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-40,-25],yr=[4000,8000],xr=[750,850], yscale='linear')


#First clue to short wavelength is that 12 and 34 see less power than 32 and 42
fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-40,-25],yr=[4000,8000],xr=[750,850], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=[-40,-25],yr=[4000,8000],xr=[750,850], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc24),vr=[-40,-25],yr=[4000,8000],xr=[750,850], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc32),vr=[-40,-25],yr=[4000,8000],xr=[750,850], yscale='linear',ax=axs[3])


tslice = 790
nsec = 20
p12avg, p12vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc12),nsec)
p13avg, p13vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc13),nsec)
p34avg, p34vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc34),nsec)
p24avg, p24vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc24),nsec)
p32avg, p32vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc32),nsec)
p41avg, p41vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc41),nsec)


fig,axs = plt.subplots(1)
axs.plot(fspec,np.sqrt(p12avg),color='black')
axs.plot(fspec,np.sqrt(p13avg),color='purple')
axs.plot(fspec,np.sqrt(p34avg),color='blue')
axs.plot(fspec,np.sqrt(p24avg),color='red')
axs.plot(fspec,np.sqrt(p32avg),color='orange')
axs.plot(fspec,np.sqrt(p41avg),color='green')
axs.set_xlim(1000,10000)




tchunk = 2  #sec
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.5,nperseg=2048)
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf13,wf32,tdat,fs,tchunk,coh_min=0.5,nperseg=2048)



#Slice the phase vs freq for particular time(s)
tz = 820
nsec = 5
pavg, phasearr = ps.slice_spectrogram(tz,tchunks,phase,nsec)
cavg, coharr = ps.slice_spectrogram(tz,tchunks,coh,nsec)
powavg, powarr = ps.slice_spectrogram(tz,tspec,np.abs(powerc13),nsec)

tr = [800,840]
yr = [4000,8000]
plot_phase_coh(xrspec=tr,yrspec=[4000,7500],
               vr1=[-50,-20],vr2=[0.4,1],vr3=[-100,100],
               xrline=yr,
               ylm1=[0,0.002],ylm2=[0,1],ylm3=[-180,180],
               tzp1=tz,tzp2=tz+nsec)


print('h')


#----------------------------------------
#Lower hybrid wave (phase and coherence test)
#----------------------------------------

fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powerc12),vr=[-40,-25],yr=[4000,12000],xr=[645,647], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc34),vr=[-40,-25],yr=[4000,12000],xr=[645,647], yscale='linear',ax=axs[1])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc24),vr=[-40,-25],yr=[4000,12000],xr=[645,647], yscale='linear',ax=axs[2])
ps.plot_spectrogram(tspec,fspec,np.abs(powerc32),vr=[-40,-25],yr=[4000,12000],xr=[645,647], yscale='linear',ax=axs[3])


#tslice = 646.75
#nsec = 0.2
tslice = 645
nsec = 2
p12avg, p12vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc12),nsec)
p13avg, p13vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc13),nsec)
p34avg, p34vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc34),nsec)
p24avg, p24vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc24),nsec)
p32avg, p32vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc32),nsec)
p41avg, p41vals = ps.slice_spectrogram(tslice,tspec,np.abs(powerc41),nsec)

fig,axs = plt.subplots(1)
axs.plot(fspec,np.sqrt(p12avg),color='black')
axs.plot(fspec,np.sqrt(p13avg),color='purple')
axs.plot(fspec,np.sqrt(p34avg),color='blue')
axs.plot(fspec,np.sqrt(p24avg),color='red')
axs.plot(fspec,np.sqrt(p32avg),color='orange')
axs.plot(fspec,np.sqrt(p41avg),color='green')
axs.set_xlim(0,12000)




tchunk = 1  #sec
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.3,nperseg=512)
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf13,wf24,tdat,fs,tchunk,coh_min=0.3,nperseg=512)

#Slice the phase vs freq for particular time(s)
tz = 644
nsec = 8
#tz = 647
#nsec = 0.3
pavg, phasearr = ps.slice_spectrogram(tz,tchunks,phase,nsec)
cavg, coharr = ps.slice_spectrogram(tz,tchunks,coh,nsec)
powavg, powarr = ps.slice_spectrogram(tz,tspec,np.abs(powerc12),nsec)


tr = [630,660]
yr = [5000,12000]
plot_phase_coh(xrspec=tr,yrspec=yr,
               vr1=[-50,-30],vr2=[0.9,1],vr3=[0,140],
               xrline=[4000,12000],
               ylm1=[0,0.002],ylm2=[0,1],ylm3=[-200,200],
               tzp1=tz,tzp2=tz+nsec)



#wf12-wf32 results
df = (12000-5000)   #Hz
dp = (9.5 - 7.1) * (3.14/180)    #rad


#phase speed in m/s
slope = dp/df  #rad * sec
d = 1.06 #separation vector in meters
vp = 2*np.pi*d/slope   #~1000 km/s
ftst = 8000 #test freq in Hz
#wave vector in rad/m 
k = (2*np.pi*ftst) / vp
#wavelength in meters
wavelength = 2*np.pi/k  #L(8000 Hz) = 140 m


#------------------------------------------------------------------------------------
#H+ fundamental  (phase and coherence test) 
# - narrowband but I'm able to squeeze out enough info to get a slope. 
#------------------------------------------------------------------------------------

tchunk = 0.5  #sec
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.3,nperseg=2048)
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf12,wf32,tdat,fs,tchunk,coh_min=0.3,nperseg=1024)
#coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf13,wf24,tdat,fs,tchunk,coh_min=0.3,nperseg=1024)
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf13,wf32,tdat,fs,tchunk,coh_min=0.3,nperseg=1024)

#Slice the phase vs freq for particular time(s)
#tz = 498
#nsec = 15
#tz = 502
#nsec = 3
tz = 505
nsec = 6
pavg, phasearr = ps.slice_spectrogram(tz,tchunks,phase,nsec)
cavg, coharr = ps.slice_spectrogram(tz,tchunks,coh,nsec)
powavg, powarr = ps.slice_spectrogram(tz,tspec,np.abs(powerc12),nsec)


tr = [450,550]
yr = [4000,4500]
plot_phase_coh(xrspec=tr,yrspec=yr,
               vr1=[-50,-30],vr2=[0.3,1],vr3=[0,140],
               xrline=[4000,4500],
               ylm1=[0,0.002],ylm2=[0,1],ylm3=[-200,200],
               tzp1=tz,tzp2=tz+nsec)

#Rough slope 
df = 4248 - 4204 
dp = 88.5 - 53.8 




#--------------------------------------
#Wave TEST - long wavelength LH 
#--------------------------------------

#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-10],yr=[6000,12000],xr=[430,432], yscale='linear')
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-20],yr=[6000,12000],xr=[800,810], yscale='linear')
fmin = 6000
fmax = 12000
wf12bp = filt.butter_bandpass_filter(wf12, fmin, fmax, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, fmin, fmax, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, fmin, fmax, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, fmin, fmax, fs, order= 10)


tr = [430.8222,430.825] 
#tr = [802.003,802.006] 
plot_wf(tr)



#--------------------------------------
#Wave 1: Berstein (~120-130 sec at 5000-8000 Hz)
#--------------------------------------

ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='linear')


fmin = 4000
fmax = 6000
wf12bp = filt.butter_bandpass_filter(wf12, fmin, fmax, fs, order= 10)
wf13bp = filt.butter_bandpass_filter(wf13, fmin, fmax, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, fmin, fmax, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, fmin, fmax, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, fmin, fmax, fs, order= 10)
wf41bp = filt.butter_bandpass_filter(wf41, fmin, fmax, fs, order= 10)

fspecz, tspecz, powercz = signal.spectrogram(wf12bp, fs, nperseg=4*2048,noverlap=2*2048,window='hann',return_onesided=True,mode='complex')
ps.plot_spectrogram(tspecz,fspecz,np.abs(powercz),vr=[-60,-30],yr=[4000,8000],xr=[120,130], yscale='log')



#Snapshot 1
tr = [122.62,122.72]  #zoomed out
#tr = [122.67718,122.68]  #zoomed in
plot_wf(tr)

plt.plot(tdat,wf12bp,tdat,-wf24bp)
plt.xlim(122.67775,122.6800)
plt.ylim(-0.05,0.05)



#tr = [122.62,122.72]  #zoomed out
tr = [124,130]  #zoomed out
S34, f = correlation_analysis.psd(wf34bp, tdat, fs, tr)
S32, f = correlation_analysis.psd(wf32bp, tdat, fs, tr)

plt.semilogy(f,S34, f, S32)
plt.ylim(1e-12, 1e-4)
#plt.xscale('log')
plt.xscale('linear')
plt.xlim(10,12000)
plt.xlim(4000,8000)


#import scipy.signal
#(f, S) = scipy.signal.periodogram(wf, fs, scaling='density')


plot_csd([tr[0]-3, tr[1]+3],nps=256)
plot_hod(tr) #Fairly circular



#Snapshot 2
tr = [123.5,123.54]
tr = [123.522,123.526] #Zoomed in
plot_wf(tr)
plot_csd([tr[0]-0.2, tr[1]+0.2])
plot_hod(tr)




plot_wf([123.00275,123.00375])
plot_wf([115.00275,128.00375])
plot_hod([123.00275,123.00375])






#NOTES: angle for Bernstein waves (wf12, wf34) is > 0 and increases towards 90. 
#Assuming wave vector propagates perp to Bo, this implies wavelengths < 5m * sin(90) = 5 m [Kintner+98]. 
#The DC power is at phase = 0, which means long wavelength (or large structures)


#Dynamic cross-spectral density
goot = np.where((tdat >= 100) & (tdat <= 200))
window = np.hanning(len(goot[0]))
wf1z = wf12[goot]*window
wf2z = wf34[goot]*window
timechunk = 1  #sec
Pxy, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wf1z,wf2z,tdat[goot],fs,timechunk,nperseg=1024,plot=True)






print("here")











#--------------------------------------
#Wave 3: 165.46 - 165.56 sec; ~100-1000 Hz
#--------------------------------------

wf12bp = filt.butter_bandpass_filter(wf12, 60, 300, fs, order= 10)
wf34bp = filt.butter_bandpass_filter(wf34, 60, 300, fs, order= 10)
wf24bp = filt.butter_bandpass_filter(wf24, 60, 300, fs, order= 10)
wf32bp = filt.butter_bandpass_filter(wf32, 60, 300, fs, order= 10)

plot_wf([165.20,165.8])
plot_hod([165.20,165.8])



#-------------------------------------------------------------------------------------
#Density structures
#--can be somewhat difficult to identify in skin data
#-------------------------------------------------------------------------------------

ps.plot_spectrogram(tspec12DC,fspec12DC,np.abs(powerc12DC),vr=[-60,10],yr=[10,1000],xr=[100,850], yscale='log')


fig, axs = plt.subplots(2)
ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[10,100],xr=[630,660], yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecs2,fspecs2,np.abs(powercs2),vr=[-60,-40],yr=[10,100],xr=[630,660], yscale='linear',ax=axs[1])
#ps.plot_spectrogram(tspec,fspec,np.abs(powerc),vr=[-60,-10],yr=[10,1000],xr=[300,340], yscale='log')
#ps.plot_spectrogram(tspecDC,fspecDC,np.abs(powercDC),vr=[-60,-10],yr=[10,1000],xr=[300,340], yscale='log')

#tr = [630,660]
tr = [100,300]
good = np.where((tdats >= tr[0]) & (tdats <= tr[1]))
goodDC = np.where((tdat12DC >= tr[0]) & (tdat12DC <= tr[1]))


w1 = filt.butter_bandpass_filter(wf1[good], 2, 80, fss, order= 10)
w2 = filt.butter_bandpass_filter(wf2[good], 2, 80, fss, order= 10)
w3 = filt.butter_bandpass_filter(wf3[good], 2, 80, fss, order= 10)
w4 = filt.butter_bandpass_filter(wf4[good], 2, 80, fss, order= 10)
w12DC = filt.butter_bandpass_filter(wf12DC[goodDC], 2, 80, fsDC, order= 10)
w34DC = filt.butter_bandpass_filter(wf34DC[goodDC], 2, 80, fsDC, order= 10)

"""
trz = [651,655]
#trz = [652.5,653] #25 Hz wiggles
fig, axs = plt.subplots(4)
ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[10,200],xr=trz, yscale='linear',ax=axs[0])
ps.plot_spectrogram(tspecs2,fspecs2,np.abs(powercs2),vr=[-60,-40],yr=[10,200],xr=trz, yscale='linear',ax=axs[1])
axs[2].plot(tdats[good],w1)
axs[2].set_xlim(trz)
axs[2].set_ylim(-0.003,0.003)
axs[3].plot(tdats[good],w1-w2)
axs[3].set_xlim(trz)
axs[3].set_ylim(-0.003,0.003)
"""


fig,axs = plt.subplots(2)
axs[0].plot(tdat12DC[goodDC],w12DC)
axs[1].plot(tdat34DC[goodDC],w34DC)
for i in range(6): axs[i].set_xlim(150,151)
for i in range(2): axs[i].set_ylim(-0.001,0.001)
for i in range(2): axs[i+2].set_ylim(-0.0006,0.0006)
axs[4].set_ylim(-0.4,0.4)
axs[5].set_ylim(-0.4,0.4)



fig,axs = plt.subplots(6)
axs[0].plot(tdats[good],w1,tdats[good],w2)
axs[1].plot(tdats[good],w3,tdats[good],w4)
axs[2].plot(tdats[good],w1-w2)
axs[3].plot(tdats[good],w3-w4)
axs[4].plot(tdat12DC[goodDC],w12DC)
axs[5].plot(tdat34DC[goodDC],w34DC)
for i in range(6): axs[i].set_xlim(150,151)
for i in range(2): axs[i].set_ylim(-0.001,0.001)
for i in range(2): axs[i+2].set_ylim(-0.0006,0.0006)
axs[4].set_ylim(-0.4,0.4)
axs[5].set_ylim(-0.4,0.4)

#Density structure events 

#150-170 - lots of few Hz structuring
#166.4 (12)

#212.5 (12-34)  - large
#240-260 (12-34) - lots of 1/3 Hz fluctuations
#260.8
#268.6 (12-34)
#286.4 (34)

#304.75 (12)
#308.1 (12-34)
#312.75 (34)
#338.75 (12)
#376.7 (12-34)
#395.5 (12)
#396.8 (12-34)
#397.3 (12)

#400-600 (intermittent ~quarter to half second turbulence - e.g. 430-434)
#647.2-647.6 (12)
#648.4 (34)
#650.2 (12-34)
#652.7 (12-34)
#692.3 (12-34)
#705.5 (12-34)
#706.5 (12-34)
#722.8 (12,34)




#Timing analysis on density structure (**only good one I've found so far**)
tshift = 0.003
plt.plot(tdats[good],w1,tdats[good]-tshift,-w2)
plt.xlim(652.4,653.0)
plt.ylim(-0.002,0.002)


plt.plot(tdats[good],w3,tdats[good],-w4)
plt.xlim(648.0,648.5)
plt.ylim(-0.001,0.001)








"""
#-----------------------------------
#Artifical signal with harmonics that slowly rise in tone over entire mission
#-----------------------------------

ps.plot_spectrogram(tspecs,fspecs,np.abs(powercs),vr=[-60,-40],yr=[80,200],xr=[380,420], yscale='linear')

tr = [388,395]
good = np.where((tdats >= tr[0]) & (tdats <= tr[1]))


w1 = filt.butter_highpass_filter(wf1[good], 1, fss, order= 10)
w2 = filt.butter_highpass_filter(wf2[good], 1, fss, order= 10)
w3 = filt.butter_highpass_filter(wf3[good], 1, fss, order= 10)

w4 = filt.butter_highpass_filter(wf4[good], 1, fss, order= 10)



plt.plot(tdats[good],w1,tdats[good],w2,tdats[good],w3,tdats[good],w4)
plt.xlim(390.4,390.5)

fig,axs = plt.subplots(2)
axs[0].plot(tdats[good],w1)
axs[1].plot(tdats[good],w2)


plot_csd_skins([390,392])

"""
