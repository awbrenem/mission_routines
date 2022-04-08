"""
Functions for plotting VISIONS-2 data

PlotParticle --> helps to plot data 
IntegrateSpectra --> Sum over energy or pitch-angles 
"""


import math 
import scipy.signal 
import numpy as np 

# Function for plotting particle data

def PlotParticle(xvals,yvals,zvals,p,axs,**PlotParams):

    pcm1 = axs[p].pcolormesh(xvals,yvals,np.transpose(np.log10(zvals)),shading='auto',cmap='turbo',vmin=PlotParams["vmin"],vmax=PlotParams["vmax"])
    axs[p].set_title(PlotParams["title"])
    axs[p].set_yscale(PlotParams["yscale"])
    axs[p].set_ylim(PlotParams["ylim"])
    axs[p].set_ylabel(PlotParams["ylabel"])
    axs[p].set_xlabel(PlotParams["xlabel"])
    if PlotParams["colorbar"] == 1:
        fig.colorbar(pcm1, ax=axs[p])



"""
Integrate generic spectra over certain y-range for each time value
Vals -> Spectra (size [xvals, yvals])
ylow, yhig -> min and max yvals to integrate over 
smoothtime -> smooth the final results by this time (sec)
norm -> set to normalize values


Examples:
smootime = 10. #smooth time in sec

vlfAmp_smoothed = IntegrateSpectra(np.transpose(Sxx), spectimes, specfreqs, 3000, 10000, smootime)
eAmpPerp_smoothed = IntegrateSpectra(elecPerp["flux"],elecPerp["times"],elecPerp["energies"],3,3000,smootime)
eAmpPar_smoothed = IntegrateSpectra(elecDowngoing["flux"],elecDowngoing["times"],elecDowngoing["energies"],3,3000,smootime)
iAmpPerp_smoothed = IntegrateSpectra(ionsPerp["flux"],ionsPerp["times"],ionsPerp["energies"],3,3000,smootime)
iAmpPar_smoothed = IntegrateSpectra(ionsDowngoing["flux"],ionsDowngoing["times"],ionsDowngoing["energies"],3,3000,smootime)
eAmpLowE_smoothed = IntegrateSpectra(elecLowE["flux"],elecLowE["times"],elecLowE["pitchangles"],-180,180,smootime)
eAmpHigE_smoothed = IntegrateSpectra(elecHigE["flux"],elecHigE["times"],elecHigE["pitchangles"],-180,180,smootime)
iAmpLowE_smoothed = IntegrateSpectra(ionsLowE["flux"],ionsLowE["times"],ionsLowE["pitchangles"],-180,180,smootime)
iAmpHigE_smoothed = IntegrateSpectra(ionsHigE["flux"],ionsHigE["times"],ionsHigE["pitchangles"],-180,180,smootime)


"""
def IntegrateSpectra(vals, xvals, yvals, ylow, yhig, smoothtime, norm=0):

    goodf = np.where((yvals >= ylow) & (yvals <= yhig))
    finVals = []
    for i in range(len(xvals)):
        finVals.append(np.sum(vals[i, goodf])/np.size(goodf))

    sr = 1/(np.median(xvals - np.roll(xvals,1)))
    npts = math.floor(sr * smoothtime)
    #require npts to be an odd number
    if npts % 2 == 0:
        npts = npts + 1

    tmp = scipy.signal.savgol_filter(finVals, npts, 5)
    if norm == 1:
        tmp = tmp/np.max(tmp)

    return tmp
