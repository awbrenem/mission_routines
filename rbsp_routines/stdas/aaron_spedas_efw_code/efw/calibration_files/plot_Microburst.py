"""Author: Karl Pederson"""

import numpy as np
import math
import matplotlib
from matplotlib import pyplot as plt

class plot_Microburst:
    
    def __init__(self):
        self.__dict__.update()
    
    """Returns f0 value (Amplitude)"""
    def getF0(self, f0):
        return f0
    
    """Returns E0 value"""
    def getE0(self, E0):
        return E0
   
    """Returns set of evals (Will be a numpy array)"""    
    def getevals(self, evals):
        return evals
    
    """Gets the final eval array and the f_0 array and plots them using
    a exp fit"""
    def getexpfit(self, f0, E0, evals):
        #f_0 = f0 * math.exp(-1 * evals/E0)
        newEvals = np.where(evals >=0)
        newEvals = list(newEvals)  
        f_0 = []
        #print(f_0)
        for e1 in newEvals:
            newEvals2 = e1
            for e2 in e1:
                f_0.append(f0 * math.exp(-1 * e2/E0))
        
        plt.axes(xlabel="Energy(keV)", ylabel="Flux", title="Graph of Flux Based on Energy(keV) (Exp Fit)")   
        plt.plot(newEvals2, f_0)
        
    """Gets the final eval array and the f_0 array and plots them using
    a ln fit"""    
    def getlogfit(self, f0, E0, evals):
        #f_0 = f0 * -1 * np.log(-1 * evals/E0) 
        newEvals = np.where(evals >= 0)
        newEvals = list(newEvals)
        f_0 = []
        for e1 in newEvals:
            newEvals2 = e1
            for e2 in e1:
                f_0.append(f0 * -1 * np.log(-1 * e2/E0))
                
        plt.axes(xlabel="Energy(keV)", ylabel="Flux", title="Graph of Flux Based on Energy(keV) (Log Fit)")   
        plt.plot(newEvals2, f_0)      
    
    """Creates a Spectrogram using the Exp fit"""    
    def getexpspecgram(self, f0, E0, evals):
        #f_0 = f0 * math.exp(-1 * evals[0]/E0)
        finalFlux = np.empty((20,1000))
        #print(finalFlux)
        timecount = 0
        for time in finalFlux:
            ecount = 0
            
            for e in time:
                finalFlux[timecount][ecount] = f0 * math.exp(-1 * evals[ecount]/E0)
                ecount+=1
            timecount +=1
        #print(finalFlux)
        
                
        
       
        plt.specgram(finalFlux, NFFT = 441, Fs = 20, noverlap = 128)
        plt.yticks(np.arange(0,10,1))
        plt.xticks(np.arange(0,20,1))
        plt.title("Spectrogram of Flux Intensity over Time from 0-1000keV (Exp Fit)")
        plt.xlabel("Time")
        plt.ylabel("Energy (keV)")
        plt.colorbar(label="Flux")
        plt.show()
        
        
        
        

    

       
"""Test"""
f0 = plot_Microburst()
f0 = f0.getF0(1100)
E0 = plot_Microburst()
E0 = E0.getE0(350)
evals = plot_Microburst()
evals = evals.getevals(np.arange(0,1000,1))
#print(evals)

test_plot = plot_Microburst()
#test_plot.getexpfit(f0, E0, evals)
#test_plot.getlogfit(f0, E0, evals)
test_plot.getexpspecgram(1100, 350, np.arange(0,1000))


