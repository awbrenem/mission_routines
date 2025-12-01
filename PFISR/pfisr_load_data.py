
"""
Load PFISR .h5 data from MADRIGAL website. 
https://data.amisr.com/database/


NOTE: see example at bottom for loading this for GIRAFF data and plotting vs mission time. 



def pfisr_load_data(filename, type, plot_test=0):
    #Load full spectrogram data [ntimes, nbeams, naltitudes]

def pfisr_slice_data(times, alt, vals, tslice, beamnumber, plot_test=0):
    #Extract PFISR values near a particular time and for a particular beam.

def pfisr_vals_to_mission_times(pfisrdat, pfisralts, mission_times, mission_altitudes, test=0):
    #Convert pfisr values vs altitude to values vs time for a particular mission


***NOTE: pfisr_load_data currently only works for density. Not yet sure how to get the altitude values for the temp and composition measurements ("MSIS")

###NOTE: see pfisr_vals_to_mission_times to convert altitude values to mission times. 


ac - alternating code - "employed to improve range resolution for detailed observations of specific regions like the D- and E-regions of the ionosphere."
lp - long pulse code - "Used for measurements requiring high signal-to-noise ratios, often in F-region studies" 
vvelsLat - resolved velocity mode

Density products [ntimes, nbeams, naltitudes]:
    Ne_mod:  is this regular density?? 
    Ne_NoTr: "processing mode that is power based, as opposed to spectra fitting, 
            electron densities with no correction for ion temperatures, i.e. Te/Ti=1"

        
#Temperature and composition:
    Tn: 
    nAr: 
    nH:
    nHe:
    nMass:
    nN:
    nN2:
    nNO:
    nO:
    nO2:

#GEOMAG values         
    maglat: 
    maglon: 


Example files
20250202.002_ac_5min-fitcal.h5
20250202.002_ac_20min-fitcal.h5
20250202.002_lp_1min-fitcal.h5
20250202.002_lp_5min-fitcal-vvelsLat-300sec.h5
20250202.002_lp_5min-fitcal.h5
20250202.002_lp_20min-fitcal.h5

"""
import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import h5py
from datetime import datetime, timezone
import numpy as np


#-----------------------------------------------------------
#Load full spectrogram data [ntimes, nbeams, naltitudes]
#-----------------------------------------------------------

def pfisr_load_data(filename, type, plot_test=0):

    #Determine which h5 group to load based on type
    if type == 'Ne_Mod' or type == 'Ne_NoTr':
            groupname = 'NeFromPower'
    if type == 'Tn' or type == 'nAr' or type == 'nH' or type == 'nHe' or type == 'nMass' or type == 'nN' or type == 'nN2' or type == 'nNO' or type == 'nO' or type == 'nO2':
            groupname = 'MSIS'

    #Read in the h5 file
    with h5py.File(filename, 'r') as f:
        print("Keys: %s" % f.keys())
        base_items = list(f.items())
        print(base_items)


        #-------------
        #Time values
        #-------------

        G1 = f.get('Time')
        #G1_items = list(G1.items())

        unix = np.array(G1.get('UnixTime')) #shape of [ntimes, 2] - start and end times in unix time

        times_utc = np.empty(unix[:,0].shape, dtype=object)
        times_local = times_utc.copy()

        for i in range(len(times_utc)):
            times_local[i] = datetime.fromtimestamp(unix[i][0])
            times_utc[i] = datetime.fromtimestamp(unix[i][0], timezone.utc)


        #---------------
        # Load desired data
        #---------------
            
        G1 = f.get(groupname)
        #G1_items = list(G1.items())

        alt = np.array(G1.get('Altitude'))/1000  #km. [nbeams, naltitudes]. (includes some negative altitudes, like seen in online spectra)
        vals = np.array(G1.get(type))



        """

        #--------------
        #Temperature and composition
        #--------------

        G1 = f.get('MSIS')
        G1_items = list(G1.items())
        G1 = f.get('NeFromPower')

        Tn = np.array(G1.get('Tn')) 
        nAr = np.array(G1.get('nAr')) 
        nH = np.array(G1.get('nH')) 
        nHe = np.array(G1.get('nHe')) 
        nMass = np.array(G1.get('nMass')) 
        nN = np.array(G1.get('nN')) 
        nN2  = np.array(G1.get('nN2')) 
        nNO = np.array(G1.get('nNO')) 
        nO = np.array(G1.get('nO')) 
        nO2 = np.array(G1.get('nO2'))    


        #GEOMAG values         
        G1 = f.get('Geomag')
        G1_items = list(G1.items())
        maglat = np.array(G1.get('MagneticLatitude'))
        maglon = np.array(G1.get('MagneticLongitude'))
        """    


        #----------------------------
        #Plot a spectrogram 
        if plot_test:
            import plot_spectrogram as ps

            tlaunch = datetime(2025,2,2,7,7,00, tzinfo=timezone.utc)  #GIRAFF 36.381 launch time

            tdelta = np.zeros(len(times_local), dtype='timedelta64[s]')
            for i in range(len(times_local)):
                tdelta[i] = times_utc[i] - tlaunch


            #extract example beam
            beamno = 12
            valsB = vals[:,beamno-1,:]
            alt12 = alt[beamno-1,:]

            timesplot = tdelta.astype(int)
            ps.plot_spectrogram(timesplot, alt12, valsB.T, vr=[1e10,1e12], yr=[0,400], 
                                zscale='linear',yscale='linear', title='PFISR Data for single beam', 
                                xlabel='Seconds', ylabel='Altitude (km)', plot_kwargs={'cmap':'turbo'})


    return times_utc, alt, vals


#---------------------------------------------------------------------
#Extract PFISR values near a particular time and for a particular beam.
#---------------------------------------------------------------------

def pfisr_slice_data(times, alt, vals, tslice, beamnumber, plot_test=0):


    #Determine closest PFISR time to desired time
    tdelta = np.zeros(len(times), dtype='timedelta64[s]')
    for i in range(len(times)):
        tdelta[i] = times[i] - tslice

    mn = min(tdelta, key=abs)
    tindex = np.flatnonzero(tdelta == mn)[0]

    valsBN = vals[tindex,beamnumber-1,:] #[ntimes, naltitudes]
    altBN = alt[beamnumber-1,:]


    if plot_test:
        import matplotlib.pyplot as plt

        plt.plot(altBN, valsBN, '.-')
        plt.show()


    return altBN, valsBN




#---------------------------------------------------------------------
#Convert PFISR model values (a function of altitude) to mission times
#---------------------------------------------------------------------

#pfisrdat - data from the PFISR model read in from pfisr_load_data.py
#mission_times - array of mission times (s)
#mission_altitudes - array of mission altitudes (km)


def pfisr_vals_to_mission_times(pfisrdat, pfisralts, mission_times, mission_altitudes, plot_test=0):
        
    import numpy as np 

    #For each time of ephem altitude data find the value of fpe from IRI at that altitude
    for i in range(len(mission_times)):
        alt = mission_altitudes[i]
        #Find nearest altitude in IRI data
        idx = (np.abs(pfisralts - alt)).argmin()
 
        if i == 0:
            pfisrdat_interp = [pfisrdat[idx]]
            pfisralts_interp = [pfisralts[idx]]
        else:
            pfisrdat_interp = np.append(pfisrdat_interp, pfisrdat[idx])
            pfisralts_interp = np.append(pfisralts_interp, pfisralts[idx])


    if plot_test:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(2, 1, figsize=(10, 8))

        axs[0].plot(mission_times, pfisralts_interp, label='alt PFISR', color='blue', linestyle='--')
        axs[0].plot(mission_times, mission_altitudes, label='alt mission', color='red', linestyle='-')
 
        axs[1].plot(mission_times, pfisrdat_interp, label='PFISR data', color='blue', linestyle='--')
        for ax in axs:
            ax.legend()
        axs[0].set_ylabel('Altitude (km)')
        axs[1].set_ylabel('PFISR Data')
        plt.xlabel('Time (s)')
        plt.ylabel('PFISR data')
        plt.title('PFISR Values at Mission Times')
        plt.legend()
        plt.show()  

    return {
        'times': np.array(mission_times),
        'pfisralts_interp': np.array(pfisralts_interp),
        'pfisrdat_interp': np.array(pfisrdat_interp)
    }




"""

#-----------------------------------------------------------------------------------------------------
#Example - load GIRAFF 36.380 density and plot fpe vs altitude and vs time for the GIRAFF flight 36.380
#----------------------------------

filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/pfisr/36.381/20250202.002_ac_1min-fitcal.h5'
type = 'Ne_Mod'  #Density in m-3
pfisrTimes, pfisrAlts, pfisrNeMod = pfisr_load_data(filename, type, plot_test=1)


tslice = datetime(2025,2,2,7,7,00, tzinfo=timezone.utc)  #GIRAFF 36.381 launch time
beamnumber = 12
pfisrAltsS, pfisrNeModS = pfisr_slice_data(pfisrTimes, pfisrAlts, pfisrNeMod, tslice, beamnumber, plot_test=0)


#load altitude data for GIRAFF 36.380
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
from gir_load_data import load_ephemeris
import matplotlib.pyplot as plt 
ephem  = load_ephemeris('380')
plt.plot(ephem['time'], ephem['altitude'])

pvst = pfisr_vals_to_mission_times(pfisrNeModS, pfisrAltsS, ephem['time'], ephem['altitude'], test=0)


#fce (Hz) vs time
plt.plot(pvst['times'],8980*np.sqrt(pvst['pfisrdat_interp']/1e6))

#compare original time vs alt with PFISR interpolated time vs altitude
plt.plot(ephem['time'],ephem['altitude'])
plt.plot(pvst['times'],pvst['pfisralts_interp'])


print('h')

"""


