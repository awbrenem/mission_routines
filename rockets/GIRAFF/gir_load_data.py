"""
Load GIRAFF data (for FIELDS data that need to be calibrated see gir_load_fields.py)
def load_EDCdespun_Henry
def load_slp
def load_zesta_mag
def load_fields_magNWU_henry
def load_ephemeris
def load_marilia_HF 
def load_igrf
def load_flh_byeye 
def load_fuh_byeye
def load_fpe_byeye


***NOT yet implemented***
def load_slp_sonogram()
def load_slp_potential(t=0)
def load_timeline()
def load_ig

"""



"""
Zonal and Meridional E-field data from Henry 
"""
def load_EDCdespun_Henry(pld):
    from scipy.io import readsav 
    import numpy as np
    
    path = "/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/efield_DC/"

    if pld == "380":
        fn = "36380_FULL_Esoln_corrected_P16.sav"
    if pld == "381":
        fn = "36381_Esoln_FULL.sav"

    vals = readsav(path + fn)

    #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
    good = np.squeeze(np.where(vals.etime >= 0.))

    d = {'Emerid':vals.emer[good],
         'Ezon':vals.ezon[good],
         'time':vals.etime[good]}

    return d




#Load lower hybrid values by-eye identified from gir_cal_determine_accurate_lower_hybrid_freq.py
def load_flh_byeye(pld):
    import pickle
    flhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/freqs_by-eye/GIRAFF_'+pld+'_lower_hybrid_freqs_byeye.pkl'
    v = pickle.load(open(flhfile,'rb'))
    return v[0][:,0], v[0][:,1]  #times, flh values

#Load upper hybrid values by-eye identified from gir_cal_determine_density_HF.py
def load_fuh_byeye(pld):
    import pickle
    fuhfile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/freqs_by-eye/GIRAFF_'+pld+'_fuh_vals_byeye.pkl'
    v = pickle.load(open(fuhfile,'rb'))
    return v[0], v[1]  #times, fuh values

#Load plasma frequency values by-eye identified from gir_cal_determine_density_HF.py
def load_fpe_byeye(pld):
    import pickle
    fpefile = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/freqs_by-eye/GIRAFF_'+pld+'_fpe_vals_byeye.pkl'
    v = pickle.load(open(fpefile,'rb'))
    return v[0], v[1]  #times, fpe values



"""
Load IGRF data 
NOTE: to get this data I had to start Python with a bash shell (not CSCode) and run pyIGRF to generate the pickle file.
    import pyIGRF
    import sys 
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
    import gir_load_data
    import numpy as np
    import pickle
    
    glat = 67
    glon = -147
    pld = '381'
    ephem = gir_load_data.load_ephemeris(pld)

    alt = ephem['altitude']  #in km
    times = ephem['time']

    BoIGRF = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2024)[6] for i in alt])

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/IGRF/'
    fn = 'giraff_'+pld+'_igrf.pkl'
    pickle.dump([times, BoIGRF], open(path+fn,'wb'))


"""
def load_igrf(pld):
    import numpy as np
    import pickle
    #import matplotlib.pyplot as plt

    fn = 'giraff_'+pld+'_igrf.pkl'
    path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/IGRF/'
    vals = pickle.load(open(path+fn,'rb'))

    alts_igrf = vals[0]
    Bo = vals[1]


    #convert alts to mission times
    ephem = load_ephemeris(pld)
    mission_times = ephem['time']
    mission_altitudes = ephem['altitude']

    for i in range(len(mission_times)):
        alt = mission_altitudes[i]
        #Find nearest altitude in IRI data
        idx = (np.abs(alts_igrf - alt)).argmin()

        if i == 0:
            alts_interp = [alts_igrf[idx]]
        else:
            alts_interp = np.append(alts_interp, alts_igrf[idx])

    #plt.plot(alts_interp,Bo)
    #plt.plot(mission_times,Bo)

    return mission_times, Bo






"""
Loads the Marilia HF spectral data pickles that I created in gir_create_marilia_HF_pickle_files.py 
Returns times, freqs, and power spectral density data.
"""
def load_marilia_HF(pld):

    import pickle

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/marilia_HF/'
    file = path + 'GIRAFF_'+pld+'_HF_spectra_timecorrected_reducedFreqRange.pkl'
    vals = pickle.load(open(file,'rb'))


    return vals[0], vals[1], vals[2]    




#Load the EFIELDS mag data rotated into NWU coordinates by Henry
def load_fields_magNWU_henry(pld):
    from scipy.io import readsav 
    import numpy as np

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/mag/'
    if pld == '380':
        fn = '36380_New_Calibrated_B_xyz_nwu_16p.sav' #From Henry
    if pld == '381':
        fn = ''

    vals = readsav(path + fn)

    #Remove negative times (starts at t=-100 sec). Not doing so messes up my spectrogram plotting routines.
    good = np.squeeze(np.where(vals.btime >= 0.))

    d = {'BNorth':vals.bn[good],
         'BWest':vals.bw[good],
         'BUp':vals.bu[good],
         'time':vals.btime[good]}

    return d




#Load the despun data from Zesta mag provided by Alex Hoffman.
#NOTE: I add a 2.67 sec timeshift to align this data with the Fields data (determined from gir_analysis_FA_currents.py). 
#"The 381 data after about 400 seconds is not real, and there are a few spikes in each that we think was caused by radiation hitting the ADC."
#There were two Zesta mags on each payload, and inner and outer. I'm only returning data from the outer as it's much better quality.
def load_zesta_mag(pld):

    import pandas as pd

    if pld == '380':
        fn = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/mag/380_zesta_ENU_magnetic_field_measurements.csv'
    if pld == '381':
        fn = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/mag/381_zesta_ENU_magnetic_field_measurements.csv'


    df = pd.read_csv(fn)
    #print(df.head())

    dfn = df.to_numpy()

    timeshift = 2.486 #Need to timeshift data to bring in align with Fields data.
    #timeshift = 2.67 #Need to timeshift data to bring in align with Fields data.

    vals = {'time':dfn[:,0]+timeshift,
            'latgeo':dfn[:,1],
            'longeo':dfn[:,2],
            'alt':dfn[:,3]/1000,
            'BEast':dfn[:,4],
            'BNorth':dfn[:,5],
            'BUp':dfn[:,6]
            }

    return vals



#----------------------------------------
#Load Langmuir probe data from E-fields package
#Values in Amps
#----------------------------------------

def load_slp(pld, t=0):

    from scipy.io import readsav 

    path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/langmuir/'
    if pld == '380':
        fn = '36380_TM1_LP_S9_amps.sav'
    if pld == '381':
        fn = '36381_TM1_LP_S9_amps.sav'

    dat = readsav(path+fn)

    #dict_keys(['tlp', 'dlp_amps', 'dlp_raw', 'dlp_rawnormal', 'dlp18bit', 'alp', 'blp', 'author', 'timenote', 
    #       'calnote', 'dataunits', 'flightno', 'format', 'link', 'samplerate', 't0', 'timetagmethod', 'timeunits'])

    times = dat['tlp']
    vals = dat['dlp_amps']



    #If specific time is requested...
    if t != 0:
        return times[times >= t][0], vals[times >= t][0]
    else:
        return times, vals




#----------------------------------------
#Load ephemeris data (single time or all)
#----------------------------------------

def load_ephemeris(pld, t=0):

    from scipy.io import readsav
    import pandas as pd


    if pld == '380':
        fn = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/ephemeris/36380_GPS2POSDAT_Main_wheader.sav'
    if pld == '381':
        fn = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/ephemeris/36381_GPS2POSDAT_Main_wheader.sav'

    #dat = readsav(fn)
    #dict_keys(['time', 'altitude', 'latitude', 'longitude', 'ewvel', 'nsvel', 'udvel', 'flight', 'location', 't0_year', 't0_day', 't0_hour', 't0_min', 't0_sec'])




        # Read the .sav file using scipy.io.readsav
    df = readsav(fn)
    df['altitude'] = df['altitude'] / 1000.0  # Convert to km

    
    #If specific time is requested...
    if t != 0:
        df = df[df['Time'] >= t]
        return df.iloc[0]
    else:
        # Convert df to a numpy dictionary
        # df = {col: df[col].to_numpy() for col in df.columns}
        return df
    