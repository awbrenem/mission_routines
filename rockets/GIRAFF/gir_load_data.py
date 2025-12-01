"""
Load GIRAFF data (for FIELDS see gir_load_fields.py)
def load_slp(t=0):




***NOT yet implemented***
def load_slp_sonogram()
def load_slp_potential(t=0)
def load_timeline()
def load_iri(alt=0)
def load_ephemeris(t=0)
def load_ig
def load_igrf

"""



#----------------------------------------
#Load swept Langmuir probe density sonogram

#First column is Time of Flight in Seconds.
#All the other columns are Power-spectral-density at specific frequency
#points. The frequency bins are written in the very first row.
#The PSD units are log10(Density^2/Hz) (Aroh email 7/19/24).

#----------------------------------------




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
    