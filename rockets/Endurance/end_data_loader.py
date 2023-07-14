"""
Load Endurance data (for FIELDS see end_fields_loader.py)

def load_slp_potential(t=0)
def load_slp(t=0):
def load_timeline()
def load_iri(alt=0)
def load_ephemeris(t=0)


"""


#----------------------------------------
#Load swept Langmuir probe SC potential 
#----------------------------------------

def load_slp_potential(t=0):

    import pandas as pd
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_Floating_Potential_02282023.csv'
    df = pd.read_csv(fn)

    #If specific time is requested...
    if t != 0:
        df = df[df['Flight Time'] >= t]
        return df.iloc[0]
    else:
        return df



#----------------------------------------
#Load swept Langmuir probe data
#----------------------------------------

def load_slp(t=0):

    import pandas as pd
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_05_07_2023_witherrorbars.csv'
    df = pd.read_csv(fn)

    #If specific time is requested...
    if t != 0:
        df = df[df['Flight Time'] >= t]
        return df.iloc[0]
    else:
        return df


#----------------------------------------
#Load timeline data of events 
#----------------------------------------

def load_timeline():
    import pandas as pd

    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/timeline/47001_Timeline.csv'
    df = pd.read_csv(fn)
    return df



#----------------------------------------
#Load IRI data (single time or all)
#----------------------------------------

def load_iri(alt=0):
    import pandas as pd

    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/models/iri/IRI_run_realistic_sunspot_F10_7.txt'
    nms = ['Height(km)','Ne(m-3)','Ti(K)','Te(K)','O_ions','H_ions','He_ions','O2_ions','NO_ions','N_ions']

    df = pd.read_csv(fn,skiprows=40,header=0,names=nms,delim_whitespace=True)

    #If specific alt is requested...
    if alt != 0:
        df = df[df['Height(km)'] >= alt]
        return df.iloc[0]
    else: 
        return df


#----------------------------------------
#Load ephemeris data (single time or all)
#----------------------------------------

def load_ephemeris(t=0):

    import pandas as pd
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/ephemeris/Endurance_GPS_velocity_position_altitude.csv'

    df = pd.read_csv(fn)

    #If specific time is requested...
    if t != 0:
        df = df[df['Flight Time'] >= t]
        return df.iloc[0]
    else:
        return df






