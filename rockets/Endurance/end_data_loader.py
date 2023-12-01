"""
Load Endurance data (for FIELDS see end_fields_loader.py)

def load_slp_potential(t=0)
def load_slp(t=0):
def load_timeline()
def load_iri(alt=0)
def load_ephemeris(t=0)
def load_ig

"""


#--------------------------------------------------------------------------------------------------------
#Load ionization gauge data (Neutral dens/temp vs alt). Note that this is a a much higher cadence than PES data
#--------------------------------------------------------------------------------------------------------

def load_ig(alt=0):

    import pandas as pd
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/IG/Endurance_AxialIG_profiles_v1.txt'
    header = ['Alt','Dens_neutral','Temp_neutral']
    #alt in km
    #neutral dens in cm-3

    df = pd.read_csv(fn, skiprows=3,names=header,delim_whitespace=True)

    #If specific time is requested...
    if alt != 0:
        df = df[df['Alt'] >= alt]
        return df.iloc[0]
    else:
        return df


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
        df = df[df['ToF [s]'] >= t]
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

    #Start and stop times of good science data. 
    goodscienceS = [125, 205, 285, 365, 445, 525, 625, 705, 785, 865]
    goodscienceE = [195, 275, 355, 435, 515, 595, 695, 775, 855, 900.6]


    return df, goodscienceS, goodscienceE






#----------------------------------------
#Load IRI data (single time or all)
#Note that the official IRI run on the Endurance Box website is set up for the location/conditions at the beginning of the mission. 
#--Year= 2022., Month= 05, Day= 11, Hour=1.5,
#--Time_type = Universal
#--Coordinate_type = Geographic
#--Latitude= 78.9235, Longitude= 11.909, Height= 100.
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


#---------------------------------------------------------
#Load PES e- data from 1 eV to 1 keV 
#Returns separate Pandas dataframes for inflow and outflow
#erange --> select specific energy range
#---------------------------------------------------------

def load_PES(erange=[0,1000]):

    import pandas as pd
    import datetime as dt
    import numpy as np 

    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/PES/Endurance_PES_flux_ESA_mode_(standard_resolution)_outflow_vs_inflow_2023-03-13.csv'

    nms =  ["direction", 0.1370,8.4128,9.1497,9.9431,11.5373,13.3872,15.5305,18.0196,20.4224,23.6898,28.1379,32.6274,37.8539,43.9113,50.9463,59.1018,68.5621,79.5376,92.2691,107.0311,123.9898,143.8542,166.9112,193.6813,224.6471,260.6623,302.3216,350.7751,406.9147,472.0391,547.5972,635.2592,736.9189,854.8783,991.6627, 1000]

    

    df = pd.read_csv(fn, skiprows=18, header=0, names=nms)

    df.index = pd.to_datetime(df.index)

    #Modify times so that t=0 corresponds to launch at 01:31 GMT
    t0 = dt.datetime(2022, 5, 11, 1, 32, 54, 338)
    deltat = df.index - t0
    deltat = deltat.seconds 
    df.index = deltat



    #remove white spaces
    df['direction'] = df['direction'].str.strip()

    df_inflow = df[df['direction'] == "inflow"]
    df_outflow = df[df['direction'] == "outflow"]


    #sum over requested energy range 
    dfkeys = df_inflow.keys()
    dfkeys = dfkeys[1:]
    good = np.where((dfkeys >= erange[0]) & (dfkeys <= erange[1]))

    df_inflow = df_inflow.iloc[:,good[0]+1]
    df_outflow = df_outflow.iloc[:,good[0]+1]


    return df_inflow, df_outflow


