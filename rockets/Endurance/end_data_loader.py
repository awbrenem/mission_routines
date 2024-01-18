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



def load_slp_fixedbias():

    import pandas as pd 
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_FixedBiasDensity_Downleg_LowAlt.txt'
    nms = ['ToF [s]','UTC Time','Altitude [km]','Normalized Density [\m3]','ACS Binary']

    df = pd.read_csv(fn,skiprows=0,header=0,names=nms,sep='\t')

    """
    #If specific time is requested...
    if t != 0:
        df = df[df['ToF [s]'] >= t]
        return df.iloc[0]
    else: 
        return df
    """
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

    badscienceS = [0] + goodscienceE
    badscienceE = goodscienceS + [1000]

    return df, goodscienceS, goodscienceE, badscienceS, badscienceE






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
    import numpy as np 
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/models/iri/IRI_run_realistic_sunspot_F10_7.txt'
    nms = ['Height(km)','Ne(m-3)','Ti(K)','Te(K)','O_ions','H_ions','He_ions','O2_ions','NO_ions','N_ions']

    df = pd.read_csv(fn,skiprows=40,header=0,names=nms,delim_whitespace=True)

    #-----------------------------------------------
    #associate upleg/downleg times based on Endurance data

    ephem, ephem2 = load_ephemeris()
    alts_end = ephem['Altitude']
    alts_iri = df['Height(km)']
    times_end = ephem['Time']
    times_iri_upleg = np.zeros(len(df['Height(km)']))
    times_iri_downleg = np.zeros(len(df['Height(km)']))


    for i in range(len(alts_iri)):
        goo = np.where(alts_end >= alts_iri[i])
        if len(goo[0]) != 0:
            times_iri_upleg[i] = times_end[goo[0][0]]

        goo = np.where(alts_end <= alts_iri[i])
        if len(goo[0]) != 0:
            maxloc = np.argmax(alts_end[goo[0]])
            times_iri_downleg[i] = times_end[goo[0][maxloc]]

    #import matplotlib.pyplot as plt
    #plt.plot(times_end, alts_end,'.')
    #plt.plot(times_iri_upleg, alts_iri,'.')
    #plt.plot(times_iri_downleg, alts_iri,'.')
    #-----------------------------------------------


    d = {'times_upleg':times_iri_upleg, 'times_downleg':times_iri_downleg}
    dfu = pd.DataFrame(data=d)
    df = pd.concat([dfu,df],axis=1,join='inner')
    #plt.plot(df2['times_upleg'],df2['Height(km)'])


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

    nms1 = ['Time','ECEF-X-VEL','ECEF-Y-VEL','ECEF-Z-VEL','Lat','Long','Altitude']
    df = pd.read_csv(fn,names=nms1,skiprows=1)


    fn2 = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/ephemeris/47.001_SynchronizedSolution.csv'
    nms = ['Time','Yaw','Pitch','Roll','RollRate','AoA_T','AngleB','a11','a12','a13','a21','a22','a23','a31','a32','a33','X_Az','X_El','Y_Az','Y_El','Z_Az','Z_El','Latgd','Long','Alt']
    df2 = pd.read_csv(fn2,names=nms,skiprows=7)



    #If specific time is requested...
    if t != 0:
        df = df[df['Time'] >= t]
        df2 = df2[df2['Time'] >= t]

        return df.iloc[0], df2.iloc[0]
    else:
        return df, df2


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


