"""
Load Endurance data (for FIELDS see end_fields_loader.py)

def load_slp_sonogram()
def load_slp_potential(t=0)
def load_slp(t=0):
def load_timeline()
def load_iri(alt=0)
def load_ephemeris(t=0)
def load_ig
def load_igrf

"""

#-----------------------------------------------------------------------
#Load the downgoing electrons file Glyn made for me to compare them to VLF modulations

#From Glyn's email on Aug 22, 2024:
#Please find enclosed some curated PES data for each of the 78 scans made. 
# C1: timestamp for the *beginning* of each PES scan. 
# C2: sum of the differential energy flux of electron precipitation. 
# C3: Then the same for electrons less than 55eV, which will capture all the photoelectron precipitation. 
# C4: total influx above the reflection potential, which perfectly cleaves polar rain 
#   from reflected photoelectrons, but is only for some of the scans.


# NOTE: To calculate the rough polar rain influx, just subtract C3 from C2 

#Timestamp (Start of scan)	Total Electron Influx	Influx below 55eV	Total influx above reflection potential
 
#-----------------------------------------------------------------------

def load_pes_electron_precipitation():
    import pandas as pd
    import datetime 
    import numpy as np
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/PES/ElectronprecipitationTotalDEF.xlsx'
    #Note that Glyn mixed up the flux < 55 and total flux columns. It's corrected below.
    header = ['time','el_influx_below_55eV','el_influx_total','total_influx_above_reflection_potential']

    df = pd.read_excel(fn, skiprows=1, names=header)


    #turn bad Nans into proper Nans
    bad = np.where(df['total_influx_above_reflection_potential'] == '             NaN')[0]
    for i in range(len(bad)):
        df['total_influx_above_reflection_potential'][bad[i]] = np.nan


    df['polar_rain'] = df['el_influx_total'] - df['el_influx_below_55eV']
    df_time = df['time']

    ttst = pd.to_datetime(df_time)

    tlaunch = datetime.datetime(2022,5,11,1,31) 

    tdelta = ttst - tlaunch

    tsec = [np.timedelta64(i, 's') for i in tdelta]


    #Remove data with low confidence 
    #after about 780 sec 
    #before about 180 sec 

    df['el_influx_total'][tsec < np.timedelta64(180)] = np.nan
    #df['el_influx_total'][tsec > np.timedelta64(780)] = np.nan
    df['el_influx_below_55eV'][tsec < np.timedelta64(180)] = np.nan
        #the influx data after 780 sec is good, according to Glyn
        #df['el_influx_below_55eV'][tsec > np.timedelta64(780)] = np.nan

    df['polar_rain'][tsec < np.timedelta64(180)] = np.nan
    #df['polar_rain'][tsec > np.timedelta64(780)] = np.nan

    df['total_influx_above_reflection_potential'][tsec < np.timedelta64(180)] = np.nan
    #df['total_influx_above_reflection_potential'][tsec > np.timedelta64(780)] = np.nan


    #df['tSinceLaunch'] = tsec
    #import matplotlib.pyplot as plt
    #plt.plot(tsec, df['total_influx_above_reflection_potential'])
    #plt.plot(df['tSinceLaunch'], df['polar_rain'])

    return df, tsec

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
#Load swept Langmuir probe density sonogram

#First column is Time of Flight in Seconds.
#All the other columns are Power-spectral-density at specific frequency
#points. The frequency bins are written in the very first row.
#The PSD units are log10(Density^2/Hz) (Aroh email 7/19/24).

#----------------------------------------

def load_slp_sonogram():


    import numpy as np
    import pandas as pd

    #Alternative to full file. Rachel Conway divided up the times, freqs, and data into separate files. 
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_PSD_frequency.txt'
    df = pd.read_csv(fn,sep='\r',skiprows=1,header=None)
    freqs = df.to_numpy()

    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_PSD_time.txt'
    df = pd.read_csv(fn,sep='\r',skiprows=1,header=None)
    times = df.to_numpy()

    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_PSD_data.txt'
    df = pd.read_csv(fn,sep='\t',header=None)
    data = df.to_numpy()


    #Load the full SLP density file. I'm not getting the times to extract properly, so I add them manually. 
    #fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_PSD.txt'
    #df = pd.read_csv(fn,sep='\t',index_col=0,skiprows=0,header=0)
    #times = [99.2635537,106.330411,111.214201,116.330558,121.214358,126.330655,131.214505,136.330751,141.214561,146.330848,151.214658,156.330944,161.214744,166.331031,171.214831,176.331118,181.214917,186.331194,191.215034,196.331281,201.215111,206.331086,211.215247,216.331162,221.215324,226.331239,231.215360,236.331305,241.215457,246.331372,251.215523,256.331438,261.215560,266.331464,271.215616,276.331551,281.215743,286.331637,291.215759,296.331664,301.215856,306.331740,311.215902,316.331827,321.215948,326.331873,331.216025,336.331919,341.216061,346.331956,351.216097,356.332022,361.216164,366.332089,371.216240,376.332155,381.216267,386.332191,391.216303,396.332218,401.216319,406.332204,411.216345,416.332220,421.216372,426.332236,431.216116,436.332353,441.216102,446.332369,451.216118,456.332385,461.216124,466.332361,471.216161,476.332397,481.216167,486.332444,491.216203,496.332440,501.216169,506.332446,511.216165,516.332442,521.216191,526.332458,531.216227,536.332454,541.216213,546.332410,551.216199,556.332396,561.216185,566.332382,571.216202,576.332096,581.216208,586.332072,591.216194,596.332048,601.216170,606.332024,611.216176,616.332030,621.216152,626.331996,631.216118,634.149013,641.216114,646.331968,651.216110,656.331994,661.216116,666.331950,671.216102,676.331986,681.215766,686.332013,691.215742,696.331667,701.215768,706.331653,711.215784,716.331669,721.215740,726.331675,731.215485,736.331691,741.215461,746.331667,751.215406,756.331653,761.215423,766.331619,771.215358,776.331233,781.215334,786.331209,791.215300,796.331144,801.215276,806.331150,811.215222,816.331096,821.215157,826.331032,831.215103,836.330978,841.215039,846.330913,891.214888]
    #goo = df.head()
    #ftmp = goo.keys()
    #freqs = [float(i[:-2]) for i in ftmp]
    #return np.asarray(times),np.asarray(freqs),df.to_numpy()


    #get rid of last time. It's sampled at much lower cadence. 
    times = times[0:-1]
    data = data[0:-1,:]


    #Times from the file are off by ~22 sec. See emails with SLP team around May 9, 2025. 
    times = times + 22

    #Because the freq cadence changes at element ~2500, the SLP sonograms aren't plotting correctly. 
    #To avoid this problem I'm dividing them up into low freq and high freq versions
    freqslow = np.squeeze(freqs[0:2500])
    freqshig = np.squeeze(freqs[2500:])
    datalow = data[:,0:2500]
    datahig = data[:,2500:]


    #There's also a time cadence issue on first and last data point that's messing up the timing of the sonograms.
    #times = times[1:-1]
    #datalow = datalow[1:-1,:]
    #datahig = datahig[1:-1,:]


    #develop new time base with even cadence
    #step = times[10] - times[9]
    #tnew = np.arange(times[0],times[-1], np.floor(step))


    #remove density data above 7.6 kHz (difficult to separate real from alised power) and after 860 sec (due to reentry)
    #(from discussions with Aroh)
    bad = np.where(freqshig >= 7500)[0]
    #bad = np.where(freqshig >= 9000)[0]
    datahig[:,bad] = np.nan #-100

    #bad2 = np.where(times >= 830.)[0]
    bad2 = np.where(times >= 850.)[0]
    datahig[bad2,:] = np.nan #-100


    return times,freqslow,datalow,freqshig,datahig

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
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/SLP/Endurance_SLP_FixedBiasDensity_Downleg_LowAlt_v2.txt'
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

    #Start and stop times of good science data (from 47001_Timeline.csv). 
    #These aren't quite correct based on observed broadband power
    #goodscienceS = [125, 205, 285, 365, 445, 525, 625, 705, 785, 865]
    #goodscienceE = [195, 275, 355, 435, 515, 595, 695, 775, 855, 900.6]
 
    #Modified good science start/stop times based on broadband power 
    goodscienceS = [125, 203,   283,   362.9, 442.8, 522.5, 620.5, 702.5, 782.8, 862.5]
    goodscienceE = [192, 272.5, 352.5, 432.5, 512.5, 592.5, 692.5, 772.5, 852.5, 900.6]

    badscienceS = [0] + goodscienceE
    badscienceE = goodscienceS + [1000]

    return df, goodscienceS, goodscienceE, badscienceS, badscienceE



"""
----------------------------------------
#Load IGRF data (single time or all) from a CSV file from 
#https://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi

NOTE: MUCH BETTER TO USE pyIGRF


import pyIGRF
glat = 78 + (55/60)
glon = 11 + (55/60)

BoIGRF = np.asarray([pyIGRF.igrf_value(glat, glon, i, 2022)[6] for i in plot_alt])

----------------------------------------
"""

"""
def load_igrf(alt=0):
    import pandas as pd
    import numpy as np 
    fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/mag/igrf_31188027945.lst.txt'
    nms = ['Height(km)','Bo','flag']
    df = pd.read_csv(fn,skiprows=0,header=0,names=nms,delim_whitespace=True)


    #-----------------------------------------------
    #associate upleg/downleg times based on Endurance data

    ephem, ephem2 = load_ephemeris()
    alts_end = ephem['Altitude']

    #array location of apogee on Endurance
    maxalt_loc = np.where(alts_end == np.max(alts_end))[0][0]


    alts_igrf = df['Height(km)']
    times_end = ephem['Time']
    times_igrf_upleg = np.zeros(len(df['Height(km)']))
    times_igrf_downleg = np.zeros(len(df['Height(km)']))

    au = np.asarray(alts_end[0:maxalt_loc])
    ad = np.asarray(alts_end[maxalt_loc:])
    tu = np.asarray(times_end[0:maxalt_loc])
    td = np.asarray(times_end[maxalt_loc:])

    for i in range(len(alts_igrf)):
        goo = np.where(au >= alts_igrf[i])
        if len(goo[0]) != 0:
            times_igrf_upleg[i] = tu[goo[0][0]]

        boo = np.where(ad <= alts_igrf[i])
        if len(boo[0]) != 0:
            maxloc = np.argmax(ad[boo[0]])
            times_igrf_downleg[i] = td[boo[0][maxloc]]
        else:
            times_igrf_downleg[i] = 0 

    times_igrf_downleg[347:] = 0

    times_igrf_downleg = times_igrf_downleg[times_igrf_downleg != 0]


    import matplotlib.pyplot as plt 
    plt.plot(times_igrf_downleg)
    #-----------------------------------------------

    goo = np.where(times_igrf_upleg != 0)[0]
    times_upleg = times_igrf_upleg[goo]
    Bo_upleg = np.asarray(df['Bo'][goo])
    #goo = np.where(times_igrf_downleg == np.min(times_igrf_downleg))[0][0]
    goo = np.where(times_igrf_downleg != 0)[0]
    times_downleg = times_igrf_downleg[goo]
    Bo_downleg = np.asarray(df['Bo'][goo])


    d = {'times_upleg':times_upleg, 'times_downleg':times_downleg,'Bo_upleg':Bo_upleg,'Bo_downleg':Bo_downleg}
    #dfu = pd.DataFrame(data=d)
    #df = pd.concat([dfu,df],axis=1,join='inner')
    #plt.plot(df2['times_upleg'],df2['Height(km)'])

    return d

#    #If specific alt is requested...
#    if alt != 0:
#        df = df[df['Height(km)'] >= alt]
#        return df.iloc[0]
#    else: 
#        return df
"""


#----------------------------------------
#Load IRI data (single time or all)
#Note that the official IRI run on the Endurance Box website is set up for the location/conditions at the beginning of the mission. 
#--Year= 2022., Month= 05, Day= 11, Hour=1.5,
#--Time_type = Universal
#--Coordinate_type = Geographic
#--Latitude= 78.9235, Longitude= 11.909, Height= 100.


#(0) Height, km
#(1) Electron_density_Ne, m-3
#(2) Ti, K
#(3) Te, K
#(4) O_ions, %
#(5) H_ions, %
#(6) He_ions, %
#(7) O2_ions, %
#(8) NO_ions, %
#(9) N_ions, %
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


