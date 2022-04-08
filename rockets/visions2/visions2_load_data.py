#Load VISIONS 2 rocket data 


#If you'd rather run the module as a script (load everything):
#%run visions2_load_data.py



from scipy.io import readsav 
from matplotlib import interactive 
interactive(True)
import numpy as np


#General data path 
path = "/Users/abrenema/Desktop/Research/Rocket_missions/VISIONS-2/VISIONS2_data/"


"""
Load Langmuir probe data

NOTE: currently only low-flyer.
"""
def load_langmuir():
    import csv

    fn = [path + "Low-Flyer/LP/35040_FixedBiasDensity_forAaron.csv",path + "High-Flyer/LP/35039_FixedBiasDensity_forAaron.csv"]

    for i in range(2):

        data = []
        alt = []
        times = []
        ne = []

        with open(fn[i], 'r') as file:
            reader = csv.reader(file, delimiter = '\t')
            header = next(reader)
            for row in reader:
                data.append(row)
                goo = row[0].split(",")
                times.append(float(goo[0]))
                alt.append(float(goo[1]))
                ne.append(float(goo[2]))

        file.close()

        if i == 0: 
            lf = {"times":times, "alt":alt, "ne":ne}
        if i == 1: 
            hf = {"times":times, "alt":alt, "ne":ne}


    return {"lowflyer":lf, "highflyer":hf}

"""
Load Efield DC data
Returns zonal and meridional E-field data and flow velocities

Currently I have only Henry's low-flyer solution (64 S/sec) from Jun 14, 2019. No high-flyer data (bad deployment)
"""

def load_efieldDC():

    edc = readsav(path + 'Low-Flyer/efield_DC/35040_unsmoothed_Esoln_v2.sav')
    vals = {"Emerid":edc["emer"], "Ezonal":edc["ezon"], "Vmerid":edc["vmer"], "Vzonal":edc["vzon"], "times":edc["time"]}
    return vals



"""
Load mag data
Returns magnetic field and various cyclotron freqs

NOTE: currently only loads Bonalsky mag data in "xyz" coord.


"""
def load_mag():
    print('LOAD MAG: Currently low flyer only!')

    #mag_bart_nwu = readsav(path + 'Mag/' + '35040_Bart_mag_NWU.sav') 
    #mag_bart_xyz = readsav(path + 'Mag/' + '35040_Bart_mag_XYZ.sav')
    #mag_bon_nwu_lmc = readsav(path + 'Mag/' + '35040_BON_mag_model_NWU_LMC.sav')
    mag_bon_xyz = readsav(path + 'Low-Flyer/Mag/' + '35040_Bonalsky_mag_XYZ.sav')
    mag_dB = readsav(path + 'Low-Flyer/Mag/' + '35040_dB_LMC_NWU_corrected.sav')
    #dict_keys(['memo', 'btime', 'dbmer', 'dbzon', 'dbpar', 'dbn', 'dbw', 'dbu'])

    times = mag_bon_xyz['t']
    timesdB = mag_dB["btime"]

    bo = np.sqrt(mag_bon_xyz['xmag']**2 + mag_bon_xyz['ymag']**2 + mag_bon_xyz['zmag']**2)
    fce = [28.*bo[i] for i in range(len(bo))]
    fcH = [fce[i]/1836. for i in range(len(fce))]
    fcO = [fcH[i]/15. for i in range(len(fce))]

    vals = {"mag_bon_xyz":mag_bon_xyz, "times":times, "bo":bo, "fce":fce, "fcH":fcH, "fcO":fcO, 
            "timesdB":timesdB, "dBmerid":mag_dB["dbmer"], "dBzonal":mag_dB["dbzon"], "dBpar":mag_dB["dbpar"],
            "dBn":mag_dB["dbn"],"dBw":mag_dB["dbw"],"dBu":mag_dB["dbu"]}

    return vals


"""
Load particle data (EEA electrons or EIA ions). 

NOTE: Currently only low-flyer!

times,energies,pitch-angles: (12480, 49, 20)
EEA (electrons):
   Energies from 3 - 30000 eV
   PA from -180 to 180 ;
EIA (ions):
   Energies from 1.7 - 16000 eV
   PA from -180 to 180


Usage:
    typeload = "eea"  or "eia"
    pdat = load_particle(typeload)

If range and sumtype are set, then code will integrate over energies or pitch angles 
e.g. 
    #integrate over pitch angles from -10 to 10 deg
    sumtype = 1
    range = [-10,10]
    pdat = load_particle(typeload,range,sumtype)

    #integrate over energies from 1000 to 2000 eV
    sumtype = 2
    range = [1000,2000]

"""
def load_particle(typeload, range="0",sumtype="0",FoldPitchangles="0"):


    print('LOAD PARTICLE - currently low flyer only ')
    if typeload == "eea":
        dataAll = readsav(path + "Low-Flyer/particles/" + 'VISIONS2EEAShiresGSFC.sav')
    if typeload == "eia":
        dataAll = readsav(path + "Low-Flyer/particles/" + 'VISIONS2EIAShiresGSFC.sav')
    if typeload != "eea" and typeload != "eia":
        print("WRONG INPUT. MUST SELECT 'eea' or 'eia'")
        return -1 

######################
    #Fold + and - pitch-angles together, if requested 
    if FoldPitchangles != "0":


        if typeload == "eea":
        #Create folded array - this has 12 unique pitch angles [0, 5, 20, 45, 70, 85, 90, 95, 110, 135, 160, 180]

            goototals = [dataAll[typeload+"s_hires_eflux_gsfc"][:,:,0], # 0 deg
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,1] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,2])/2.,  #5 and -5 
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,3] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,4])/2.,  #20 and -20 
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,5] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,6])/2.,  #45 and -45
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,7] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,8])/2.,  #70 and -70
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,9],  #85
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,10] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,11])/2.,  #90 and -90
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,12],  #95
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,13] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,14])/2.,  #110 and -110
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,15] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,16])/2.,  #135 and -135
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,17] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,18])/2.,  #160 and -160
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,19]]  #180

            #Construct new dataAll array with the summed pitch-angles
            #The energy channels will be the same 
            dataAll[typeload+"s_hires_eflux_gsfc"] = np.transpose(goototals,axes=[1,2,0])
            dataAll[typeload+"s_hires_pitchnom_gsfc"] = [0,5,20,45,70,85,90,95,110,135,160,180] #Define new pitch angles


        if typeload == "eia":
        #Create folded array - this has 12 unique pitch angles [0, 5, 20, 45, 70, 85, 90, 95, 110, 135, 160, 180]


            goototals = [dataAll[typeload+"s_hires_eflux_gsfc"][:,:,0], # 0 deg
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,1] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,2])/2.,  #10 and -20
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,3] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,4])/2.,  #30 and -45 
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,5] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,6])/2.,  #60 and -70
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,7],  #80
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,8] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,9])/2.,  #-90 and 90
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,10],  #100
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,11] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,12])/2.,  #110 and -115
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,13],  #120
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,14],  #130
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,15],  #140
                (dataAll[typeload+"s_hires_eflux_gsfc"][:,:,16] + dataAll[typeload+"s_hires_eflux_gsfc"][:,:,17])/2.,  #155 and -155
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,18],  #170
                dataAll[typeload+"s_hires_eflux_gsfc"][:,:,19]]  #180

            #Construct new dataAll array with the summed pitch-angles
            #The energy channels will be the same 
            dataAll[typeload+"s_hires_eflux_gsfc"] = np.transpose(goototals,axes=[1,2,0])
            dataAll[typeload+"s_hires_pitchnom_gsfc"] = [0,15,25,37.5,65,80,90,100,112.5,120,130,140,155,170,180] #Define new pitch angles


######################


    #integrate over energies or pitch-angles, if requested
    if sumtype != "0":

        goo = []
        datar = []
        axis = 0

        #sum over pitch-angles
        if sumtype == 1:
            goo = np.where((np.asarray(dataAll[typeload+"s_hires_pitchnom_gsfc"]) >= range[0]) & (np.asarray(dataAll[typeload+"s_hires_pitchnom_gsfc"]) <= range[1]))
            #Get reduced array with only desired elements to be summed over. This allows us to simply sum over an entire axis
            datar = dataAll[typeload+"s_hires_eflux_gsfc"][:,:,goo]
            axis = 2
            dataAll[typeload+"s_hires_pitchnom_gsfc"] = "NaN"

        #sum over energies
        if sumtype == 2:
            goo = np.where((np.asarray(dataAll[typeload+"s_hires_energy_gsfc"]) >= range[0]) & (np.asarray(dataAll[typeload+"s_hires_energy_gsfc"]) <= range[1]))
            datar = dataAll[typeload+"s_hires_eflux_gsfc"][:,goo,:]
            axis = 1
            dataAll[typeload+"s_hires_energy_gsfc"] = "NaN"


        if np.size(goo) != 1:
            datar = np.squeeze(datar)
        if datar.ndim == 4:
            datar = np.squeeze(datar,3)

        #Sum over the desired axis 
        dataAll[typeload+"s_hires_eflux_gsfc"] = np.squeeze(np.sum(datar,axis=axis, keepdims=True))


        if sumtype == 1:
            dataAllFin = {"flux":dataAll[typeload+"s_hires_eflux_gsfc"],
                "times":dataAll[typeload+"s_hires_ftime_gsfc"],
                "energies":dataAll[typeload+"s_hires_energy_gsfc"]}
        if sumtype == 2:
            dataAllFin = {"flux":dataAll[typeload+"s_hires_eflux_gsfc"],
                "times":dataAll[typeload+"s_hires_ftime_gsfc"],
                "pitchangles":dataAll[typeload+"s_hires_pitchnom_gsfc"]}

 
    if sumtype != 1 and sumtype != 2:
        dataAllFin = {"flux":dataAll[typeload+"s_hires_eflux_gsfc"],
            "times":dataAll[typeload+"s_hires_ftime_gsfc"],
            "energies":dataAll[typeload+"s_hires_energy_gsfc"],
            "pitchangles":dataAll[typeload+"s_hires_pitchnom_gsfc"]}


    return dataAllFin


"""
Load VLF data
"""
def load_vlf():
    print('LOAD VLF: Currently low flyer only')
    vlf12_lf = readsav(path + "Low-Flyer/efield_VLF/" + '35040_LFDSP_S1_VLF12_mvm_AaronB.sav') #Low flyer (35.040)
    vlf34_hf = readsav(path + "High-Flyer/efield_VLF/" + '35039_Main_VLF34B_VLF34Boosted.sav') 

    data = {"vlf12_lf":vlf12_lf, "vlf34_hf":vlf34_hf}
    return data

"""
Low flyer:
dict_keys(['author', 'calibrations', 'datafile', 'dvlf12', 'flight', 'format', 't0', 'tvlf12', 'units'])
High flyer:
dict_keys(['channel', 'chdesc', 'dvals', 'sfidvals', 'tvals', 'samplerate', 'timetagmethod', 'flight', 'datafile', 'inputfile', 'format', 'link', 't0', 't0units', 'author', 'date', 'dataunits', 'timeunits', 'firstsfid'])

np.shape(vlf12_lf["dvlf12"])
(15918022,)

np.shape(vlf34_hf["dvals"])
(7990127, 6)
np.shape(vlf34_hf["sfidvals"])
(7990127,)
np.shape(vlf34_hf["tvals"])
(7990127,)

"""


#allow to be run as a script
if __name__ == '__main__':
    print("<> Running vision2_load_data.py as a script! <>")
    #x = load_particle("eea")
    #print(x)
    x = load_mag()
    x.keys()
    print("here_fin")


