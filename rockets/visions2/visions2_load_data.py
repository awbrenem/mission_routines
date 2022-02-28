#Load VISIONS 2 rocket data 


#If you'd rather run the module as a script (load everything):
#%run visions2_load_data.py



from scipy.io import readsav 
from matplotlib import interactive 
interactive(True)
import numpy as np


#General data path 
path = "/Users/abrenema/Desktop/Research/Rocket_missions/VISIONS-2/VISIONS2_data/"



def load_langmuir():
    import csv
    
    fn = path + "Low-Flyer/LP/35040_FixedBiasDensity_for Aaron.csv"

    header = []
    data = []
    alt = []
    times = []
    ne = []

    with open(fn, 'r') as file:
        reader = csv.reader(file, delimiter = '\t')
        header = next(reader)
        for row in reader:
            data.append(row)
            goo = row[0].split(",")
            times.append(float(goo[0]))
            alt.append(float(goo[1]))
            ne.append(float(goo[2]))

    file.close()

    #       time, alt, e- density 
    #fn = "High_Flyer/LP/35039_FixedBiasDensity_for Aaron.csv"

    return {"times":times, "alt":alt, "ne":ne}


def load_mag():
    print('LOAD MAG: Currently low flyer only!')

    #mag_bart_nwu = readsav(path + 'Mag/' + '35040_Bart_mag_NWU.sav') 
    #mag_bart_xyz = readsav(path + 'Mag/' + '35040_Bart_mag_XYZ.sav')
    #mag_bon_nwu_lmc = readsav(path + 'Mag/' + '35040_BON_mag_model_NWU_LMC.sav')
    mag_bon_xyz = readsav(path + 'Low-Flyer/Mag/' + '35040_Bonalsky_mag_XYZ.sav')
    times = mag_bon_xyz['t']

    bo = np.sqrt(mag_bon_xyz['xmag']**2 + mag_bon_xyz['ymag']**2 + mag_bon_xyz['zmag']**2)
    fce = [28.*bo[i] for i in range(len(bo))]
    fcH = [fce[i]/1836. for i in range(len(fce))]
    fcO = [fcH[i]/15. for i in range(len(fce))]

    vals = {"mag_bon_xyz":mag_bon_xyz, "times":times, "bo":bo, "fce":fce, "fcH":fcH, "fcO":fcO}

    return vals


"""
Load particle data. 
If range and sumtype are set, then code will integrate over energies or pitch angles 
e.g. 
    #integrate over pitch angles from -10 to 10 deg
    sumtype = 1
    range = [-10,10]

    #integrate over energies from 1000 to 2000 eV
    sumtype = 2
    range = [1000,2000]

NOTE: need to finish for EEI data
"""
def load_particle(range="0",sumtype="0"):

    print('LOAD PARTICLE - currently low flyer only ')
    eea_data = readsav(path + "Low-Flyer/particles/" + 'VISIONS2EEAShiresGSFC.sav')
    eia_data = readsav(path + "Low-Flyer/particles/" + 'VISIONS2EIAShiresGSFC.sav')


    #integrate over energies or pitch-angles, if requested
    if sumtype != "0":

        goo = []
        datar = []
        axis = 0

        #sum over pitch-angles
        if sumtype == 1:
            goo = np.where((np.asarray(eea_data["eeas_hires_pitchnom_gsfc"]) >= range[0]) & (np.asarray(eea_data["eeas_hires_pitchnom_gsfc"]) <= range[1]))
            #Get reduced array with only desired elements to be summed over. This allows us to simply sum over an entire axis
            datar = eea_data["eeas_hires_eflux_gsfc"][:,:,goo]
            axis = 2
            eea_data["eeas_hires_pitchnom_gsfc"] = "NaN"

        if sumtype == 2:
            goo = np.where((np.asarray(eea_data["eeas_hires_energy_gsfc"]) >= range[0]) & (np.asarray(eea_data["eeas_hires_energy_gsfc"]) <= range[1]))
            datar = eea_data["eeas_hires_eflux_gsfc"][:,goo,:]
            axis = 1
            eea_data["eeas_hires_energy_gsfc"] = "NaN"


        if np.size(goo) != 1:
            datar = np.squeeze(datar)
        if datar.ndim == 4:
            datar = np.squeeze(datar,3)

        #Sum over the desired axis 
        eea_data["eeas_hires_eflux_gsfc"] = np.squeeze(np.sum(datar,axis=axis, keepdims=True))




    return {"eea_data":eea_data, "eia_data":eia_data}


def load_vlf():
    print('LOAD VLF: Currently low flyer only')
    vlf12_data = readsav(path + "Low-Flyer/efield_VLF/" + '35040_LFDSP_S1_VLF12_mvm_AaronB.sav') #Low flyer (35.040)
    #vlf34_HighFlyer_data = readsav(path + "High-Flyer/efield_VLF/" + '35039_Main_VLF34B_VLF34Boosted.sav') 
    return vlf12_data



def load_extra():
    print('Hi from within load_extra!')


#allow to be run as a script
if __name__ == '__main__':
    print("<> Running vision2_load_data.py as a script! <>")
    x = load_langmuir()
    print(x)
