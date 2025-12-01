#Convert IRI model values to mission times


#dict_keys(['H_km', 'Ne_cm-3', 'Ne_NmF2', 'Tn_K', 'Ti_K', 'Te_K', 'O+', 'N+', 'H+', 'He+', 'O2+', 'NO+', 'Clust', 'TEC_1E16m-2', 't_%'])

def iri_vals_to_mission_times(iridat, mission_times, mission_altitudes):
        
    import numpy as np 

    #For each time of ephem altitude data find the value of fpe from IRI at that altitude
    for i in range(len(mission_times)):
        alt = mission_altitudes[i]
        #Find nearest altitude in IRI data
        idx = (np.abs(iridat['H_km'] - alt)).argmin()
 
        if i == 0:
            H_km_interp = [iridat['H_km'][idx]]
            Ne_cm3_interp = [iridat['Ne_cm-3'][idx]]
            Ne_NmF2_interp = [iridat['Ne_NmF2'][idx]]
            Tn_K_interp = [iridat['Tn_K'][idx]]
            Ti_K_interp = [iridat['Ti_K'][idx]]
            Te_K_interp = [iridat['Te_K'][idx]]
            Oplus_interp = [iridat['O+'][idx]]
            Nplus_interp = [iridat['N+'][idx]]
            Hplus_interp = [iridat['H+'][idx]]
            Heplus_interp = [iridat['He+'][idx]]
            O2plus_interp = [iridat['O2+'][idx]]
            NOplus_interp = [iridat['NO+'][idx]]
        else:
            H_km_interp = np.append(H_km_interp, iridat['H_km'][idx])
            Ne_cm3_interp = np.append(Ne_cm3_interp, iridat['Ne_cm-3'][idx])
            Ne_NmF2_interp = np.append(Ne_NmF2_interp, iridat['Ne_NmF2'][idx])
            Tn_K_interp = np.append(Tn_K_interp, iridat['Tn_K'][idx])
            Ti_K_interp = np.append(Ti_K_interp, iridat['Ti_K'][idx])
            Te_K_interp = np.append(Te_K_interp, iridat['Te_K'][idx])
            Oplus_interp = np.append(Oplus_interp, iridat['O+'][idx])
            Nplus_interp = np.append(Nplus_interp, iridat['N+'][idx])
            Hplus_interp = np.append(Hplus_interp, iridat['H+'][idx])
            Heplus_interp = np.append(Heplus_interp, iridat['He+'][idx])
            O2plus_interp = np.append(O2plus_interp, iridat['O2+'][idx])
            NOplus_interp = np.append(NOplus_interp, iridat['NO+'][idx])


    plt.plot(mission_times, Oplus_interp, label='O+ from IRI model', color='blue', linestyle='--')
    plt.xlabel('Time (s)')



#Load sample IRI data
filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/IRI/iri_output_36381.txt'
iridat = iri_vals_to_mission_times(filename)
fpeiri = 8980.*np.sqrt(iridat['Ne_cm-3'])  #in Hz

