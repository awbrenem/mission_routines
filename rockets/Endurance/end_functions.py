"""
Various useful functions for plotting/analyzing Endurance sounding rocket data
"""


#Return harmonics of an array of cyclotron frequencies 
#fc --> array of cyclotron frequencies 
#nvals --> desired harmonics 

#Return is size [fc, nvals]

def harmonics(fc, nvals):
    v = np.zeros((len(fc), len(nvals)))
    for cnt, n in enumerate(nvals):
        v[:,cnt] = [f*n for f in fc]

    return v


