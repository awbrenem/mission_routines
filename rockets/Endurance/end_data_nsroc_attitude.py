"""
Read in NSROC's Endurance attitude solution file.

Of relevance:
X_Az --> azimuth relative to geographic North of the nosecone
Z_Az --> azimuth relative to geographic North of the rail

"""


import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

fn = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/ephemeris/47.001_SynchronizedSolution.csv'
nms = ['Time','Yaw','Pitch','Roll','RollRate','AoA_T','AngleB','a11','a12','a13','a21','a22','a23','a31','a32','a33','X_Az','X_El','Y_Az','Y_El','Z_Az','Z_El','Latgd','Long','Alt']
dat = pd.read_csv(fn,skiprows=7,header=0,names=nms)


last = 43800
times = np.asarray(dat['Time'][0:last])

#Some string values corresponding to improperly labeled NaNs. Turn these to zeros
azim_nose = np.asarray(pd.to_numeric(dat['X_Az'][0:last]))
goo = np.where(azim_nose < 0)
azim_nose[goo[0]] = 360 + azim_nose[goo[0]]

azim_rail = np.asarray(pd.to_numeric(dat['Z_Az'][0:last]))


plt.plot(times,azim_rail)
plt.ylim(-180,180)


plt.plot(times,azim_nose)
plt.ylim(0,360)




