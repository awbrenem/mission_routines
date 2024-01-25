"""
Rotate Endurance data to geophysical system (W, N, U) - aka West, North, Up.
NOTE: This is similar ECEF (Earth-centered, Earth-fixed) coord - up to some minus signs - that the Endurance rocket velocity data come in. 
NOTE: Coefficients Henry uses come from NSROC attitude solution (47.001_SynchronizedSolution.csv)

-----------------------------
NOTE on testing:  (see end_data_nsroc_attitude.py for comparison with direct NSROC values)
What I know works: 
    Comparison with the NSROC azimuth of the rail wrt geographic N (end_data_nsroc_attitude.py) and my version of this is good. 
    The same comparison but with V1 shows that V1 is always off by 20 deg. 
    So, I know that the probe rotations into WNU coordinates are working. 

What kinda works:
    The elevation angle. I only get positive values. When I try to use the signed values something doesn't work out right. 

What doesn't work: 
    recreation of Henry's "Azimuth of the x-axis" plot. I.e. the nose azimuth wrt geographic N.
-----------------------------



(based on Henry Freudenreich's endurattfrag.pro)

input: vx,vy,vz - vector components in rocket coordinates, defined as: 
x-hat = nose cone 
y-hat = completes RH coord
z-hat = rail 
(see coord system in Endurance pfr47001_report_v6.pdf for proof that z-axis is along rail)


The rocket probes are shifted by 22.5 deg to the rail (unusual for a rocket mission).
For example, V1 is 22.5 deg counterclockwise to the rail. 

------------------------------------
Useful coordinates on launcher are: 

angle = 22.5 deg 

#rail
    bx = 0
    by = 0 
    bz = 1
#nose
    bx = 1
    by = 0 
    bz = 0
#V1 probe vector
    bx = 0 
    by = -np.sin(angle)
    bz = np.cos(angle)
#V2 probe vector
    bx = 0
    by = np.sin(angle)
    bz = -np.cos(angle)
#V3 probe vector 
    bx = 0 
    by = np.cos(angle)
    bz = np.sin(angle)
#V4 probe vector 
    bx = 0 
    by = -np.cos(angle)
    bz = -np.sin(angle)

#DESI boom 1 
    bx = 0 
    by = np.sin(angle)
    bz = np.cos(angle)
#DESI boom 2 
    bx = 0 
    by = np.cos(angle)
    bz = -np.sin(angle)
#DESI boom 3 
    bx = 0 
    by = -np.sin(angle)
    bz = -np.cos(angle)
#DESI boom 4 
    bx = 0 
    by = -np.cos(angle)
    bz = np.sin(angle)

#rail vector in WNU 
    rail_wnu = [0.682, -0.731, 0]




------------------------------------



Launcher Settings:
==================
Launcher:               U4
Latitude:    78.931500 deg
Longitude:   11.850400 deg
Altitude:         111.5 ft
Elevation:        84.3 deg
Azimuth:         222.2 deg  (angle relative to North direction in WNU coord)
Bank:              0.0 deg
 

"""



#def end_cal_rotate_to_geophysical(vx,vy,vz):

from os.path import dirname, join as pjoin
from scipy.io import readsav
import numpy as np
import matplotlib.pyplot as plt
import math

#Load Henry's rotation matrices (he gets these from the NSROC solution)
path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/data/rotation_matrices/'
sav_fname = pjoin(path,'47001_orig_att.sav')
dat = readsav(sav_fname)


#angle b/t V1 and rail
angle = np.radians(22.5)

#TOF of interest
#t=0
t = np.arange(0,900,1)

#signed or unsigned angles?
signed = True

#-------------------------------------------------------
#Define the vectors corresponding to the probe directions.
#Used to determine their WNU coordinates vs time
#-------------------------------------------------------

#rail
#bx = 0
#by = 0 
#bz = 1

#V1
bx = 0 
by = -np.sin(angle)
bz = np.cos(angle)

#V2
#bx = 0
#by = np.cos(np.radians(90)-angle)
#bz = -np.cos(angle)

#V3
#bx = 0 
#by = np.cos(angle)
#bz = np.sin(angle)

#V4
#bx = 0 
#by = -np.cos(angle)
#bz = -np.sin(angle)



#--------------------------------------------
#test to see if this is defined correctly
#vtst = [bx,by,bz]
#mag = np.linalg.norm(vtst)
#vtst = vtst / mag
#vtst2 = [0,0,1]
#vtmp = [vtst[i]*vtst2[i] for i in range(len(vtst))]
#degv = np.degrees(np.arccos(vtmp))  #22.5 deg 
#--------------------------------------------

bw = np.zeros(len(t))
bn = np.zeros(len(t))
bu = np.zeros(len(t))

elev = np.zeros(len(t))
azi_rail = np.zeros(len(t))
azi_rail2 = np.zeros(len(t))
azi_north = np.zeros(len(t))

atime = dat['atime']
a11 = dat['a11']
a12 = dat['a12']
a13 = dat['a13']
a21 = dat['a21']
a22 = dat['a22']
a23 = dat['a23']
a31 = dat['a31']
a32 = dat['a32']
a33 = dat['a33']


for i in range(len(t)):

    #Interpolate coefficients to time of interest
    c11=np.interp(t[i],atime,a11)
    c12=np.interp(t[i],atime,a12)
    c13=np.interp(t[i],atime,a13)
    c21=np.interp(t[i],atime,a21)
    c22=np.interp(t[i],atime,a22)
    c23=np.interp(t[i],atime,a23)
    c31=np.interp(t[i],atime,a31)
    c32=np.interp(t[i],atime,a32)
    c33=np.interp(t[i],atime,a33)

    #plt.plot(atime,a11)
    #plt.plot(t,c11,'*')
    #plt.show()

    # normalize
    r = np.sqrt(c11**2 + c21**2 + c31**2)
    c11 /= r 
    c21 /= r 
    c31 /= r 
    r = np.sqrt(c12**2 + c22**2 + c32**2)
    c12 /= r 
    c22 /= r 
    c32 /= r 
    r = np.sqrt(c13**2 + c23**2 + c33**2)
    c13 /= r 
    c23 /= r 
    c33 /= r 

    bw[i] = -(c11*bx+c12*by+c13*bz)
    bn[i] = c21*bx+c22*by+c23*bz
    bu[i] = c31*bx+c32*by+c33*bz



    vtst = [bw[i],bn[i],bu[i]] / np.linalg.norm([bw[i],bn[i],bu[i]])

    #calculate elevation angle 
    vtst2 = [bw[i],bn[i],0] 
    vtst2 = vtst2 / np.linalg.norm(vtst2)
    elev[i] = np.degrees(np.arccos(np.dot(vtst, vtst2)))


    #calculate signed (-180,180) azimuth angle relative to the rail
    rail_wnu = [0.682, -0.731, 0]
    if signed == True:
        cosTh = np.dot(rail_wnu,vtst)
        sinTh = np.cross(rail_wnu,vtst)
        tmp = np.rad2deg(np.arctan2(sinTh,cosTh))
        azi_rail[i] = tmp[2]
    else:
        azi_rail[i] = np.degrees(np.arccos(np.dot(vtst, rail_wnu)))

    #calculate signed (-180,180) azimuth angle relative to geographic North
    if signed == True:
        cosTh = np.dot([0,1,0],vtst)
        sinTh = np.cross([0,1,0],vtst)
        tmp = np.rad2deg(np.arctan2(sinTh,cosTh))
        azi_north[i] = tmp[2]
    else:
        azi_north[i] = np.degrees(np.arccos(np.dot(vtst, [0,1,0])))



#plt.plot(t,elev)


#plt.plot(t,azi_rail)
#plt.plot(t,azi_rail2)

print('h')

#fig, axs = plt.subplots(3)
#axs[0].plot(t,azi_north)
#axs[1].plot(t,azi_north2)
#axs[2].plot(t,elev)
##axs[0].set_xlim(0,200)
#axs[2].set_ylim(-90,90)

#Test angle relative to some reference 
#plt.plot(t,azi_rail)
#plt.plot(times,azim_nose)
plt.plot(times,azim_rail)
plt.plot(t,azi_north)
plt.ylim(-180,180)








