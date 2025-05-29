#wave polarization test
import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_hodogram_dynamic as hod
import numpy as np
import matplotlib.pyplot as plt

nperiods = 10
dt = 0.01

times = np.arange(0,20,dt)

f1 = 1 
f2 = 0.9*f1
wt1 = 2*3.14*f1 * times 
wt2 = 2*3.14*f2 * times 



#v = np.arange(0,6.28*nperiods,dt)
#v2 = np.arange(0,6*nperiods,dt)
#al = 2 #amplitude of linear waves
ac = 1 #amplitude of circular waves
phase_l = 45. * (3.14/180) #phase of wave 1 relative to x-axis
phase_l2 = 230. * (3.14/180) #phase of wave 2 relative to x-axis
phase_t = 5. * (3.14/180) #phase b/t x and y components of each wave (zero for linear)
#phase_t = 2. * (3.14/180) #phase b/t x and y components of each wave (zero for linear)
a1 = 1
a2 = 0.2
#al12 = 1 
#al34 = 1
#bl12 = 0.5
#bl34 = 0.5

#Define waves in geophysical coord then project onto sc coord

wx1 = a1*np.cos(phase_l)*np.cos(wt1)
wx2 = a1*np.cos(phase_l2)*np.cos(wt2)

wy1 = a2*np.sin(phase_l)*np.cos(wt1 + phase_t)
wy2 = a2*np.sin(phase_l2)*np.cos(wt2 + phase_t)


plt.plot(wx1,wy1)
plt.plot(wx2,wy2)

hod.plot_hodogram_dynamic(wx1,wy1, npts=2, gap=5,pauseT=0.01)
hod.plot_hodogram_dynamic(wx1 + wx2,wy1+ wy2, npts=2, gap=5,pauseT=0.05)




#wx = np.cos(phase_l)*np.cos(v)
#wtst = np.cos(phase_l)*np.cos(v2)

#wy = np.sin(phase_l)*np.cos(v)
#wx2 = np.cos(phase_l2)*np.cos(v + phase_t)
#wy2 = np.sin(phase_l2)*np.cos(v + phase_t)

#hod.plot_hodogram_dynamic(wx + wx2,wy + wy2, npts=20, gap=5,pauseT=0.005)
#hod.plot_hodogram_dynamic(wx,wy, npts=20, gap=5,pauseT=0.005)
#hod.plot_hodogram_dynamic(wx2,wy2, npts=20, gap=5,pauseT=0.005)


#Q1
wf12 = 1*(wx1+wx2)
wf34 = 1*(wy2+wy2)
hod.plot_hodogram_dynamic(wf12,wf34, npts=20, gap=5,pauseT=0.005)

#Q2
wf12 = 1*(wy1+wy2)
wf34 = -1*(wx1+wx2)
hod.plot_hodogram_dynamic(wf12,wf34, npts=20, gap=5,pauseT=0.005)

#Q3
wf12 = -1*(wx1+wx2)
wf34 = -1*(wy1+wy2)
hod.plot_hodogram_dynamic(wf12,wf34, npts=20, gap=5,pauseT=0.005)

#Q4
wf12 = -1*(wy1+wy2)
wf34 = 1*(wx1+wx2)
hod.plot_hodogram_dynamic(wf12,wf34, npts=20, gap=5,pauseT=0.005)
























#wave 1
l12_q1 = al12*np.cos(v)
l12_q2 = al12*np.cos(v)
l12_q3 = al12*-1*np.cos(v)
l12_q4 = al12*-1*np.cos(v)
l34_q1 = al34*np.cos(v)
l34_q2 = al34*-1*np.cos(v)
l34_q3 = al34*-1*np.cos(v)
l34_q4 = al34*np.cos(v)

#wave 2
l12_r1 = -1*bl12*np.sin(v+phase_l)
l34_r1 = bl34*np.sin(v+phase_l)









hod.plot_hodogram_dynamic(wx,wy, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(w2x,w2y, npts=20, gap=5,pauseT=0.005)


hod.plot_hodogram_dynamic(l12_q1,l34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(l12_r1,l34_r1, npts=20, gap=5,pauseT=0.005)


plt.plot(l12_q1,l34_q1)
plt.plot(l12_r1,l34_r1)


supl12_q1 = l12_q1 + l12_r1
supl34_q1 = l34_q1 + l34_r1






hod.plot_hodogram_dynamic(supl12_q1,supl34_q1, npts=20, gap=5,pauseT=0.005)







c12_q1 = ac*np.sin(v)
c12_q2 = ac*np.cos(v)
c12_q3 = ac*-1*np.sin(v)
c12_q4 = ac*-1*np.cos(v)

c34_q1 = ac*np.cos(v)
c34_q2 = ac*-1*np.sin(v)
c34_q3 = ac*-1*np.cos(v)
c34_q4 = ac*np.sin(v)

s12_q1 = l12_q1 + c12_q1
s12_q2 = l12_q2 + c12_q2
s12_q3 = l12_q3 + c12_q3
s12_q4 = l12_q4 + c12_q4

s34_q1 = l34_q1 + c34_q1
s34_q2 = l34_q2 + c34_q2
s34_q3 = l34_q3 + c34_q3
s34_q4 = l34_q4 + c34_q4

hod.plot_hodogram_dynamic(l12_q1,l34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(l12_q2,l34_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(l12_q3,l34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(l12_q4,l34_q4, npts=20, gap=5,pauseT=0.005)

hod.plot_hodogram_dynamic(c12_q1,c34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(c12_q2,c34_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(c12_q3,c34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(c12_q4,c34_q4, npts=20, gap=5,pauseT=0.005)

hod.plot_hodogram_dynamic(s12_q1,s34_q1, npts=20, gap=5,pauseT=0.005) #elliptical RH
hod.plot_hodogram_dynamic(s12_q2,s34_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(s12_q3,s34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(s12_q4,s34_q4, npts=20, gap=5,pauseT=0.005)

#------------------------
#Plots in geophysical coordinates (preferred coord of the waves)
#All of these should show the same exact polarization!!!
hod.plot_hodogram_dynamic(l12_q1,l34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*l34_q2,l12_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*l12_q3,-1*l34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(l34_q4,-1*l12_q4, npts=20, gap=5,pauseT=0.005)

hod.plot_hodogram_dynamic(c12_q1,c34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*c34_q2,c12_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*c12_q3,-1*c34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(c34_q4,-1*c12_q4, npts=20, gap=5,pauseT=0.005)

hod.plot_hodogram_dynamic(s12_q1,s34_q1, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*s34_q2,s12_q2, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(-1*s12_q3,-1*s34_q3, npts=20, gap=5,pauseT=0.005)
hod.plot_hodogram_dynamic(s34_q4,-1*s12_q4, npts=20, gap=5,pauseT=0.005)



