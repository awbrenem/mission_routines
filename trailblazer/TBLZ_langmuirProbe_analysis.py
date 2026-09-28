# **********************************
"""
Simulate expected I-V curves for GDC Trailblazer
Use to estimate max and min currents for amplifier selection
Use to estimate sampling conops


GDC Trailblazer conditions
altitudes from 350-500 km

We're on the hook to measure densities from 10^2 - 10^7 cm-3 and Te from 300 - 10000 K.

"""

import sys

sys.path.append(
    "/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/mission_routines/trailblazer/"
)
sys.path.append("/Users/abrenema/Desktop/code/Aaron/github/antenna_analysis/")
import matplotlib.pyplot as plt
import numpy as np
from langmuirProbes import *
from matplotlib.collections import LineCollection

# **********************************

# SLP Probe properties
betaVal = 0.5  # cylindrical probe
r = 0.32  # cm (cylindrical probe radius)
l = 5  # cm (probe length)
## Effective area of cylindrical probe - i.e. rectangular cross section.
Aeff = 2 * r * l


# C/NOFS example
# betaVal = 1  # sphere
# r = 2.5 / 2  # cm (spherical probe radius)
# Aeff = np.pi * r**2  # sunlit area only


# Fixed bias probe (patch probe)
# betaVal = 0  # planar probe
# r = 3.95 / 2.0  # cm (planar probe radius)
## r = 6 / 2.0  # cm (planar probe radius)
# Aeff = np.pi * r**2  # cm^2 (probe area)


B = 50000  # nT
VelSC = 7600  # m/sec
# VelSC = 0  # m/sec


# Range of densities we're required to work over is 10^2 - 10^7 cm-3
ne = 1e7  # cm^-3
ni = ne

# Range of electron temps from 400-10000 K (0.03-0.86 eV)
Te = 0.86  # eV
# Te = 0.03  # eV
# Te = 0.3
Ti = Te * 0.8
# Ti = 0.05  # eV
amu = 16  # Pure O+


# See Asenovski+18 (https://doi.org/10.1016/j.jastp.2017.12.013) for a discussion of photoemission current densities
# Iph0 = 3  # nA/cm^2 (current density estimate for Gold [e.g. Hartzell+23])
Iph0 = 6  # nA/cm^2 (estimate for TiN [Laakso+26; Pederson+08])


# Debye length (few cm)
lambda_D = debye_length(ne, Te)  # m


Vf = potentialFloating_NoPhotoelectrons(Te, amu)


# Thermal electron and ion current (Barjatya +2007 thesis Eqn 1.2)
Ithe = currentThermalElectron(ne, Te, Aeff)
Ithi = currentThermalIon(ni, Ti, Aeff, amu)


dVr = np.arange(
    -6, 12, 0.001
)  # probe potential relative to spacecraft potential (V - Vp)

Iph = np.empty(len(dVr))
Ie_ret = np.empty(len(dVr))
Isat_e = np.empty(len(dVr))
Isat_i = np.empty(len(dVr))
Isat_eMeso = np.empty(len(dVr))
Isat_iMeso = np.empty(len(dVr))
Iram_ion = np.empty(len(dVr))
# Iram_el = np.empty(len(dVr))


for i in range(len(dVr)):
    # Relevant b/t Vf and Vp
    Ie_ret[i] = currentElectrons_ElectronRetardation(ne, Te, dVr[i], Aeff)

    # Photoelectron current
    TePhoto = 5  # (1-10 eV typically)
    Iph[i] = currentPhotoelectron(Iph0, TePhoto, dVr[i], Aeff)  # nA

    # Electron and ion saturation currents
    Isat_e[i] = currentElectronSaturation(ne, Te, dVr[i], Aeff, betaVal)
    Isat_eMeso[i] = currentMesothermalElectronSaturation(ne, Te, VelSC, dVr[i], Aeff)

    Isat_i[i] = currentIonSaturation(ni, Ti, dVr[i], Aeff, betaVal, amu)
    Isat_iMeso[i] = currentMesothermalIonSaturation(ni, Ti, VelSC, dVr[i], Aeff, amu)

    Iram_ion[i] = ramCurrent_ion(ni, Aeff, VelSC)  # nA

# -----------------------------------------------------------------------------------------------------


# sum all the currents
Itot = np.empty(len(dVr))
for i in range(len(Itot)):
    Itot[i] = np.nansum([Ie_ret[i], Isat_e[i], Isat_i[i], Iph[i], Iram_ion[i]])


# Change color of total current based on which sub current dominates
current_magnitudes = np.column_stack(
    [np.abs(Ie_ret), np.abs(Isat_e), np.abs(Isat_i), np.abs(Iph), np.abs(Iram_ion)]
)
dominant_current = np.nanargmax(current_magnitudes, axis=1)
dominant_colors = np.array(["green", "red", "blue", "goldenrod", "black"])[
    dominant_current
]


# Find max/min values for some of the plots.
minv = np.nanmin([Ie_ret, Isat_e, Iram_ion, Iph, Isat_i, Itot]) / 1000
maxv = np.nanmax([Ie_ret, Isat_e, Iram_ion, Iph, Isat_i, Itot]) / 1000
minvp = np.nanmin([Ie_ret, Isat_e]) / 1000
maxvp = np.nanmax([Ie_ret, Isat_e]) / 1000
minvm = np.nanmin([Iram_ion, Iph, Isat_i]) / 1000
maxvm = np.nanmax([Iram_ion, Iph, Isat_i]) / 1000


# Plot individual currents
fig, axs = plt.subplots(4)
axs[0].plot(dVr, Ie_ret / 1000, label="Ie_ret", color="green")
axs[0].set_ylabel("Ie_ret (uA)")
axs[0].plot(dVr, Isat_e / 1000, label="Isat_e", color="red")
axs[0].set_ylabel("+ currents (uA)")
axs[1].plot(dVr, np.abs(Isat_i) / 1000, label="|Isat_i|", color="blue")
axs[1].set_ylabel("-Isat_i (uA)")
axs[1].plot(dVr, np.abs(Iph) / 1000, label="|Iph|", color="goldenrod")
axs[1].set_ylabel("-Iph (uA)")
axs[1].plot(dVr, np.abs(Iram_ion) / 1000, label="|Iram ion|", color="black")
axs[1].set_ylabel("- currents (uA)")
itot_values = np.abs(Itot) / 1000
itot_points = np.column_stack([dVr, itot_values])
itot_segments = np.stack([itot_points[:-1], itot_points[1:]], axis=1)
itot_line = LineCollection(itot_segments, colors=dominant_colors[:-1], linewidths=1.5)
axs[2].add_collection(itot_line)
axs[2].plot([], [], color="green", label="Ie_ret")
axs[2].plot([], [], color="red", label="Isat_e")
axs[2].plot([], [], color="blue", label="|Isat_i|")
axs[2].plot([], [], color="goldenrod", label="|Iph|")
axs[2].plot([], [], color="black", label="|Iram ion|")
axs[2].set_ylabel("Itot (uA)")
axs[2].set_xlabel("volts")
axs[3].plot(dVr, Ie_ret / 1000, color="green", label="Ie_ret")
axs[3].plot(dVr, Isat_e / 1000, color="red", label="Isat_e")
axs[3].plot(dVr, Iram_ion / 1000, color="black", label="Iram ion")
axs[3].plot(dVr, Iph / 1000, color="goldenrod", label="Iph")
axs[3].plot(dVr, Isat_i / 1000, color="blue", label="Isat_i")
axs[3].plot(dVr, Itot / 1000, color="purple", linestyle="--", linewidth=2, label="Itot")
for i in range(4):
    axs[i].set_xlim(np.min(dVr), np.max(dVr))
    axs[i].set_yscale("log")
    axs[i].legend()
axs[0].set_ylim(minvp / 100, maxvp * 100)
axs[1].set_ylim(minvm / 100, maxvm * 100)
axs[2].set_ylim(minv / 100, maxv * 100)
axs[3].set_yscale("linear")
axs[3].set_ylim(minv + 0.4 * minv, maxv + 0.4 * maxv)
axs[3].set_ylabel("Currents (uA)\nLinear scale")
axs[3].set_xlabel("Probe Potential Relative to Spacecraft")


# Save total current as pickle file so I can use it in the TBLZ_langmuirProbe_IV_extract.py script to test these curves.
# i.e. make sure that I'm extracting the same ne, Te, ni, Vp that I'm inputting.
import pickle

with open("/Users/abrenema/Desktop/TBLZ_LangmuirProbe_Itot.pkl", "wb") as f:
    pickle.dump((Itot, dVr), f)


# Max current needed in uA
print(np.nanmax(Itot) / 1000)


tst = np.where(dVr > 0)[0][0]


# Relevant b/t Vf and Vp
Ie_ret = currentElectrons_ElectronRetardation(ne, Te, dV, Aeff)


# Electron and ion saturation currents
Isat_e = currentElectronSaturation(ne, Te, dV, Aeff, betaVal)
Isat_eMeso = currentMesothermalElectronSaturation(ne, Te, VelSC, dV, Aeff)

Isat_i = currentIonSaturation(ni, Ti, dV, Aeff, betaVal, amu)
Isat_iMeso = currentMesothermalIonSaturation(ni, Ti, VelSC, dV, Aeff, amu)

# --------------------------------------------------------------------------------
# Find voltage limits required for a range of temperatures and densities.


dVr = np.arange(-1, 1, 0.05)


Ie_ret400k = np.empty(len(dVr))
Ie_ret10000k = np.empty(len(dVr))
for i in range(len(dVr)):
    Ie_ret400k[i] = currentElectrons_ElectronRetardation(4000, 0.03, dVr[i], Aeff)
    Ie_ret10000k[i] = currentElectrons_ElectronRetardation(4000, 0.86, dVr[i], Aeff)


# Plot individual currents
fig, axs = plt.subplots(2)
axs[0].plot(dVr, Ie_ret400k)
axs[0].set_ylim(0, np.nanmax(Ie_ret400k))
axs[1].plot(dVr, Ie_ret10000k)
axs[1].set_ylim(0, np.nanmax(Ie_ret10000k))

for i in range(2):
    axs[i].set_xlim(-1.5, 1.5)
    axs[i].set_yscale("log")
    axs[i].set_ylim(0.1, 1e14)

    # axs[i].set_xlim(np.nanmin(dVr), np.nanmax(dVr))

plt.show()
plt.show()


# Henry example from C/NOFS: 100 pA for a few*10^2/cc for a 2.6 cm diameter sphere and -5V bias
dV = -5
ni = 300
Aeff = 4 * np.pi * (1.3) ** 2  # cm^2 (sphere radius = 1.3 cm)
betaVal = 0.5
amu = 16
Ti = 1000 / 11600
Isat_i = currentIonSaturation(ni, Ti, dV, Aeff, betaVal, amu)
