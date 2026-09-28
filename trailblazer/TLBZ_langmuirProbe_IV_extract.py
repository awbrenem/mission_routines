"""
Analyze the artificially-generated Langmuir Probe sweep for Trailblazer (TBLZ_langmuirProbe_analysis.py) to see if the
values of ne, ni, Te, Vp I determine here match those I input into the artificially-generated sweep.





e.g. Barjatya thesis Chapter 3.3.2

Step 1 (rough Vp and rough ni):
1. [Vf accurate] Vf where total current crosses zero
2. [Ie] Fit a line in ion saturation region and subtract this from total current --> gives us the electron collection current (+I values on I-V curve only)
3. [Vp rough (aka plasma potential)] Ideally, this is the inflection point where the I-V curve deviates from an exponential.
    Since this can be difficult to see, we can take first derivative of e- current wrt voltage and look for the maxima.
    NOTE: that if Te < 5000 K we can limit the search for this maxima to within 0.5 eV of Vf.
    NOTE: the dI/dV curve tends to get very noisy at around the inflection point, making it difficult to identify.
4. [ni rough] Ion density (first order approx.) determined by equating the value of the ion saturation current at Vp to the ion ram current.
    Isat(Vp) = Iram

Step 2:
1. [Vp accurate (aka plasma potential)] Do linear least squares fit of the total collected current to Eqn 3.7 (Barjatya thesis).
    Use rough ni from Step 1 as for Eqn 3.7 and use the least squares fit to solve for Te and Vp.
    Do the fit only b/t values Vf-0.35 eV and Vf+0.5 eV. For larger + values the e- current will deviate from exponential,
    and for larger - values the ion current will dominate.
2. [ni accurate] Determine by evaluating ion saturation current at Vp accurate.
3. [Te accurate] Determine by evaluating the slope of the log of the e- current for V<Vp. I.e. Te ~ 1/slope of ln(I) vs V for V<Vp in electron retardation region

------------------------------


"""

import pickle
import sys

import numpy as np
from matplotlib import pyplot as plt
from scipy.io import readsav

sys.path.append("/Users/abrenema/Desktop/code/Aaron/github/antenna_analysis/")
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

k = 1.38e-23  # Boltzmann constant (m^2 kg s^-2 K^-1)
q = 1.6e-19  # elementary charge (C)
me = 9.11e-31  # electron mass (kg)
mi = 1.67e-27  # proton mass (kg)
eo = 8.85e-12  # permittivity of free space (s^4 A^2 m^-3 kg^-1)


# Units check out. NOTE: curve_fit always expects independent variable to be first argument!
def iv_model(V, Te_eV, Vp):
    # ion ram current (negative sign because ion current flows TO the probe)
    Iram_uA = -ni_rough * q * AreaProbe * VelSC * 1e6

    # thermal electron current amplitude in microamps
    Te_K = Te_eV * 11604.5
    Ithe_uA = ni_rough * q * AreaProbe * np.sqrt(k * Te_K / (2.0 * np.pi * me)) * 1e6

    return Iram_uA + Ithe_uA * np.exp(q * (V - Vp) / (k * Te_K))


# Find values of ne and beta that allow the e- saturation current to match the observed data
# NOTE: curve_fit always expects independent variable to be first argument
def iv_model_electronSaturation(V, ne, betaOLM):

    # thermal electron current amplitude in microamps
    Te_K = Te_fit * 11604.5
    Ithe_uA = ne * q * AreaProbe * np.sqrt(k * Te_K / (2.0 * np.pi * me)) * 1e6

    return Ithe_uA * (1 + q * (V - Vp_fit) / (k * Te_K)) ** betaOLM


def segments_from_mask(mask, time):
    change = np.diff(mask.astype(int))
    starts = np.where(change == 1)[0] + 1
    ends = np.where(change == -1)[0]
    if mask[0]:
        starts = np.concatenate(([0], starts))
    if mask[-1]:
        ends = np.concatenate((ends, [len(mask) - 1]))
    return [(time[s], time[e], s, e) for s, e in zip(starts, ends)]


def select_positive_plateaus(periods, voltage, voltage_min=1.5, max_samples=12000):
    selected = []
    for t0, t1, s, e in periods:
        mean_v = np.mean(voltage[s : e + 1])
        n = e - s + 1
        if mean_v > voltage_min and n < max_samples:
            selected.append((t0, t1, s, e, mean_v, n))
    return selected


# ---------------------------------------------------------------------------


with open(
    "/Users/abrenema/Desktop/code/Aaron/github/mission_routines/trailblazer/TBLZ_LangmuirProbe_Itot.pkl",
    "rb",
) as f:
    loaded_arr = pickle.load(f)

iz = loaded_arr[0] / 1000
vz = loaded_arr[1]


plt.plot(vz, iz)


VelSC = 7600  # m/s


##upsweep
# ss = 336.876
# se = 338.12
##downsweep
# ss = 337.627
# se = 337.876


# plot of an I-V curve for the selected sweep
fig, axs = plt.subplots(2, figsize=(10, 8))
axs[0].plot(vz, iz)
axs[0].set_yscale("linear")
axs[0].set_xlabel("Sweep Voltage (V)")
axs[0].set_ylabel("Sweep Current (micro Amps)")
axs[0].set_title("Langmuir Probe I-V Curve")
axs[0].axhline(y=0, color="red", linestyle="--", linewidth=0.5, label="Threshold")
axs[1].plot(vz, iz)
axs[1].set_yscale("linear")
axs[1].set_xlabel("Sweep Voltage (V)")
axs[1].set_ylabel("Sweep Current (micro Amps)")
axs[1].set_title("Langmuir Probe I-V Curve")
axs[0].set_ylim(-1.2, 60)
for i in range(2):
    axs[1].set_xlim(-1.5, 3)
plt.show()


# ---------------------------------------------------
# Step 1 [Barjatya thesis]
# ---------------------------------------------------

# ******
# Part 1 (Vf) - find floating potential (where current crosses zero)
# This will be the lower bound of the electron retardation region
# ******

goo = np.argmin(np.abs(iz))
Vf = vz[goo]

# tmp = np.where(iz > 0)[0]
# Vf = vz[tmp[0]]  # floating potential in volts


# ******
# Part 2 - fit a line to the ion saturation region and subtract this from the total current to get the electron collection current
# ******

# We can fit the line to the region where V < Vf - 0.5 eV (i.e. where the current is dominated by ions).
Vf_offset = 0.5  # eV. Change this value to get a better fit on the ion saturation region. Good guess is 0.5 eV

tmp = np.where(vz < (Vf - Vf_offset))[0]
# Fit a line to this region
p = np.polyfit(vz[tmp], iz[tmp], 1)  # linear fit (y = mx + b) to ion saturation region
iz_fit = np.polyval(p, vz)  # evaluate the fit at all voltages
iz_electron = iz - iz_fit  # subtract the ion fit from the total current to get


# Check goodness of fit for the ion saturation region
fig, axs = plt.subplots(3, figsize=(10, 8))
axs[0].plot(vz, iz)
axs[0].plot(vz, iz_fit, linestyle="--")
axs[0].set_yscale("linear")
axs[0].set_ylim(-10, 60)
axs[0].set_xlabel("Sweep Voltage (V)")
axs[0].set_ylabel("Sweep Current (micro Amps)")
axs[0].set_title("Jets I-V Curve goodness of ion sat region fit")
axs[0].scatter(Vf, 0, color="red", label="Floating Potential")
axs[0].axvline(
    x=Vf, color="red", linestyle="--", linewidth=0.5, label="Floating Potential"
)
axs[0].axhline(y=0, color="red", linestyle="--", linewidth=0.5, label="Threshold")
axs[1].plot(vz, iz)
axs[1].plot(vz, iz_fit, linestyle="--")
axs[1].set_yscale("linear")
axs[1].set_xlabel("Sweep Voltage (V)")
axs[1].set_ylabel("Sweep Current (micro Amps)")
axs[1].set_ylim(-0.8, 0)
axs[2].plot(vz, iz_electron)
axs[2].axhline(y=0, color="red", linestyle="--", linewidth=0.5, label="Threshold")
axs[2].axvline(
    x=Vf, color="red", linestyle="--", linewidth=0.5, label="Floating Potential"
)
axs[2].set_yscale("linear")
axs[2].set_ylim(-2, 60)
axs[2].set_xlabel("Sweep Voltage (V)")
axs[2].set_ylabel("Sweep Current (micro Amps)")
axs[2].set_title("Jets I-V Curve w/ Ion Current Removed")


# ******
# Part 3 (Vp_rough) - find rough plasma potential by taking the first derivative of the electron current
# wrt voltage and finding the location of the maxima.
#    Note that if Te < 5000 K we can limit the search for this maxima to within 0.5 eV of Vf.
# ******

dIedV = np.gradient(iz_electron, vz)

# Find the location of the maxima in the smoothed derivative.
# NOTE: for Te<5000 K this should be within 0.5 eV of Vf
# goodv = np.where((vz > (Vf - 0.5)) & (vz < (Vf + 0.5)))[0]
goodv = np.where((vz > (Vf - 0.5)) & (vz < (Vf + 1)))[0]
if len(goodv) == 0:
    raise ValueError(
        "No finite points found in the search window around Vf for Vp_rough."
    )
max_idx = np.argmax(dIedV[goodv])
Vp_rough = vz[goodv][max_idx]  # rough plasma potential in volts


fig, axs = plt.subplots(3, figsize=(10, 8))
axs[0].plot(vz, iz_electron)
axs[0].set_yscale("linear")
axs[0].set_xlabel("Sweep Voltage (V)")
axs[0].set_ylabel("Electron current (micro Amps)")
axs[0].set_title("Electron current after ion subtraction")
axs[0].axhline(y=0, color="red", linestyle="--", linewidth=0.5, label="Threshold")
axs[1].plot(vz, dIedV)
axs[1].set_yscale("linear")
axs[1].set_xlabel("Sweep Voltage (V)")
axs[1].set_ylabel("dIe/dV")
axs[1].set_title("Smoothed dIe/dV")
axs[1].set_ylim(-0.1, 120)
axs[2].plot(vz, np.log(np.abs(iz_electron)))
axs[2].set_yscale("linear")
axs[2].set_xlabel("Sweep Voltage (V)")
axs[2].set_ylabel("ln(|I|)")
axs[2].set_title("Log electron current")
axs[2].set_ylim(-10, 4)
for i in range(3):
    axs[i].set_xlim(-2, 3)
    axs[i].axvline(
        x=Vf, color="red", linestyle="--", linewidth=0.5, label="Floating Potential"
    )
    axs[i].axvline(
        x=Vp_rough, color="red", linestyle="--", linewidth=0.5, label="Vp rough"
    )


# ******
# Part 4 (ni_rough) - find rough ion density by equating the value of the ion saturation current ("Iz_fit") at Vp to the ion ram current.
# Note that this really only works when heavy ions dominate and we can thus ignore thermal motions.
# ******

Iz_fit_Vp = np.polyval(p, Vp_rough)  # value of the ion saturation current at Vp
e = 1.6e-19  # elementary charge in Coulombs

r = 0.32  # cm (cylindrical probe radius)
l = 5  # cm (probe length)
A = 2 * np.pi * r * l / (1e4)  # m^2 (probe area)
AreaProbe = A * 0.5  # sunlit area only


ni_rough = np.abs(
    (Iz_fit_Vp / 1e6) / (e * AreaProbe * VelSC)
)  # ion density in m^-3, using the formula Iram = ni * e * A * sqrt(2*m_i*Vp) where A is the probe area

print(ni_rough / 1e6)  # cm-3


# ---------------------------------------------------
# Step 2 [Barjatya thesis]
# ---------------------------------------------------

# ******
# Part 1 [Te, Vp] Do linear least squares fit of the total collected current to Eqn 3.7 (Barjatya thesis).
# Use rough ni from Step 1 as for Eqn 3.7 and use the least squares fit to solve for Te and Vp.
# Do the fit only b/t values Vf-0.35 eV and Vf+0.08 eV. For larger + values the e- current will deviate from exponential,
# and for larger - values the ion current will dominate. Note that the resulting value of Vp will likely be closer to 2 V.

# NOTE: I'm using scipy.optimize.curve_fit to do this fit, since we have two free parameters (Te and Vp) and the equation is non-linear.
# We want to compare "iz" to "Iram" + "Ithe * Ie_ret" where "Iram" is the ion ram current, "Ithe" is the thermal electron current, and "Ie_ret" is the exponential term that accounts for the retarding potential.


# Need to determine voltage range for the fit.
# Barjatya uses Vf - 0.35 eV to Vf + 0.08 eV
# ******


# goodvRet = np.where((vz > (Vf - 0.35)) & (vz < (Vf + 0.08)))[0]
# goodvRet = np.where((vz > (Vf)) & (vz < (Vf + 0.3)))[0]
goodvRet = np.where((vz > (-0.75)) & (vz < (-0.5)))[0]
# goodvRet = np.where((vz > (Vf - 0.35)) & (vz < (Vf + 0.5)))[0]
# goodvRet = np.where((vz > (Vf - 0.35)) & (vz < Vp_rough))[0]
# goodvRet = np.where((vz > Vf) & (vz < Vp_rough))[0]
# goodvRet = np.where((vz > Vf) & (vz < 1.4))[0]


# Temperature bounds
# Upper bound < 5000 K (0.43 eV) to avoid the electron current deviating from exponential


# tst = iv_model(1.5, 0.1, Vp_rough)
init_guess = [0.2, Vp_rough]  # [Te, Vp]
# bounds = ([0.01, Vf - 0.5], [2, Vf + 0.8])
# bounds = ([0.01, Vf - 0.5], [2, Vp_rough + 0.1])
bounds = ([0.01, Vf], [0.43, Vp_rough + 0.2])


# init_guess = [0.1, Vf + 0.2]
# bounds = ([0.01, Vf], [0.43, 2])


print(init_guess)  # [Te, Vp]
print(bounds)

# Put current in terms of microAmps b/c this is what the output of curve_fit will have
# Curve_fit will adjust the parameters of iv_model to minimize the difference between the model and the data in a least-squares sense.
# popt, pcov = np.polyfit(iv_model, vz[goodvRet], iz[goodvRet], 3)
# , p0=init_guess, bounds=bounds

popt, pcov = curve_fit(
    iv_model, vz[goodvRet], iz[goodvRet], p0=init_guess, bounds=bounds
)


Te_fit, Vp_fit = popt
Te_err, Vp_err = np.sqrt(np.diag(pcov))

print(f"Te fit = {Te_fit:.3f} eV ± {Te_err:.3f} eV")
print(f"Vp fit = {Vp_fit:.4f} V ± {Vp_err:.4f} V")


# To see how fitting went, plot the log value. The Te should be fit to a region with a straight line.
plt.plot(vz, np.log(np.abs(iz)), ".", label="data")
plt.scatter(Vf, 0, color="red", label="Floating Potential")
plt.scatter(
    Vp_fit,
    np.log(np.abs(iv_model(Vp_fit, *popt))),
    color="blue",
    label="Fitted Plasma Potential",
)
plt.axvline(x=Vf, color="red", linestyle="--", linewidth=0.5, label="Vf")
plt.axvline(
    x=Vp_rough,
    color="red",
    linestyle="--",
    linewidth=0.5,
    label="Rough Plasma Potential",
)
plt.plot(vz[goodvRet], np.log(np.abs(iv_model(vz[goodvRet], *popt))), "-", label="fit")
plt.xlabel("Sweep Voltage (V)")
plt.ylabel("Sweep Current (ln(|micro Amps|))")
plt.axhline(y=0, color="red", linestyle="--", linewidth=0.5)
plt.title("Langmuir Probe I-V Fit")
plt.xlim(-2, 2)
plt.yscale("linear")
plt.ylim(-10, 3)
plt.legend()
plt.show()


"""
plt.plot(vz, iz, ".", label="data")
plt.scatter(Vf, 0, color="red", label="Floating Potential")
plt.scatter(
    Vp_fit, iv_model(Vp_fit, *popt), color="blue", label="Fitted Plasma Potential"
)
plt.axvline(x=Vf, color="red", linestyle="--", linewidth=0.5, label="Vf")
plt.axvline(
    x=Vp_rough,
    color="red",
    linestyle="--",
    linewidth=0.5,
    label="Rough Plasma Potential",
)
plt.plot(vz[goodvRet], iv_model(vz[goodvRet], *popt), "-", label="fit")
plt.xlabel("Sweep Voltage (V)")
plt.ylabel("Sweep Current (micro Amps)")
plt.axhline(y=0, color="red", linestyle="--", linewidth=0.5)
plt.title("Langmuir Probe I-V Fit")
plt.xlim(-2, 0)
plt.yscale("log")
plt.ylim(0.1, 100)
plt.legend()
plt.show()
"""

# ********************************************
# Part 2 [ni] - refine ni by evaluating ion saturation current line fit at accurate Vp
# ********************************************

Iz_fit_Vp = np.polyval(p, Vp_fit)  # value of the ion saturation current at Vp
ni = np.abs(
    Iz_fit_Vp / 1e6 / (e * AreaProbe * VelSC)
)  # ion density in m^-3, using the formula Iram = ni * e * A * sqrt(2*m_i*Vp) where A is the probe area (assumed to be 1 cm^2 here)


# ********************************************
# Part 3 [ne, beta] - Determine ne and check accuracy of fit values to OLM theory
# ********************************************


######################################
#############################################################
# MANUALLY TEST FIT PARAMETERS TO ELECTRON SATURATION DATA

# ni and beta=1.2 match almost perfectly

# Electron saturation region
goodvEsat = np.where(vz > Vp_fit)[0]

Ie_tst = np.empty(len(vz[goodvEsat]))

for i in range(len(goodvEsat)):
    # Ie_tst[i] = iv_model_electronSaturation(vz[goodvEsat][i], 26738 * 1e6, 1.5)
    Ie_tst[i] = iv_model_electronSaturation(vz[goodvEsat][i], 0.6 * ni / 1e6, 1.5)


plt.plot(vz, iz)
plt.plot(vz[goodvEsat], Ie_tst)
plt.yscale("linear")
plt.xlabel("Sweep Voltage (V)")
plt.ylabel("Sweep Current (micro Amps)")
plt.title("Langmuir Probe I-V Curve")
plt.axvline(x=Vf, color="blue", linestyle="--", linewidth=0.5, label="Vf")
plt.axvline(x=Vp_fit, color="blue", linestyle="--", linewidth=0.5, label="Vp fit")
plt.xlim(0, 3)
plt.ylim(-0.5, 20)
plt.show()

#############################################################


# put density in log values, otherwise it dominates and messes up the fit.
init_guess = [np.log10(ni), 0.5]
bounds = ([np.log10(ni) - 0.01, 0], [np.log10(ni) + 0.01, 1.5])

print(init_guess)  # [ne, betaOLM]
print(bounds)


# Put current in terms of microAmps b/c this is what the output of curve_fit will have
# popt, pcov = curve_fit(iv_model, vz[goodv3], iz[goodv3] * 1e6, p0=init_guess, bounds=bounds)
popt, pcov = curve_fit(
    iv_model, vz[goodvEsat], iz[goodvEsat], p0=init_guess, bounds=bounds
)


ne_fitLog, beta_fit = popt
ne_errLog, beta_err = np.sqrt(np.diag(pcov))

ne_fit = 10**ne_fitLog
ne_err = 10**ne_errLog
print(f"ne fit = {ne_fit / 1e6:.3f} cm-3 ± {ne_err / 1e6:.3f} cm-3")
print(f"Beta fit = {beta_fit:.4f}  ± {beta_err:.4f} ")


plt.plot(vz, iz, ".", label="data")
plt.axvline(x=vz[goodvEsat][0], color="red", linestyle="--", linewidth=2)
plt.axvline(x=vz[goodvEsat][-1], color="red", linestyle="--", linewidth=2)
# plt.plot(vz[goodvEsat], iv_model_electronSaturation(vz[goodvEsat], *popt), '-', label='fit')
plt.plot(
    vz[goodvEsat],
    iv_model_electronSaturation(vz[goodvEsat], ne_fit / 1e6, beta_fit),
    "-",
    label="best fit with beta=" + str(round(beta_fit, 2)),
)
plt.plot(
    vz[goodvEsat],
    iv_model_electronSaturation(vz[goodvEsat], ne_fit / 1e6, 0.5),
    "--",
    label="OLM, beta=0.5",
)
plt.plot(
    vz[goodvEsat],
    iv_model_electronSaturation(vz[goodvEsat], ne_fit / 1e6, 1),
    "--",
    label="OLM, beta=1",
)
plt.xlabel("Sweep Voltage (V)")
plt.ylabel("Sweep Current (micro Amps)")
plt.title("Langmuir Probe I-V Fit")
plt.legend()
plt.show()


print("h")
print("h")
