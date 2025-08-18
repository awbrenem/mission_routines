#Estimate the cyclotron/Landau damping energy for Bernstein waves based on the Doppler-shifted cyclotron resonance condition. 
#Makes nice contour plots of Ez vs f and the perpendicular wavelength. 




#t=840 sec is when WHAMP was run. 
#Primary wave power is at 5.5, 6.3, 7.0 kHz, with the 6.3 kHz band being the largest
#fcH = 754 Hz
#fcH*6 = 4525
#fcH*7 = 5279
#fcH*8 = 6033
#fcH*9 = 6787 
#fcH*10 = 7541


import numpy as np
import matplotlib.pyplot as plt

#WHAMP perpendicular wavelength estimated to be 0.5 - 2 meters
#and phase velocities from 3 - 10 km/sec 
#Solutions exist for theta_kb > 85 deg only. 





# Parameter ranges
f_vals = np.linspace(4000, 7500, 100)           # Hz
lambda_perp_vals = np.linspace(0.5, 2, 100)     # meters

F, LAMBDA = np.meshgrid(f_vals, lambda_perp_vals)


mp = 1.27e-27   #proton mass (kg)
m_ion = 1  #multiplicative factor for ion mass (e.g., 1 for protons, 16 for O+)
e1eV = 1.6e-19  #Joules/eV 


theta_kb = 89.95   #[-->increasing increases Ez]
#theta_kb = 89   #[-->increasing increases Ez]
#lambda_perp = 2  #meters  (0.5-2m predicted from WHAMP)  [-->increasing increases Ez]
#f = 7500. #wave freq
fc = 754. #proton cyclotron freq

max_energy = 1e6  #eV (max energy for cyclotron damping color scale to plot - this energy represents the maximum color)

harmonics = [6, 7, 8, 9]
fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
axes = axes.flatten()

for idx, n in enumerate(harmonics):
    w = 2 * np.pi * F
    wc = 2 * np.pi * fc
    kperp = 2 * np.pi / LAMBDA
    kpar = kperp * np.cos(theta_kb * np.pi / 180)
    Vz = (w - n * wc) / kpar
    #Ez = 0.5 * mp * Vz**2 / e1eV
    Ez = 0.5 * m_ion * mp * Vz**2 / e1eV

    #Ez_masked = np.ma.masked_where((Ez <= 0) | (Ez >= max_energy), Ez)
    Ez_masked = np.ma.masked_where((Ez <= 0) | (Ez >= 1e8), Ez)

    contour = axes[idx].contourf(
        F, LAMBDA, np.log10(Ez_masked), levels=50, cmap='viridis', vmin=-1.5, vmax=np.log10(max_energy)
    )
    # Add contour lines for 100 eV and 500 eV
    cs = axes[idx].contour(
        F, LAMBDA, Ez, levels=[100, 1000, 10000], colors=['white', 'black', 'black'], linewidths=1.5
    )
    axes[idx].clabel(cs, fmt={100: '100 eV', 1000: '1000 eV', 10000:'10 keV'}, colors=['white', 'black', 'black', 'black'], fontsize=10)

    axes[idx].set_xlabel('f (Hz)')
    axes[idx].set_ylabel('lambda_perp (m)')
    axes[idx].set_title(
        f'Harmonic={n}, theta_kb={theta_kb}°, fcH+={fc} Hz'
    )
    axes[idx].vlines([5500, 6300, 7000], ymin=0.5, ymax=2, colors='red', linestyles='dashed')
    axes[idx].legend(['5.5 kHz', '6.3 kHz', '7.0 kHz'], loc='upper left', bbox_to_anchor=(1.02, 1))

# Add a single colorbar for all subplots
fig.colorbar(contour, ax=axes, orientation='vertical', label='log10(Ez) (eV)')
fig.suptitle('Cyclotron damping @ 840 sec: log10(Ez) vs f and lambda_perp\nRed lines = Bernstein power peaks\nmax energy plotted = '+str(max_energy)+' eV', fontsize=16)
plt.show()


















"""
harmonics = [6, 7, 8, 9]
#harmonics = [0]
fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
axes = axes.flatten()

for idx, n in enumerate(harmonics):
    w = 2 * np.pi * F
    wc = 2 * np.pi * fc
    kperp = 2 * np.pi / LAMBDA
    kpar = kperp * np.cos(theta_kb * np.pi / 180)
    Vz = (w - n * wc) / kpar
    Ez = 0.5 * mp * Vz**2 / e1eV

    Ez_masked = np.ma.masked_where((Ez <= 0) | (Ez >= max_energy), Ez)

    contour = axes[idx].contourf(F, LAMBDA, np.log10(Ez_masked), levels=50, cmap='viridis', vmin=-1.5, vmax=np.log10(max_energy))
    axes[idx].set_xlabel('f (Hz)')
    axes[idx].set_ylabel('lambda_perp (m)')
    axes[idx].set_title(
        f'Harmonic={n}, theta_kb={theta_kb}°, fcH+={fc} Hz'
    )
    axes[idx].vlines([5500, 6300, 7000], ymin=0.5, ymax=2, colors='red', linestyles='dashed')
    axes[idx].legend(['5.5 kHz', '6.3 kHz', '7.0 kHz'], loc='upper left', bbox_to_anchor=(1.02, 1))

# Add a single colorbar for all subplots
fig.colorbar(contour, ax=axes, orientation='vertical', label='log10(Ez) (eV)')
fig.suptitle('Cyclotron damping @ 840 sec: log10(Ez) vs f and lambda_perp\nRed lines = Bernstein power peaks\nmax energy plotted = '+str(max_energy)+' eV', fontsize=16)
plt.show()
"""
