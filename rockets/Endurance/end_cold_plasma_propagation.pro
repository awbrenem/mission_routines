;Cold plasma dispersion for Endurance rocket 
;
.compile /Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/cold_plasma_dispersion.pro


;t=840 ;sec 

epol = 0.77
freq = 5000.
dens = 246920.0
Bo = 47669.  
pH = 0.02
pO = 0.98
pHe = 0.



x = cold_plasma_dispersion(epol=epol,freq=freq,dens=dens,Bo=Bo,H_plus=pH,He_plus=pHe,O_plus=pO)

  
  R               FLOAT           2961.94
   L               FLOAT          -3008.46
   D               FLOAT           2985.20
   S               FLOAT          -23.2593
   P               FLOAT          -796503.
   P_S             FLOAT           34244.5
   RLMINUSPS       FLOAT      -2.74370e+07
   A               FLOAT              -NaN
   B               FLOAT              -NaN
   C               FLOAT       7.09755e+12
   F               FLOAT              -NaN
   RL_S            FLOAT           383111.
   R_CUTOFF        FLOAT       5.18029e+06
   R_CUTOFF_ELECTRONS_ONLY
                   FLOAT       5.17925e+06
   L_CUTOFF        FLOAT       3.84529e+06
   L_CUTOFF_ELECTRONS_ONLY
                   FLOAT       3.84452e+06
   P_CUTOFF        FLOAT       4.46236e+06
   P_CUTOFF_ELECTRONS_ONLY
                   FLOAT       4.46226e+06
   X_CUTOFF        FLOAT       5.17925e+06
   Z_CUTOFF        FLOAT       3.84452e+06
   S_RESONANCE     FLOAT       4.65760e+06
   R_RESONANCE     FLOAT       1.33473e+06
   Z_OBLIQUE_RESONANCE
                   FLOAT              -NaN
   FUH_RESONANCE   FLOAT       4.65760e+06
   
   
     
;  freq in Hz
;  kmag in 1/km
;  theta_kb in deg
;  n=index of refraction
;  wavelength in km
;  resangle in deg
;  phasevel in km/sec
;  bpol is magnetic polarization ratio
;  energy in eV for nonrelativistic energies
;  resonance and cutoff values in Hz, calculated without ion contributions
;  cBw_to_Ew is the ratio cB/E


 