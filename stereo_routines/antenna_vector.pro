 pro antenna_vector,sc,antV;   Version = '02 July 2017';   ;  file: \laptop documents;	calculates effective length vector using Oswald tables, rheometry
;	the vectors are in spaceccraft coordinates, approx RTN coordinates
;	except they are in the order R +-N T
;	On 30 June 2017, values were changed to match Oswald et al, Adv. Space;	Res. 43,355, (2009);	X direction (SC, not antenna)  is toward sun, so all antennas have 
;	negative X component. In terms of coords RTN, approximately +Y (S/C) is -N;	 for A, N for B, approx +Z (S/C) is T for B, -T for A;	the first index is the spacecraft X,Y or Z direction, 
;	the second is 0=Z, 1=Y, 2=X antenna. 0=Z corresponds to Oswald original choice,;	 but ours is 0 = X, 2 = z;	theta is the angle to the +X axis (toward sun), phi the angle in 
;		the YZ plane, measured from -Z direction;;	The drawing, e.g. on p 102 of N.B. V 8, is looking toward the sun, so;		phi is positie counterclockwise (Oswald 2009 just below fig 5);; S.D. Bale ·et al, The Electric Antennas for the STEREO/WAVES Experiment ;Space Sci Rev (2008) 136: 529–547  DOI 10.1007/s11214-007-9251-x;;W.Macher, T H Oswald, G Fischer and H O Rucker;Meas. Sci. Technol. 18 (2007) 3731–3742 doi:10.1088/0957-0233/18/12/008;Rheometry of multi-port spaceborne antennas including mutual antenna;capacitances and application toSTEREO/WAVES; 
antV = fltarr(3,3) 
degrad = !pi/180.
;;  print,'degree to radians,inv ',degrad
;  print,'!radeg ',!radeg;  print,'spacecraft ',sc;;	Now same lengths for A and B;	Now (7 Jan 2007) changed to "final report' and average HGa;	and open input;	Now changed to Oswald et al Adv. Space Res. 2009 (30 Jun 2017);	note Z antenna is Oswald's 1;;  print,'got to here ';    phirot = 0.       ; phirot was originally to account for position;							in flight but now use attitude matrices   
 ; 	X antenna       theta = 120.5       phi   = -135.5 + phirot
     antV(0,0) = cos(degrad*theta); radial (S/C X) component of Z antenna
     antV(1,0) = sin(degrad*theta)*sin(degrad*phi)
     antV(2,0) =-sin(degrad*theta)*cos(degrad*phi);	y antenna     theta = 114.9     phi = 127.2 + phirot
    antV(0,1) = cos(degrad*theta)	    ; radial component of Y antenna
     antV(1,1) = sin(degrad*theta)*sin(degrad*phi)
     antV(2,1) =-sin(degrad*theta)*cos(degrad*phi)
; Z antenna
     theta = 125.2     phi = 15.9 + phirot
      antV(0,2) = cos(degrad*theta)	     ; radial component of X antenna
     antV(1,2) = sin(degrad*theta)*sin(degrad*phi)
     antV(2,2) =-sin(degrad*theta)*cos(degrad*phi);  
return
end