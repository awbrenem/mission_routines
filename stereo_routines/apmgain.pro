pro APMGAIN, F,GAPMCPI,GAPMPI,PHASE
;
;	GAIN OF ANTENNA POTENTIAL MONITOR CIRCUITS, INCLUDING BOTH PREAMP
;		BOX AND MAIN BOX
;	 ANTENNA IMPEDANCE NOT TAKEN INTO ACCOUNT
;
;       IN THE FOLLOWING:
;               R2, CIN ARE THE RESISTANCE FROM DEVICE + INPUT TO
;                       CHASSIS AND CCP IS COUPLING CAPACITOR
;               R1 IS RESISTANCE FROM ANTENNA TO S/C GROUND
;
;	DATA R1,R2,R3,R4,C4 /1.5E8,33.E6,33.E6,150.E6,0.E-10/

	R1=1.5E8
	R2=33.E6
	R3=33.E6
	R4=150.E6
	C4=0.
;	DATA CANT /53.E-12/		! from Bob Manning, for Stacer
;	DATA CB /50.E-12/
;	DATA CC /240.E-12/	;! 8 FT., 30 PF/FT
;	DATA CC /150.E-12/	;! 5 FT., 30 PF/FT
	CCABLE = 134.7E-12		;! average of measured values
	C3 = 100.E-12

	R5 = 39.2E3
	R6 = 1.E5
	C5 = .1E-6
	R7 = 6.81E3
	C7 = 1.E-7 
        TWOPI = 6.28318531
;
        W = TWOPI*F
;
;       GAINS
;
        Y1 = COMPLEX((1./R3+1./R4), W*(C3+CCABLE+C4))
	GAPMCPI = 1./(1. + R2*Y1)
;
;	CIRCUIT IN MAIN BOX
;
	GAPMCPI = GAPMCPI/COMPLEX(R6/R5,W*C5*R6)
	GAPMCPI = -GAPMCPI/COMPLEX(1.,W*R7*C7)
;	apm inverts signa, but above makes + in be + out
        GAPMPI = ABS(GAPMCPI)
	PHASE = 57.2957795*ATAN(IMAGINARY(GAPMCPI),FLOAT(GAPMCPI))
;
	RETURN
	END
