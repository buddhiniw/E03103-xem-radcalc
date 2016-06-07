      SUBROUTINE w1w2(qsq,w,w1,w2,inp)
C...
C...THIS FITS ARE FOR HYDROGEN AND DEUTERIUM.
C...USES ATWOOD'S RESONANCE FITS AND CALCULATE VW2 AS:
C...VW2=B*F2 WHERE B IS BACKGROUND WITH RESONANCES AND F2 IS
C...A UNIVERSAL FUNCTION TO THE DEEP INELASTIC
C... 
C...VW2,W2,W1 AND XSECTN ARE THE FITTED STRUCTURE
C...FUNCTIONS AND THE FITTED CROSS SECTION RESPECTEVILY. (OUTPUT)
C...XMOTT IS THE MOTT CROSS SECTION                      (OUTPUT)
C...X IS THE PARTON MOMENTUM FRACTION, X = QSQ/2MV       (OUTPUT)
C...
      REAL*8 C(24),CF(11),cn(11),CD(24),CFD(11),F2,B,QSQ,WW
      DATA PMSQ/.880324/,TPM/1.876512/
      DATA R/.18/,ALPHAX/137.0388/
C...
C...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...CROSS SECTION FIT FOR HYDROGEN OR DEUTERIUM
C...USES ATWOOD'S RESONANCE FIT WITH COEFFICIENTS FROM FITTING
C...E87,E49A,E49B.COEFFICIENTS SUPLIED BY ARIE BODEK FOR E139
C...FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE
C...CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59
C...THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)
C...
C...C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE
C... TERMS)
C...
      DATA   C(1) / 0.10741163D 01/,  C(2) / 0.75531124D 00/,
     *       C(3) / 0.33506491D 01/,  C(4) / 0.17447015D 01/,
     *       C(5) / 0.35102405D 01/,  C(6) / 0.10400040D 01/,
     *       C(7) / 0.12299128D 01/,  C(8) / 0.10625394D 00/,
     *       C(9) / 0.48132786D 00/,  C(10)/ 0.15101467D 01/,
     *       C(11)/ 0.81661975D-01/,  C(12)/ 0.65587179D 00/,
     *       C(13)/ 0.17176216D 01/,  C(14)/ 0.12551987D 00/,
     *       C(15)/ 0.74733793D 00/,  C(16)/ 0.19538129D 01/,
     *       C(17)/ 0.19891522D 00/,  C(18)/-0.17498537D 00/,
     *       C(19)/ 0.96701919D-02/,  C(20)/-0.35256748D-01/,
     *       C(21)/ 0.35185207D 01/,  C(22)/-0.59993696D 00/,
     *       C(23)/ 0.47615828D 01/,  C(24)/ 0.41167589D 00/
C...
C...CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)
C...OMEGAW FIT
C...
      DATA CF(1) / 0.25615498D 00/,  CF(2) / 0.21784826D 01/,
     *     CF(3) / 0.89783738D 00/,  CF(4) /-0.67162450D 01/,
     *     CF(5) / 0.37557472D 01/,  CF(6) / 0.16421119D 01/,
     *     CF(7) / 0.37635747D 00/,  CF(8) / 0.93825625D 00/,
     *     CF(9) / 0.10000000D 01/,  CF(10)/ 0.0           /,
     *     CF(11)/ 0.50000000D 00/
      DATA CN(1) / 0.06400000D 00/,  CN(2) / 0.22500000D 00/,
     *     CN(3) / 0.41060000D 01/,  CN(4) /-0.70790000D 01/,
     *     CN(5) / 0.30550000D 01/,  CN(6) / 0.16421119D 01/,
     *     CN(7) / 0.37635747D 00/,  CN(8) / 0.93825625D 00/,
     *     CN(9) / 0.10000000D 01/,  CN(10)/ 0.0           /,
     *     CN(11)/ 0.50000000D 00/
C...
C...
C... FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE
C... CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96
C... THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)
C...
C... CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONAN
C... TERMS)
C...
      DATA  CD(1) / 0.10521935D 01/, CD(2) / 0.76111537D 00/,
     *      CD(3) / 0.41469897D 01/, CD(4) / 0.14218146D 01/,
     *      CD(5) / 0.37119053D 01/, CD(6) / 0.74847487D 00/,
     *      CD(7) / 0.12399742D 01/, CD(8) / 0.12114898D 00/,
     *      CD(9) / 0.11497852D-01/, CD(10)/ 0.14772317D 01/,
     *      CD(11)/ 0.69579815D-02/, CD(12)/ 0.12662466D 00/,
     *      CD(13)/ 0.15233427D 01/, CD(14)/ 0.84094736D-01/,
     *      CD(15)/ 0.74733793D 00/, CD(16)/ 0.19538129D 01/,
     *      CD(17)/ 0.19891522D 00/, CD(18)/-0.24480414D 00/,
     *      CD(19)/ 0.14502846D-01/, CD(20)/-0.35256748D-01/,
     *      CD(21)/ 0.35185207D 01/, CD(22)/-0.21261862D 00/,
     *      CD(23)/ 0.69690531D 01/, CD(24)/ 0.40314293D 00/
C...
C...CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)
C....OMEGAW FIT
C...
      DATA CFD(1) / 0.47708776D 00/, CFD(2) / 0.21601918D 01/,
     *     CFD(3) / 0.36273894D 01/, CFD(4) /-0.10470367D 02/,
     *     CFD(5) / 0.49271691D 01/, CFD(6) / 0.15120763D 01/,
     *     CFD(7) / 0.35114723D 00/, CFD(8) / 0.93825625D 00/,
     *     CFD(9) / 0.10000000D 01/, CFD(10)/ 0.0           /,
     *     CFD(11)/ 0.50000000D 00/
C...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...
C...COMPUTE SOME KINEMATIC QUANTITIES
C...
	ww = w
	v = (w*w + qsq - pmsq)*.5326
	vsq = v*v
	bres = b(ww,qsq,c)
!	bres =1.
	if (inp .eq. 1) then
	univ =f2(ww,qsq,cf)
	else
	univ = f2(ww,qsq,cn)
	endif
c...compute vw2,w2,w1,and the cross section
C...
      !write(6,*) 'In w1w2 with w and Q', w, qsq
      VW2=UNIV*bres	!to take out the background and resonance
      W2=VW2/V		! set bres to 1.0
      W1=(1.0+VSQ/QSQ)/(V*(1.0+R))*VW2
      RETURN
      END
C...
