       DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL,modelhack)
C
C      INTEGRATION ROUTINE FROM SLIBSOUR
C
       implicit none

       external funct
       REAL*8 FUNCT,LOWER,UPPER,EPSLON
       INTEGER NLVL
       INTEGER LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I
       REAL*8 VALINT(50,2),MX(50),RX(50),FMX(50),FRX(50),FMRX(50),
     1 ESTRX(50),EPSX(50)
       REAL*8 R,FL,FML,FM,FMR,FR,EST,ESTL,ESTR,ESTINT,L,AREA,ABAREA,M,
     1 COEF,ROMBRG,EPS
	logical modelhack
C
C
       LEVEL = 0
       NLVL = 0
       ABAREA = 0.0
       L = LOWER
       R = UPPER
       FL = FUNCT(L,modelhack)
       FM = FUNCT(0.5*(L+R),modelhack)
       FR = FUNCT(R,modelhack)
       EST = 0.0
       EPS = EPSLON
100    LEVEL = LEVEL + 1
       M = 0.5*(L+R)
       COEF = R-L
       IF(COEF.NE.0) GO TO 150
       ROMBRG = EST
       GO TO 300
150    FML = FUNCT(0.5*(L+M),modelhack)
       FMR = FUNCT(0.5*(M+R),modelhack)
       ESTL = (FL+4.0*FML+FM)*COEF
       ESTR = (FM+4.0*FMR+FR)*COEF
       ESTINT = ESTL+ESTR
       AREA= DABS(ESTL)+DABS(ESTR)
       ABAREA = AREA+ABAREA-DABS(EST)
       IF(LEVEL.NE.MAXLVL) GO TO 200
       NLVL=NLVL+1
       ROMBRG = ESTINT
       GO TO 300
200    IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1 (LEVEL.LT.MINLVL)) GO TO 400
       ROMBRG = (1.6D1*ESTINT-EST)/15.0D0
300    LEVEL = LEVEL - 1
       I = RETRN(LEVEL)
       VALINT(LEVEL,I) = ROMBRG
       GO TO (500,600),I
400    RETRN(LEVEL) = 1
       MX(LEVEL) = M
       RX(LEVEL) = R
       FMX(LEVEL) = FM
       FMRX(LEVEL) = FMR
       FRX(LEVEL) = FR
       ESTRX(LEVEL) = ESTR
       EPSX(LEVEL) = EPS
       EPS = EPS/1.4
       R = M
       FR = FM
       FM = FML
       EST = ESTL
       GO TO 100
500    RETRN(LEVEL) = 2
       L = MX(LEVEL)
       R = RX(LEVEL)
       FL = FMX(LEVEL)
       FM = FMRX(LEVEL)
       FR = FRX(LEVEL)
       EST = ESTRX(LEVEL)
       EPS = EPSX(LEVEL)
       GO TO 100
600    ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)
       IF(LEVEL.GT.1) GO TO 300
       QUADMO = ROMBRG/12.0D0
       RETURN
       END
