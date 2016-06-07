C=======================================================================
 
      SUBROUTINE QUASIY8(E01,EP1,TH1,ia1, iz1, amuM1, avgM1,SIGMA_dble)                                
                                                                        
C Calculates quasielastic cross section 
c using Donnelly/Sick super scaling model
c but Paris w.f. for deuteron
      IMPLICIT NONE     
      REAL E0,EP,TH,SIGMA
      REAL  CHBAR,ALPHA,PM,PI
      real*8 e01, ep1, th1,ia1,iz1,amuM1,avgM1,sigma_dble
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828)                                       
      PARAMETER (PI    = 3.1415927)                                     
      REAL FF/1./  !phony initialization to please compiler: SER 4/15/93
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM 
      COMMON     /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                      
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL*4 PAULI_SUP1,PAULI_SUP2
      REAL*8 GEP,GEN,GMP,GMN
      real*8 qv,a,b,c,amf2,ams,disc,y,dydep,en,emp,denominator
      REAL THR,SINSQ,COSTH,QSQ,TAU,W1,W2,CMOTT,CSMOTT,RECOIL,DSIGE,FY,SD
      real*8 alfa,c1,c2
      real*8 kappa,lam,lamp,taup,squigglef,psi,psip,nuL,nut
      real*8 kf,es,GM2bar,GE2bar,W1bar,W2bar,Delta,GL,GT
      common/testing/prttst
      logical prttst
! 100*prob(Pz) from Paris wave function in 10 MeV bins from -1 to 1 GeV
! integeral is 1. Integral is 100.
      integer izz
      real*8 depdpz,pz
! using just w*w for d-state (tail too big)
      real*8 fydr(200)/
     > 0.00000,0.00002,0.00003,0.00005,0.00006,0.00009,0.00010,0.00013,
     > 0.00015,0.00019,0.00021,0.00026,0.00029,0.00034,0.00038,0.00044,
     > 0.00049,0.00057,0.00062,0.00071,0.00078,0.00089,0.00097,0.00109,
     > 0.00119,0.00134,0.00146,0.00161,0.00176,0.00195,0.00211,0.00232,
     > 0.00252,0.00276,0.00299,0.00326,0.00352,0.00383,0.00412,0.00447,
     > 0.00482,0.00521,0.00560,0.00603,0.00648,0.00698,0.00747,0.00803,
     > 0.00859,0.00921,0.00985,0.01056,0.01126,0.01205,0.01286,0.01376,
     > 0.01467,0.01569,0.01671,0.01793,0.01912,0.02049,0.02196,0.02356,
     > 0.02525,0.02723,0.02939,0.03179,0.03453,0.03764,0.04116,0.04533,
     > 0.05004,0.05565,0.06232,0.07015,0.07965,0.09093,0.10486,0.12185,
     > 0.14268,0.16860,0.20074,0.24129,0.29201,0.35713,0.44012,0.54757,
     > 0.68665,0.86965,1.11199,1.43242,1.86532,2.44703,3.22681,4.24972,
     > 5.54382,7.04016,8.48123,9.40627,9.40627,8.48123,7.04016,5.54382,
     > 4.24972,3.22681,2.44703,1.86532,1.43242,1.11199,0.86965,0.68665,
     > 0.54757,0.44012,0.35713,0.29201,0.24129,0.20074,0.16860,0.14268,
     > 0.12185,0.10486,0.09093,0.07965,0.07015,0.06232,0.05565,0.05004,
     > 0.04533,0.04116,0.03764,0.03453,0.03179,0.02939,0.02723,0.02525,
     > 0.02356,0.02196,0.02049,0.01912,0.01793,0.01671,0.01569,0.01467,
     > 0.01376,0.01286,0.01205,0.01126,0.01056,0.00985,0.00921,0.00859,
     > 0.00803,0.00747,0.00698,0.00648,0.00603,0.00560,0.00521,0.00482,
     > 0.00447,0.00412,0.00383,0.00352,0.00326,0.00299,0.00276,0.00252,
     > 0.00232,0.00211,0.00195,0.00176,0.00161,0.00146,0.00134,0.00119,
     > 0.00109,0.00097,0.00089,0.00078,0.00071,0.00062,0.00057,0.00049,
     > 0.00044,0.00038,0.00034,0.00029,0.00026,0.00021,0.00019,0.00015,
     > 0.00013,0.00010,0.00009,0.00006,0.00005,0.00003,0.00002,0.00000/
!*** NO D-sate (tail too small)
       real*8 fydn(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00030,
     > 0.00033,0.00038,0.00041,0.00047,0.00050,0.00057,0.00061,0.00068,
     > 0.00073,0.00081,0.00087,0.00095,0.00102,0.00111,0.00119,0.00128,
     > 0.00137,0.00147,0.00156,0.00167,0.00176,0.00189,0.00198,0.00210,
     > 0.00220,0.00232,0.00242,0.00254,0.00264,0.00276,0.00284,0.00296,
     > 0.00303,0.00313,0.00319,0.00327,0.00331,0.00338,0.00340,0.00344,
     > 0.00344,0.00347,0.00345,0.00348,0.00346,0.00349,0.00350,0.00358,
     > 0.00366,0.00386,0.00412,0.00454,0.00510,0.00594,0.00704,0.00859,
     > 0.01060,0.01334,0.01688,0.02156,0.02763,0.03544,0.04561,0.05878,
     > 0.07561,0.09741,0.12537,0.16167,0.20811,0.26892,0.34765,0.45096,
     > 0.58609,0.76539,1.00428,1.32177,1.75207,2.33167,3.11006,4.13194,
     > 5.42556,6.92172,8.36256,9.28786,9.28786,8.36256,6.92172,5.42556,
     > 4.13194,3.11006,2.33167,1.75207,1.32177,1.00428,0.76539,0.58609,
     > 0.45096,0.34765,0.26892,0.20811,0.16167,0.12537,0.09741,0.07561,
     > 0.05878,0.04561,0.03544,0.02763,0.02156,0.01688,0.01334,0.01060,
     > 0.00859,0.00704,0.00594,0.00510,0.00454,0.00412,0.00386,0.00366,
     > 0.00358,0.00350,0.00349,0.00346,0.00348,0.00345,0.00347,0.00344,
     > 0.00344,0.00340,0.00338,0.00331,0.00327,0.00319,0.00313,0.00303,
     > 0.00296,0.00284,0.00276,0.00264,0.00254,0.00242,0.00232,0.00220,
     > 0.00210,0.00198,0.00189,0.00176,0.00167,0.00156,0.00147,0.00137,
     > 0.00128,0.00119,0.00111,0.00102,0.00095,0.00087,0.00081,0.00073,
     > 0.00068,0.00061,0.00057,0.00050,0.00047,0.00041,0.00038,0.00033,
     > 0.00030,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/
! using 1.5 * (1-cos**2)**2 for dtate (better than either of above)
       real*8 fyd(200)/
     > 0.00000,0.00002,0.00002,0.00004,0.00004,0.00006,0.00007,0.00010,
     > 0.00011,0.00014,0.00015,0.00018,0.00020,0.00024,0.00026,0.00031,
     > 0.00033,0.00038,0.00042,0.00048,0.00052,0.00058,0.00063,0.00070,
     > 0.00076,0.00084,0.00090,0.00099,0.00107,0.00117,0.00125,0.00136,
     > 0.00145,0.00158,0.00168,0.00181,0.00192,0.00207,0.00219,0.00235,
     > 0.00249,0.00265,0.00280,0.00297,0.00313,0.00332,0.00347,0.00368,
     > 0.00385,0.00406,0.00424,0.00447,0.00467,0.00491,0.00513,0.00541,
     > 0.00565,0.00597,0.00626,0.00665,0.00703,0.00751,0.00803,0.00868,
     > 0.00938,0.01030,0.01135,0.01265,0.01421,0.01616,0.01848,0.02143,
     > 0.02496,0.02942,0.03489,0.04170,0.05015,0.06060,0.07374,0.09018,
     > 0.11064,0.13651,0.16893,0.21017,0.26205,0.32886,0.41420,0.52472,
     > 0.66768,0.85544,1.10347,1.43058,1.87107,2.46124,3.25006,4.28232,
     > 5.58540,7.08967,8.53679,9.46493,9.46493,8.53679,7.08967,5.58540,
     > 4.28232,3.25006,2.46124,1.87107,1.43058,1.10347,0.85544,0.66768,
     > 0.52472,0.41420,0.32886,0.26205,0.21017,0.16893,0.13651,0.11064,
     > 0.09018,0.07374,0.06060,0.05015,0.04170,0.03489,0.02942,0.02496,
     > 0.02143,0.01848,0.01616,0.01421,0.01265,0.01135,0.01030,0.00938,
     > 0.00868,0.00803,0.00751,0.00703,0.00665,0.00626,0.00597,0.00565,
     > 0.00541,0.00513,0.00491,0.00467,0.00447,0.00424,0.00406,0.00385,
     > 0.00368,0.00347,0.00332,0.00313,0.00297,0.00280,0.00265,0.00249,
     > 0.00235,0.00219,0.00207,0.00192,0.00181,0.00168,0.00158,0.00145,
     > 0.00136,0.00125,0.00117,0.00107,0.00099,0.00090,0.00084,0.00076,
     > 0.00070,0.00063,0.00058,0.00052,0.00048,0.00042,0.00038,0.00033,
     > 0.00031,0.00026,0.00024,0.00020,0.00018,0.00015,0.00014,0.00011,
     > 0.00010,0.00007,0.00006,0.00004,0.00004,0.00002,0.00002,0.00000/

       th=real(th1)
       e0=real(e01)
       ep=real(ep1)
       avgM=real(avgM1)
       iz=int(iz1)
       ia=int(ia1)
       amuM=real(amuM1)
       
       SIGMA = 0.
       avgN=iA-iZ
       avgA=real(iA)

       IG=12
       IDUT=0
       INEL_MODEL=1
       PAULI_MODEL=1
       NUC_METHOD=0
       NUC_MODEL=0

!       write(*,*) 'start peter', E01,EP1,TH1,ia1, iz1, amuM1, avgM1
!       write(*,*) 'convert peter', E0,EP,TH,ia, iz, amuM, avgM
       IF (iA.EQ.1) RETURN
        
       
       SIGMA = 0.
       DYDEP = 0.
       AMS   = avgM-PM
       THR   = TH*PI/180.
       SINSQ = SIN(THR/2.)**2
       COSTH = COS(THR)
       QSQ   = 4.*E0*EP*SIN(THR/2.)**2
       TAU   = QSQ/4.0/PM**2
        
       CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)
        
c     if(prttst)
!       write(*,*) 'various param',   ig,e0,ep,th,qsq,gep,gmp,gmn,gen
       
       
!     Pauli suppression model
       CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,real(QSQ),real(E0),
     >   PAULI_SUP1,PAULI_SUP2)
!       write(*,*) 'supression factors', pauli_sup1, pauli_sup2
       W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)
     +   
       W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU)  
c     if(prttst) write(8,'(10f7.3)') pauli_sup2,gep**2,gen**2,w1,tau,w2
       
       CMOTT  = CHBAR**2*0.001*ALPHA**2/4.
       CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2
       RECOIL = PM/(PM+E0*(1.-COSTH))
       DSIGE  = (W2+2.0*TAN(THR/2.)**2*W1)*CSMOTT*RECOIL
!       write(*,*) 'modd and stuff', cmott, csmott, recoil, dsige,w1,w2
       
!     Modifed to use Superscaling from Sick, Donnelly, Maieron,
!     nucl-th/0109032
       if(IA.eq.2) kf=0.085
       if(iA.eq.2) Es=0.0022
       if(IA.eq.3) kf=0.150
       if(iA.eq.3) Es=0.010 
       if(IA.eq.4) kf=0.200
       if(iA.eq.4) Es=0.015 
       if(IA.gt.4) kf=0.165
       if(iA.gt.4) Es=0.015 
       if(IA.gt.7) kf=0.228
       if(iA.gt.7) Es=0.020 
       if(IA.gt.16) kf=0.230
       if(iA.gt.16) Es=0.025 
       if(IA.gt.25) kf=0.236
       if(iA.gt.25) Es=0.018 
       if(IA.gt.38) kf=0.241
       if(iA.gt.38) Es=0.028 
       if(IA.gt.55) kf=0.241
       if(iA.gt.55) Es=0.023 
       if(IA.gt.60) kf=0.245
       if(iA.gt.60) Es=0.028 
       if(ep.ge.e0) then 
         write(6,'(1x,''error,e0,ep='',2f8.3)') e0,ep
         return
       endif
       QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)
     +   
       kappa = qv / 2. / pm
       lam = (e0 - ep) / 2. / pm
       if(abs(kappa**2 - lam**2 - tau).gt.0.01) then
         write(6,'(1x,''error,tau='',3f8.3)') kappa**2, lam**2,tau
         return
       endif
       lamp = lam - Es / 2. / pm
       taup = kappa**2 - lamp**2
       squigglef = sqrt(1. + (kf/pm)**2) -1.
       if(1.+lamp.le.0.) then
         write(6,'(1x,''error,lamp='',3f8.3)') lam,lamp
         return
       endif
       if(taup * (1. + taup).le.0.) then
         write(6,'(1x,''error,taup='',3f8.3)') kappa**2, lam**2,tau
         return
       endif
       psi =  (lam  - tau ) / sqrt(squigglef) /
     >   sqrt((1.+lam )* tau + kappa * sqrt(tau * (1. + tau)))
       psip = (lamp - taup) / sqrt(squigglef) / 
     >   sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))
       nuL = (tau / kappa**2)**2
       nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
       GM2bar = Pauli_sup1 * (iZ * GMP**2 + avgN * GMN**2)  
       GE2bar = Pauli_sup2 * (iZ * GEP**2 + avgN * GEN**2) 
       W1bar = tau * GM2bar
       W2bar = (GE2bar + tau * GM2bar) / (1. + tau)
       Delta = squigglef * (1. - psi**2) * (
     >   sqrt(tau * (1.+tau)) / kappa + squigglef/3. *
     >   (1. - psi**2) * tau / kappa**2)
       GL = kappa**2 / tau * (GE2bar + Delta * W2bar) / 
     >   2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
       GT = (2. * tau * GM2bar + Delta * W2bar) /
     >   2. / kappa / (1. + squigglef * (1. + psi**2) / 2.)
!     from Maria Barbaro: see superfy.m1
       FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))
       sigma = csmott * FY * (nuL * GL + nuT * GT) / kf 
       if(sigma.lt.0.) write(6,'(''ERROR, sigma...='',10f8.1)')
     >   sigma,csmott,fy,nuL,GL,nuT,GT
       if(prttst) write(8,'(/1x,12f7.3)') e0,ep,th,W1,W2,csmott/10000.,
     >   fy,nul,gl,nut,gt,sigma/100.
       
!     Use PWIA and Paris W.F. for deuteron
       if(IA.eq.2) then
         pz = (2.*pm*(e0-ep) - qsq) / 2. / qv
         izz = int((pz + 1.0) / 0.01) + 1
         izz = min(200,max(1,izz))
         depdpz = ep * qv / pm / e0
!     Facotr of 100 and 0.01 cance to give:
         sigma = dsige * fyd(izz) / depdpz  
         if(prttst) write(8,'(i4,4f7.3)') izz,pz,depdpz,
     >     fyd(izz),sigma/100.
       endif
!       write(*,*) 'got sigma of ', sigma
       sigma_dble=dble(sigma)
       return
       
! Are we changing quasi-elastic cross section due to Asymmetry?         
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0)                                
c    >     .AND.(INDEX(TARGET,'_P').GT.0) )                             
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >  CSMOTT*RECOIL,DSIGE)                                            
                                                                        
       QV    = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH)
     +   
       AMF2  = PM**2
     +   
       A     = 4.*QV**2-4.*(E0-EP+avgM)**2
     +   
       B     = -4.*QV*(AMS**2-AMF2-QV**2+(E0-EP+avgM)**2)
     +   
       C     = (E0-EP+avgM)**4+(AMS**2-AMF2-QV**2)**2-2.*(E0-EP+avgM)**2
     +   *(AMF2+QV**2+AMS**2)                    
       DISC  = B**2-4.*A*C
     +   
       IF (A.EQ.0..OR.DISC.LT.0.) RETURN
     +   
       Y     = (-B-SQRT(DISC))/2./A
     +   
       
C     JACOBIAN DY/DEP
C     
       
       EN    = SQRT(AMF2+(QV+Y)**2)
     +   
       EMP   = SQRT(AMS**2+Y**2)
     +   
       DENOMINATOR = (Y/EMP+(QV+Y)/EN)
     +   
       IF (DENOMINATOR.LE.1.E-30) RETURN
     +   
       DYDEP = (1.+(QV+Y)/EN*(EP-E0*COSTH)/QV)/DENOMINATOR
     +   
       
C     GAUSSIAN SMEARING FUNCTION
C     
       
       FY = 0.
     +   
       IF (iA.LE.4) THEN
     +   
!     *** 034 way too small for 4He: swith to 100 
c**   SD = 0.034                                         
         SD = 0.100
       ELSE
     +     
         SD = 0.119                                                   
       END IF
     +     
c     special case for C
c***  if(IA.eq.12) SD=0.09
       
       fy=0.
       IF ((Y**2/2./SD**2).LT.40.)FY = 1./SD/SQRT(2.*PI)*EXP(-Y**2/2./SD
     +   **2)                   
       
!     Modified to use fit to FY from nucl-th/9812078 (C. delgi Atti...)
!     Modified 11/05 to get better agreement with Naida Fromin's data
       if(IA.eq.2) then
         a = 6.1 - 0.50
         alfa = .045 
         c1 = 0.018 
         c2 = 0.25 
       endif
       if(IA.eq.3) then
         a = 7.1 - 1.00
         alfa = .083
         c1 = 0.041 * 1.00
         c2 = 0.33 * 1.20
       endif
       if(IA.eq.4) then
         a = 6.8 - 1.00
         alfa = .167
         c1 = 0.106 * 1.20
         c2 = 0.65 * 1.20
       endif
       if(IA.gt.4.and.IA.le.30) then
         a = 5.1
         alfa = .166
         c1 = 0.083
         c2 = 0.57 * 1.20
       endif
       if(IA.gt.30) then
         a = 4.6
         alfa = .138
         c1 = 0.058
         c2 = 0.62 * 1.30
       endif
!***  peb pyb changed exponent from -6 to -8 10/05
!***  to try to get better agreement with data
!     changed to 8 on 11/15/05
       fy = c1 * exp(-1.0 * a**2 * y**2) / (alfa**2 + y**2) +
     >   c2 * exp(-8.0 * abs(y))
       
       SIGMA = DSIGE*FY*DYDEP
        
       sigma_dble=dble(sigma)
      
       RETURN
        
       END        


C====================================================================== 
                          

      SUBROUTINE PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0, PAULI_SUP1,PAULI_SUP2)
!-----------------------------------------------------------------------
! Gets Pauli Suppression factor for quasi-elastic scattering from
! Several different Models.
! The VanOrden Model is very Slow.
! Used by both INTERNAL and EXTERNAL
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*4 QSQ,E0,PAULI_SUP1,PAULI_SUP2,FDEL
      REAL Q,TAU,FSHELL,FGAUSS,FF_BESSEL,XX
      INTEGER PAULI_MODEL,NUC_MODEL,iA
      REAL FF
      REAL*4 MP24/3.52/   
      LOGICAL OUT_OF_RANGE
      REAL*4 P_FERMI/.25/

!      write(*,*) 'suppressing pauli' , pauli_model, nuc_model, iA, qsq, e0

      IF(PAULI_MODEL.EQ.0) THEN
        FF = 0.                                                           
        IF (iA.EQ.2.AND.QSQ.LE.8.)  FF = FDEL(QSQ) 
        IF(IA.GT.2) THEN
         OUT_OF_RANGE =.FALSE.
         IF(NUC_MODEL.EQ.1) THEN
          FF = FF_BESSEL(REAL(QSQ),OUT_OF_RANGE) !only for some Nuclei,feor limited Q2
         ENDIF
         IF(OUT_OF_RANGE.OR.(NUC_MODEL.EQ.0)) THEN !use if FF_BESSEL out of range
          IF (iA.LE.20) THEN
            FF  = FSHELL(QSQ)
          ELSE     !ia >20
            FF  = FGAUSS(QSQ) 
          ENDIF
         ENDIF  
        ENDIF                           
        Pauli_sup2 = (1.-FF**2)    !old model from Stein
        PAULI_sup1 =1.
      ELSE
        IF(PAULI_MODEL.EQ.1) THEN !Tsai RMP 46,816(74) eq.B54
!         write(*,*) 'woo-hoo'
          TAU   = QSQ/MP24     
          Q=SQRT(QSQ*TAU+QSQ)
          IF((Q .GT. 2.*P_FERMI).OR.(iA.EQ.1)) THEN
            PAULI_SUP2 =1.0
          ELSE
            PAULI_SUP2=3.0*Q*(1.0-0.08333*(Q/P_FERMI)**2)
          ENDIF
!          write(*,*) 'got factor 2 ', pauli_sup2
        ELSEIF(PAULI_MODEL.EQ.3) THEN
          CALL  Q_E_VANORDEN(QSQ,E0,PAULI_SUP2)
        ELSEIF(PAULI_MODEL.EQ.2) THEN ! no suppression
          PAULI_SUP2 =1.0       ! No suppression
        ELSE
          PAULI_SUP2 =1.0       ! No suppression
        ENDIF
        PAULI_SUP1= Pauli_sup2
      ENDIF
      IF(PAULI_MODEL.EQ.4) THEN ! Louk's formula for Deuterium
        TAU   = QSQ/MP24     
        Q=SQRT(QSQ*TAU+QSQ)
        xx = q * 1000. / 100.
        Pauli_sup2 = (1.  - 1. / (1. + xx * (0.0060 + xx * (0.3882 + 
     >    xx * 0.2477))))
        PAULI_SUP1= Pauli_sup2
      endif
           
!     write(*,*) 'about to return ', pauli_sup1, pauli_sup2
      RETURN
      END

C====================================================================== 
   
C=======================================================================
                                                                        
      SUBROUTINE SECNUCLW(E0,EP,TH,SIGMA)                               
!--------------------------------------------------------------------
! SIGMA is cross section nb/GeV/str/nucleus.                                       
C Bastardization of Javier's SECNUC, inspired by patches, patches, and  
C patches splattered all over XSECFIT, INEFT, and this routine. Target  
C information is passed through common /targt/.                         
                                                                        
C SECNUCLW calculates the nuclear cross section in the deep inelastic   
C region using a variety of fits.                                  
C Added Fermi smearing Sept. 2005 P. Bosted
C Changed arguement of INELAST to WSQ rather than W
c Changed 8/06 to use better smearing pyb
!--------------------------------------------------------------------
      implicit none
      real e0,ep,th,sigma,sinsq,cossq,tansq,qsq,wsq,w,csmott,f1c,x
      real avgn,avga,avgM,amuM,r,dr,nu,eps,kappa,sigres,flux
      REAL*8 W1,W2,sigt,rc
      INTEGER ISM
      REAL*8 W1p,W2P,DW2DPF,WSQP,pf,kf,es,dw2des
      logical goodfit

! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/

      integer iz,ia,i
      COMMON      /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                     
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/
      common/testing/prttst
      logical prttst,first/.true./
      real*8 xval1(50),xvall(50),temp(4)
! This if Fermi smearing for deuteron based n Pairs w.f.
! total probability is 10. Bins are 40 MeV wide in pz
      integer izz
      real*8 pz,qv
! use w**2 for d-state
      real*8 fydr(15)/0.0192,0.0303,0.0538,0.1091,0.2544,0.6857,2.0061,
     >  3.5775,2.0061,0.6857,0.2544,0.1091,0.0538,0.0303,0.0192/
! *** with NO D-sate!
       real*8 fydn(15)/0.0040,0.0102,0.0277,0.0764,0.2150,0.6410,1.9589,
     >  3.5301,1.9589,0.6410,0.2150,0.0764,0.0277,0.0102,0.0040/
! using 1.5 * (1-cos**2)**2 for dtate
       real*8 fyd(15)/ 0.0094,0.0187,0.0411,0.0970,0.2462,0.6866,2.0207,
     >  3.6003,2.0207,0.6866,0.2462,0.0970,0.0411,0.0187,0.0094/
                                                                        
C Calculate QSQ and WSQ for this kinematic point:                       
                                                                        
      SINSQ  = SIN(TH*3.1415927/180./2.)**2                             
      COSSQ  = 1.0-SINSQ                                                
      TANSQ  = SINSQ/COSSQ                                              
      QSQ    = 4.0*E0*EP*SINSQ                                          
      WSQ    = PM*PM+2*PM*(E0-EP)-QSQ                                   
      NU     = E0 - EP
      sigma = 0.

      CSMOTT = (.001)*(19732./(2.0*ALPHA*E0*SINSQ))**2*COSSQ            
                                                                        

! Apply Fermi smearing to inelatic, except for H2
      IF(IA .eq. 1) THEN
! arguements of INELAST are real*8
        CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1,W2)!get s.f. /nucleon 
      endif
      if(IA.ge.2) then
        W1 = 0.
        W2 = 0.
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
        if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
        if(IA.eq.3) kf=0.150
        if(iA.eq.3) Es=0.010 
        if(IA.eq.4) kf=0.200
        if(iA.eq.4) Es=0.015 
        if(IA.gt.4) kf=0.165
        if(iA.gt.4) Es=0.015 
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
! adjust pf to give right width based on kf
        pf = 0.5 * kf 
! assume this is 2 * pf * qv

        qv = SQRT(E0*E0 + EP*EP - 2. * E0 * EP *
     >     COS(TH * 3.1415927 / 180.))
        DW2DPF = 2. * qv
        dw2des = 2. * (e0 + PM - ep) 
        if(prttst) write(8,'(1x,''pf,dw2dpf'',4f10.3)') pf,
     >    dw2dpf,qv,qv - sqrt(nu**2 + qsq)
        do ism = 1,15
           WSQP = WSQ + XX(ISM) * PF * DW2DPF - es * dw2des
           if(IA.eq.2) wsqp = wsq + (-0.3 + 0.6 * 
     >        (float(ism)-0.5)/15.) * dw2dpf
           IF(WSQP.GT. 0.) THEN
            CALL INELAST(DBLE(QSQ),WSQP,W1P,W2P)!get s.f. /nucleon 
            if(IA.eq.2) then
              W1 = W1 + W1P * FYD(ISM)/10.
              W2 = W2 + W2P * FYD(ISM)/10.
              if(prttst) write(8,'(1x,''insm d'',10f8.3)') 
     >          qsq, wsq,wsqp,w1,w1p,fyd(ism)/10.
            else
              W1 = W1 + W1P * FY(ISM)
              W2 = W2 + W2P * FY(ISM)
              if(prttst) write(8,'(1x,''insm'',10f8.3)') 
     >          qsq, wsqp,w1,w1p,fy(ism)
            endif
          ENDIF
        ENDDO
        CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1P,W2P)!get s.f. /nucleon 
        if(prttst) write(98,'(1x,''inelas sm'',6f10.3)') 
     >     qsq, w,w1,w1p,w2,w2p
      ENDIF
c      if(prttst) write(*,'(1x,''inelast got w1,w2='',2f7.3)') w1,w2                                                                  
      if(w1.lt.0.0.or.w2.lt.0.0) then
        write(6,'(1x,''error, e,ep,th,q2,w,w1,w2='',6f8.3)') e0,ep,
     >    th,qsq,w,w1,w2
      endif

      SIGMA = (W2+2.*TANSQ*W1)*CSMOTT
C             (per nucleon, including emc effect (which includes neutron
C              excess correction), if any)                              
                                                                        
      SIGMA  = avgA*SIGMA  !per nucleus                                               
                                                                        
! change cross section for polarized beam and target?                   
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0) )
c    >     .AND.( (INDEX(TARGET,'_P').GT.0).OR.                         
c    >            (INDEX(TARGET,'_A').GT.0) )   )  !inelastic only      
c    > CALL ASYM_INE(E0,EP,TANSQ,QSQ,TARGET,SIGMA)                      

      if(sigma.lt.0.0) write(6,'(1x,''error, sigma='',e12.4)') sigma


      RETURN                                                            
      END                                                               
         
C====================================================================== 
                                                                        
      SUBROUTINE INELAST(QQ,WSQ,W1,W2)                                    
! Choose which inelastic model is to be used for resonance and for DIS  
! Returns structure functions per NUCLEON (11/3/95 SER)
C CHANGED TO USE WSQ RATHER THAN W peb 10/05

      Implicit NONE
      REAL  avgN, avgA, avgM, amuM 
      COMMON /TARGT/ iZ, iA, avgN, avgA, avgM, amuM 
      INTEGER IZ,IA
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL   
      REAL*8 QQ,WSQ,W,X,W1,W2,F2_E665,R8,EMCFAC,F1C,RC
      REAL*8 MP2/0.8803/,MP/.93828/,F2,F2nF2p
      LOGICAL FIRST/.TRUE./                                             
      REAL Q_RES_LIMIT/.2/                                              
      REAL QQ4,X4,F2_4,W1_4,W2_4,FITEMC_N,R4,DR4,NU,F2Allm,FNP_NMC
      REAL rnphi,rnp,rnplo,fact
      LOGICAL GD
                                                                        
      W1 = 0.
      W2 = 0.
      if(WSQ.lt.(0.93828 + 0.13957)**2) return
      W = SQRT(WSQ)

      IF(INEL_MODEL.EQ.0) THEN                                          
       CALL INEFT(QQ,W,W1,W2)   !Bodek fit                              

! 8/06 changed to use proton fit for all nuclei too!
! modified by simple n/p model
      ELSEIF(INEL_MODEL.EQ.3) THEN 
c     CALL CHRISTY705(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
       CALL CHRISTY31606(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
       NU = (WSQ +QQ -MP2)/(2.*MP)
       W1 = F1C / MP
       W2 = W1 /(1.0 + NU*NU/QQ) * (1.0 + RC)
! Simple model for n/p. Constnat in resonance region,
! and DIS at high W. Smoothly join from W=1.7 to 2.3
! Fixed to work for any Z,A 8/06 pyb
       if(IA.gt.1) then
        rnplo = 0.66
        X4 = QQ/(WSQ -MP2 +QQ) 
        qq4 = qq
        rnphi = FNP_NMC(X4,QQ4)
        fact = max(0., min(1., (wsq-1.7)/0.6))
        rnp = rnplo * (1.-fact) + fact * rnphi
        W1 = W1 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
        W2 = W2 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
       endif

! Simona fit to deuteron
      ELSEIF(INEL_MODEL.EQ.4.and.IA.eq.2) THEN 
       CALL ResCsSim(WSQ,QQ,F1C) ! Simona H2 resonance fit
       X = QQ/(WSQ -MP2 +QQ) 
       CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)       
       RC = R4 ! She used R1998 and just fit F1 (same as sigt)
       NU = (WSQ +QQ -MP2)/(2.*MP)
       W1 = F1C / MP / 2.0
       W2 = W1 /(1.0 + NU*NU/QQ) * (1.0 + RC)

      ELSEIF((INEL_MODEL.EQ.9).OR.(INEL_MODEL.EQ.12)) THEN              
       X = QQ/(WSQ -MP2 +QQ)                                           
          IF(FIRST) THEN                                                
           FIRST=.FALSE.                                                
           WRITE(10,'('' LOWER Q**2 CUT OFF FOR STUART RESONACE FIT='', 
     >      F5.2)')Q_RES_LIMIT                                          
          ENDIF                                                         
       IF(QQ.GT.Q_RES_LIMIT) THEN
! Only H2 model. Calls F2GLOB_PETER for DIS.  
! Remember that F2GLOB_PETER returns structure functions/nucleon for H2 and D2
! Convert to Real*4 for this call
        QQ4=QQ
        X4=X                                       
        CALL INELSTU(QQ4, X4, F2_4, W1_4, W2_4, INEL_MODEL) !stuart for res, F2_Glob DI
        W2=W2_4
        W1=W1_4 
       ELSE  !out of range of inelstu.                                  
        CALL INEFT(QQ,W,W1,W2)   !Bodek fit                             
       ENDIF

      ELSEIF(INEL_MODEL.EQ.1) THEN  !NMC model for DIS, INEFT for resonance
        X4 = QQ/(WSQ -MP2 +QQ) 
        QQ4=QQ
        CALL NMCDIS(QQ4,X4,W1_4,W2_4,amuM)
        W2=W2_4
        W1=W1_4
      ELSEIF(INEL_MODEL.EQ.2) THEN   ! E665 Fit (includes resonances
        X = QQ/(WSQ -MP2 +QQ) 
        IF(amuM.LT.1.5) THEN
         CALL  F2PGLO(QQ,X,F2_E665,R8) ! Prot E665
         EMCFAC=1.
        ELSE  ! Deuteron or heavier 
         CALL  F2DGLO(QQ,X,F2_E665,R8) ! Deut E665 
         EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD)) !with neutron excess 8/19/98
        ENDIF
         CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
         NU = (W**2 +QQ -MP2)/(2.*MP)
         W2 = F2_E665/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
      ELSEIF(INEL_MODEL.EQ.6) THEN   ! F2allm
        X4 = QQ/(WSQ -MP2 +QQ) 
        qq4 = qq
        F2_4 = F2ALLM(X4,QQ4)   ! proton fit
        F2 = F2_4
        x = x4
        IF(amuM.GT.1.5) THEN
! Changed to FNP_NMC 7/19/06: seems to work much better!
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))    
          F2 = (FNP_NMC(X4,QQ4) + 1.) / 2.  * F2
        endif
        if(amuM.gt.2.5) then
          EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD))
          F2 = F2 * EMCFAC
        endif
        CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
        NU = (W**2 +QQ -MP2)/(2.*MP)
        W2 = F2/NU 
        W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
      ELSEIF(INEL_MODEL.GE.100) THEN !joining a DIS model with resonance model
        X4 = QQ/(W**2 -MP2 +QQ) 
        QQ4=QQ
        CALL F2_JOIN(QQ4,X4,amuM,iZ,INEL_MODEL,W1_4,W2_4,iA)
        W2=W2_4
        W1=W1_4
      ELSE                                                              
       WRITE(6,'('' NOT GOOD H2 INELASTIC MODEL='',I3)')INEL_MODEL        
       STOP                                                             
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

!-------------------------------------------------------------------------
      SUBROUTINE F2_JOIN(QQ,X,amuM,iZ,INEL_MODEL,W1,W2,iA) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon by joining the resonance
!  model with the DIS model
! amuM=atomic weight
! INEL_MODEL =rrdd where rr is resonance model number and dd is DIS model 
!             It must be >=100
! DIS_MODEL= 1 ineft
!            2 f2nmc
!            3 f2nmc95
!            4 SMC
!            5 E665
!            6 f2allm (HERA fit, obtained from Antje Burrel
!            9 f2glob model 9
!           12 f2glob model 12
!           33 rock_res9
! RES_MODEL  1 ineft
!            2 H2Model
!            3 rock_res9
!            5 E665
!            
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
! EMC effect corrected for neutron excess 8/19/98
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,ERR_LO,ERR_HI,NU,WSQ,f2allm
!      real fnp_nmc,corr
      real w1model,w2model,w1in,w2in,fnp_nmc
      INTEGER IZ,iA
      REAL MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC_N,FRAC,amuM,W1_DIS,W2_DIS,DSLOPE,ST,SLOPE,SY
      real f2nf2p,emcfacp,ROCK_RES9
      REAL*8 W18,W28,QQ8,W8,amuM8,F28,RDUM
      INTEGER INEL_MODEL,IHD,RES_MODEL,DIS_MODEL
      LOGICAL GD
      CHARACTER*1 HD(2)/'H','D'/
!** changed from 0.3 to 0.0 to avoid discontinuities. H2model
!** is well enough behaved that this should be fine
      REAL Q_RES_LIMIT/0.0/  
      logical ReallyJoin
      common/testing/prttst
      logical prttst

      RES_MODEL =INEL_MODEL/100
      DIS_MODEL =INEL_MODEL -RES_MODEL*100

      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      QQ8=QQ
      W8= SQRT(max(0.,WSQ))
      amuM8 =amuM 
      IHD =MIN(IA,2)   !no more than deuteron

      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC_N(X,REAL(IA),REAL(IZ),GD)) !with neutron excess 8/19/98
!*** If using F2ALL, use R1990 because this used to get these F2 values
!*** changed 1/19/02
      IF(DIS_MODEL.EQ.6) THEN
        CALL R1990F(X,QQ,R,DR)
      ELSE
        CALL R1998(X,QQ,R,DR,GD)
      ENDIF

!* If A>2 and doing Anje fit, only use F2ALL
      ReallyJoin=.true.
      if(IA.gt.2.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 
!* If RES_MODEL and DIS_Model both 6, don't join
      if(RES_MODEL.eq.6.and.DIS_MODEL.EQ.6) ReallyJoin=.false. 

!** Don't join if using rock_res9: should be good everywhere
      if(DIS_MODEL.EQ.33) ReallyJoin=.false.

      IF(WSQ.GT.3.5.or.(.not.ReallyJoin)) THEN        !DIS
	IF(DIS_MODEL.EQ.1) THEN   !Bodek fit  REAL*8      
          CALL INEFT(QQ8,W8,W18,W28) !already has EMC correction
          W1=W18
          W2=W28
        ENDIF
        IF((DIS_MODEL.GE.2.AND.DIS_MODEL.LE.6).or.
     >     DIS_MODEL.EQ.33 ) THEN
         IF(DIS_MODEL.EQ.2) CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
         IF(DIS_MODEL.EQ.3) CALL F2NMC_NEW(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.4) CALL F2SMC98(IHD,X,QQ,F2,ERR_LO,ERR_HI)
         IF(DIS_MODEL.EQ.5)  THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          F2 = F28
         ENDIF
         IF(DIS_MODEL.EQ.33) F2=ROCK_RES9(ihd,x,qq)
         W2 = F2/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF
!** Changed 1/19/02 to use Anje's n/p fit and nuclear dep.
!** Changed 7/20/06 back to NMC for n/p (works MUCH better!)
        IF(DIS_MODEL.EQ.6) THEN  !(HERA fit, obtained from Antje Burrel
          F2 = F2ALLM(X,QQ)   ! proton fit
          IF(iHD.EQ.2)  F2 = (FNP_NMC(X,QQ)+1.)/2. *F2
! This is what Antje uses:
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))
c          IF(IHD.EQ.2)  F2 = (F2NF2P+1.)/2. *F2
!* For C, Al, Cu, Au, override EMCFAC
          emcfacp=emcfac
          if(ia.eq.12) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** use carbon for Al since Antje's seems to low at high x
          if(ia.eq.27) emcfacp=(0.926-0.400*x-0.0987*exp(-27.714*x)+  
     >                 0.257*x**0.247)
! *** turn off Antje's Al corr., 
!         if(ia.eq.27) emcfacp=(0.825-0.46*x-0.19*exp(-21.8*x)+             !!!ALUMINUM  !!!
!    >                    (0.34*x**(-4.91))*x**5.0)
          if(ia.eq.65) emcfacp=(1.026-0.56*x-0.34*exp(-45.7*x)+             !!!COPPER   !!!
     >                    (0.26*x**(-4.41))*x**5.0) 
          if(ia.eq.197) emcfacp=(0.970-1.433*x-0.334*exp(-54.53*x)+          !!!GOLD    !!!
     >                  1.074*x**0.711)
!**** changed to use Rock emc, not Antje!
!***          W2 = F2/NU *EMCFACP
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
          if(prttst) write(*,'(1x,''W2,F2,R,emcfacp='',4f10.4)')
     >      W2,F2,R,emcfacp  
        ENDIF

        IF((DIS_MODEL.EQ.9).OR.(DIS_MODEL.EQ.12)) THEN
          CALL F2GLOB_PETER(X,QQ,HD(IHD),DIS_MODEL,F2,ST,SY,SLOPE,DSLOPE,GD)
          W2 = F2/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R)
        ENDIF
        if(prttst) write(*,'(1x,''model,x,q2,F2,r='',i3,5f10.4)')
     >    DIS_MODEL,x,qq,f2,r,emcfac  
      ENDIF

      IF(WSQ.LT.4.3.and.ReallyJoin) THEN !Resonance
        W1_DIS =W1
        W2_DIS =W2
        IF(RES_MODEL.EQ.1) THEN  !old Bodek resonance fit: EMC corrected
          CALL INEFT(QQ8,W8,W18,W28)
          W1=W18
          W2=W28
        ELSEIF(RES_MODEL.EQ.5) THEN
          IF(IHD.EQ.1) CALL  F2PGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Prot E665
          IF(IHD.EQ.2) CALL  F2DGLO(DBLE(QQ),DBLE(X),F28,RDUM) ! Deut E665 
          W2 = F28/NU *EMCFAC
          W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
        ENDIF 
        IF(RES_MODEL.EQ.2) THEN !Stuart Fit to Resonance, valid forQ2>.75
          CALL INEFT(QQ8,W8,W18,W28) !already had EMC, neutron corrected
          W1=W18
          W2=W28
          w1in = w1
          w2in = w2
c old way          IF((QQ.LT.10.).AND.(QQ.GE.Q_RES_LIMIT).AND.iA.lt.3) THEN
          IF(QQ.GE.0.3.AND.iA.lt.3) THEN
            if(ia.eq.1) CALL H2MODEL(QQ,WSQ,W1model,W2model)
            if(ia.eq.2) CALL D2MODEL_IOANA(QQ,WSQ,W1model,W2model)
            if(prttst) WRITE(*,'('' H2/D2Model,q2,wsq,w1,w2='',
     >        5f8.3)') QQ,WSQ,W1model,W2model,EMCFAC 
            if(qq.gt.0.5.and.qq.lt.15.) then
              w1=w1model
              w2=w2model
            endif
            if(qq.ge.0.3.and.qq.le.0.5) then
              frac=(qq-0.3)/(0.5-0.3)
              w1 = w1model * frac + w1in *(1.0-frac)
              w2 = w2model * frac + w2in *(1.0-frac)
            endif
            if(qq.ge.10.0.and.qq.le.15.) then
              frac=(qq - 10.0)/(15.0 - 10.0)
              w1 = w1in * frac + w1model *(1.0-frac)
              w2 = w2in * frac + w2model *(1.0-frac)
            endif
          ENDIF
        ENDIF 

        IF(WSQ.GT.3.5) THEN   !interpolate
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           if(prttst) WRITE(*,'(1x,''joining w2'',3f10.3)') 
     >       w2_dis,w2,frac
           W1 = W1_DIS*FRAC + W1*(1.0 - FRAC)
           W2 = W2_DIS*FRAC + W2*(1.0 - FRAC)
        ENDIF
      ENDIF
!**** If using Antje model, then also use Blok correction factors
!* taken out 1/19/02
!      if(DIS_MODEL.EQ.6) THEN
!        if(ihd.eq.1) call f2corr1h(x,corr)
!        if(ihd.eq.2) call f2corr2h(x,corr)
!        w1=w1*corr
!        w2=w2*corr
!      endif
      RETURN
      END

C====================================================================== 

      SUBROUTINE INELSTU(X,q2,F2,W1,W2,IMOD)
! New call to Linda Stuarts fit to resonance Larry's DIS fit
! Steve Rock 10/14/93

*******************************************************************************
* This subroutine returns W1 and W2, the inelastic structure functions
*
* 2/93, LMS.
******************************************************************************
       IMPLICIT NONE

       INTEGER IMOD
       REAL X, Q2, NU, W1, W2,  MP/0.93828/, MP2/.8803/,
     >        WSQ, WW1, WW2, FRAC
       REAL   F2,R,DR,ST,SY,SLOPE,DSLOPE
       LOGICAL GOOD

       W1 = 0.0
       W2 = 0.0

       WSQ = MP2   + Q2*(1./X - 1. )
       NU = (WSQ +Q2 -MP2)/(2.*MP)
       IF(WSQ.LT.1.16) RETURN

! Get proton structure functions
       CALL R1998(X,Q2,R,DR,GOOD)
       IF(WSQ.GT.3.5) THEN
         CALL F2GLOB_PETER(X,Q2,'H',IMOD,F2,ST,SY,SLOPE,DSLOPE,GOOD)
         W2 = F2/NU
         W1 = (1.0 + NU*NU/Q2)*W2/(1.0 + R)
       ENDIF
! H2 model only officially works above q2=.75
       IF(Q2.LT.10.0.AND.WSQ.LT.4.3) THEN
         IF(WSQ.LE.3.5) THEN
           CALL H2MODEL(Q2,WSQ,W1,W2)
           F2 = NU*W2
         ELSE
           CALL H2MODEL(Q2,WSQ,WW1,WW2)
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + WW1*(1.0 - FRAC)
           W2 = W2*FRAC + WW2*(1.0 - FRAC)
         ENDIF
       ENDIF
       F2 = NU*W2

       RETURN
       END

C====================================================================== 

!-----------------------------------------------------------------------
      SUBROUTINE NMCDIS(QQ,X,W1,W2,amuM) 
!----------------------------------------------------------------------
! all arguements are Real*4  
! Calculates structure functions per nucleon. 
!  NMC fit for DIS and Old Bodek fit for Resonance.
! IA=atomic weight
! Makes EMC effect correction for heavy nuclei
!                 Steve Rock 11/95
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL QQ,X,W1,W2,R,DR,F2,ERR_F2,NU,WSQ, MP/0.93828/, MP2/.8803/
      REAL EMCFAC,FITEMC,FRAC,amuM
      REAL*8 W18,W28,QQ8,W8,amuM8
      INTEGER IA,IHD
      LOGICAL GOOD


      W1 = 0.0
      W2 = 0.0

      WSQ = MP2   + QQ*(1./X - 1. )
      IF(WSQ.LT.1.16) RETURN
      NU = (WSQ +QQ -MP2)/(2.*MP)
      EMCFAC = ABS(FITEMC(X,amuM,GOOD)) !isoscaler nucleus
      CALL R1998(X,QQ,R,DR,GOOD)
      IA = amuM+.5  !make into integer
      amuM8 =amuM
      IF(WSQ.GT.3.5) THEN
        
        IHD =MIN(IA,2)   !no more than deuteron
        CALL F2NMC(IHD,X,QQ,F2,ERR_F2)
        W2 = F2/NU *EMCFAC
        W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R) 
      ENDIF

      IF(WSQ.LT.4.3) THEN
         QQ8=QQ
         W8= SQRT(max(0.,WSQ))
         IF(WSQ.LE.3.5) THEN
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit: EMC corrected
           W1=W18
           W2=W28
         ELSE  !interpolate
           CALL INEFT(QQ8,W8,W18,W28) !old resonance fit       
           FRAC = (WSQ - 3.5)/(4.3 - 3.5)
           W1 = W1*FRAC + W18*(1.0 - FRAC)
           W2 = W2*FRAC + W28*(1.0 - FRAC)
         ENDIF
      ENDIF
      F2 = NU*W2
      RETURN
      END


      SUBROUTINE H2MODEL(QQ4,WW4,W1,W2)

********************************************************************************
*
* This subroutine calculates model cross sections for inclusive electron-proton
* scattering in the resonance region. The cross section is returned in 
* nanobarns/sr GeV. This fit is a modified version of Linda Stuart's 8/91 fit. 
* One major difference is that the coefficients represent fit results from a 
* substantially expanded data set, including all inclusive SLAC data in the 
* range 1.1 < W^2 < 5. There are other differences; for a complete discussion, 
* see Keppel's Ph.D. thesis. 2/94 CEK
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
* SIG1     = Cross section in nb/sr/GeV**2 (DSIG/DOMEGA/DW**2)
* SIG2     = Cross section in nb/sr/GeV (DSIG/DOMEGA/DEP)
* SIG_NRES = Non-resonant contribution to SIG1.
* SIG_RES  = Resonant contribution to SIG1.
* SIG_NR   = Non-resonant contribution to SIG2.
* SIG_R    = Resonant contribution to SIG2.
* SIGroper = Possible Roper contribution to SIG2.
* goroper  = Logical variable, set true if including possible Roper strength
*
* SIG1, SIG2, SIG_NRES, and SIG_RES are 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
        IMPLICIT NONE

        logical goroper 
        logical goodfit 
        INTEGER I
        REAL*4  QQ4,WW4,W1,W2
        REAL*8  SIN2, SIG1(2), SIG2(2), SIG_RES(2), SIG_R(2), 
     >    SIG_RES1(2), SIG_RES2(2), SIG_NRES(2), SIG_NR(2), 
     >    DIPOLE, COS2, TAU, EPS, K, DW2DEP, DEPCONV, 
     >    PI, AM, ALPHA,NU,CONV,SIGT,SIGL,qq,ww,
     >    RADCON, R_NRES, sigroper(2), x, r, dr,fac,rlog
        
        qq = qq4
        ww = ww4
        
        RADCON = 3.141592654/180.
        PI = 3.14159265
        AM = .9382727
        ALPHA = 0.00729735
        CONV = 0.0025767

C        SIN2 = SIN(TH*RADCON/2.0)*SIN(TH*RADCON/2.0)
C        COS2 = 1 - SIN2
C        Q2 = 4.0*E*EP*SIN2
C        W2 = AM*AM + 2.0*AM*(E - EP) - 4.0*E*EP*SIN2

         goroper = .true.

C        write(6,*) q2,w2

        IF(WW.LT.1.15) THEN            ! Below pion threshold
          DO I = 1,2
            SIG1(I) = 0.0
            SIG2(I) = 0.0
            SIG_NRES(I) = 0.0
            SIG_RES(I) = 0.0
           SIG_NR(I) = 0.0
            SIG_R(I) = 0.0
            sigroper(I) = 0.0
          ENDDO
          RETURN
        ENDIF

        K = (WW - AM*AM)/(2.0*AM)
        NU = K + QQ/(2.0*AM)
        TAU = NU*NU/QQ
        x = qq/2./AM/nu
C        EPS = 1.0/(1.0 + 2.0*(1.0 + TAU)*SIN2/COS2)
        DIPOLE = 1.0/(1.0 + QQ/0.71)**2
c        DEPCONV = ALPHA*K*EP/(4.0*PI*PI*Q2*E)*(2.0/(1.0 - EPS))*1000.
c        DW2DEP = 2.0*AM + 4.0*E*SIN2

! H2MOD_FIT returns cross sections in units of microbarns/(dipole FF)**2
        CALL H2MODEL_FIT(QQ,WW,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)

         R_NRES = 0.25/SQRT(QQ)


        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
C        SIG_NRES(2) = SIG_NRES(2)*DIPOLE*DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1))
     >                 *DIPOLE*DIPOLE


         SIGT = SIG_NRES(1)+SIG_RES(1)
         SIGL = R_NRES*SIG_NRES(1)
         W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
         W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END
C====================================================================== 


      SUBROUTINE H2MODEL_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J, KK
      REAL*8 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     WR, KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, GAM, 
     >     GR, RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER, ERRCHECK

      REAL*8 PI, ALPHA, AM, MPPI, MDELTA, MPI, GAM1WID, GAM2WID, MASS1,
     >     ROPERWID, MASSROPER, MASS2, DELWID, FIT(30), SQRTWDIF,
     >     XR, MQDEP

       REAL*8 RCOEF(27), COEF(22), RERR1(200),rerr2(125), ERR1(200), 
     >        ERR2(53)
    
      LOGICAL FIRST

      FIRST = .TRUE.

      
      pi = 3.14159265
      alpha = 0.00729735
      am = .9382727
      MPPI = 1.07783
      MDELTA = 1.229 
      MPI = 0.13957 
      GAM1WID = 0.080
      GAM2WID = 0.090
      MASS1 = 1.5062 
      ROPERWID = 0.0500 
      MASSROPER = 1.4000 
      MASS2 = 1.6810
      DELWID = 0.120
      XR = 0.1800 
      MQDEP = 3.40 

      DATA RCOEF/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >   4.6273E-01,   1.7844E-01,   -1.5416E+02, -1.4891E+02,
     >   2.4102E+02,   2.5823E+00,    7.1004E+00, -8.9771E+00,
     >   1.3744E+00,  -1.2085E+00,    1.1218E-01/

      DATA COEF/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/
 

! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(25) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
        FIT(I + 12) = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 18)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 21)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO
      if (goroper) FIT(26)  = FIT(25)/sqrt(Q2)
      if (goroper) FIT(27)  = FIT(25)/Q2

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,27
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12.OR.J.EQ.17
     > .OR.J.EQ.18 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
          elseIF(j.ge.25.and.j.le.27) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
          ENDIF
        ENDDO
      endif



      RETURN
      END

C
      Subroutine F2NMC(T,x,qsq,F2,err_f2) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in Phy Lett. b295 159-168 Proton and Deuteron F2 Structure Functions  
! in deep inelastic muon scattering   P. Amaudruz et. al.               
!                                                                       
      
                                                                        
      real a(7,2)/-0.1011,2.562,0.4121,-0.518,5.967,-10.197,4.685,      
     >            -0.0996,2.489,0.4684,-1.924,8.159,-10.893,4.535/      
      real b(4,2)/0.364,-2.764,0.0150,0.0186,                           
     >            0.252,-2.713,0.0254,0.0299/                           
      real c(4,2)/-1.179,8.24,-36.36,47.76,                             
     >            -1.221,7.50,-30.49,40.23/                             
      real Lam/0.25/                                                    
      real Qsqo/20.0/                                                   
      real log_term                                                     
      real A_x,B_x,C_x,F2,err_f2,x,qsq                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
                                                                        
      A_x=x**a(1,t)*(1-x)**a(2,t)*(a(3,t)+a(4,t)*(1-x)+a(5,t)*(1-x)**2  
     >  +a(6,t)*(1-x)**3+a(7,t)*(1-x)**4)                             
      B_x=b(1,t)+B(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = alog(qsq/lam**2)/alog(qsqo/lam**2)                     
      F2=A_x*(log_term)**B_x*(1+C_x/qsq)                                
      
      return                                                            
      end                                                               

C====================================================================== 

      SUBROUTINE F2GLOB_PETER(X,Q2,Target,MODEL,F2,ST,SY,slope,dslope,GOODFIT
     +  )
                                                                        
!     Returns F2 and related quantities from the either the LAMBDA12 model  
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.           
!                                                                       
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2  
!                                                                       
! Further, program returns uncertainty in F2 based on both statistics   
! (ST) and systematic effects due to our choice of models (SY).         
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical
! uncertainty in this slope.                                            
!                                                                       
! Best model is LAMBDA12.  SY is estimated to be the difference between 
! the two models.                                                       
!                                                                       
! Errors due to overall normalization are not included in program output
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.      
! Systematic errors due to radiative corrections are shown in Reference 
! to be very kinematic independent, and are everywhere <.5%, and thus,  
! are ignored by this routine (also see documentation to dFRC in file   
! HELP.DOCUMENT).                                                       
!                                                                       
! Coefficients and correlation matrix elements are from File            
! E.13 F1990.MATRICES, dated 27Jan90.                                   
                                                                        
      IMPLICIT NONE                                                     
      LOGICAL GOODFIT,FIRST/.TRUE./                                     
      REAL   X, Q2, F2,ST,SY,SLOPE,DSLOPE                               
      REAL   XPP,XP,Y,POLY,F2B,POL1,STB,Q,QTH,DQ,F2TH,QUAD,SCB,F2L      
      REAL   BINDING,STL                                                
      INTEGER MODEL, I, J, K                                            
      CHARACTER*1 TARGET                                                
      REAL*8    B(2,9),L(2,9,9) ! OMEGA9 variables 
      REAL*8    C(2,12),M(2,12,12) ! LAMBDA12 variable
      REAL*8 V(9),Z(12),U(12),LIN                                       
      
      
!     Model #9  27 Jan 90.                                                  
      DATA (B(1,J),J=1,9)/      !  HYDROGEN Ci         
     >  0.7338659870D0,  11.0245522588D0,   2.6185804129D0,            
     >  4.0956321483D0,   0.1206495422D0,   1.9714128709D0,            
     >  3.8893348719D0, -14.0507358314D0,   8.8080576075D0/            
      DATA ((L(1,J,K),K=1,9),J=1,9)/ !  HYDROGEN MijDiDj   
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,             
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,             
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,             
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,             
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,             
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,             
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,             
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,             
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,             
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,             
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,             
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,             
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,             
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,             
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,             
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,             
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,             
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,             
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,             
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,             
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,             
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,             
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,             
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,             
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,             
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,             
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/             
                                                                        
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci       
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,             
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,             
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/             
                                                                        
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj  
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,             
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,             
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,             
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,             
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,             
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,             
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,             
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,             
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,             
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,             
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,             
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,             
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,             
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,             
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,             
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,             
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,             
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,             
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,             
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,             
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,             
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,             
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,             
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,             
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,             
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,             
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/             
                                                                        
!    MODEL #12:    27Jan90.                                             
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci            
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,             
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,             
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,             
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/             
      DATA ((M(1,J,K),K=1,12),J=1,12)/     !     HYDROGEN MijDiDj       
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,             
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,             
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,             
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,             
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,             
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,             
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,             
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,             
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,             
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,             
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,             
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,             
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,             
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,             
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,             
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,             
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,             
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,             
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,             
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,             
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,             
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,             
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,             
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0,             
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,             
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,             
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,             
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,             
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,             
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,             
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,             
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,             
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,             
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,             
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,             
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,             
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,             
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,             
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,             
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,             
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,             
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,             
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,             
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,             
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,             
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,             
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,             
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/             
                                                                        
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci           
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,             
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,             
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,             
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/             
      DATA ((M(2,J,K),K=1,12),J=1,12)/     !     DEUTERIUM MijDiDj      
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,             
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,             
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,             
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,             
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,             
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,             
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,             
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,             
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,             
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,             
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,             
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,             
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,             
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,             
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,             
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,             
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,             
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,             
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,             
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,             
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,             
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,             
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,             
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0,             
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,             
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,             
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,             
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,             
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,             
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,             
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,             
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,             
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,             
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,             
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,             
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,             
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,             
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,             
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,             
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,             
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,             
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,             
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,             
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,             
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,             
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,             
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,             
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/               
      !---------------------------------------------------------------- 
      !---------------------------------------------------------------- 
                                                                        
                                                                        
                                                                        
      i = 1                                                             
      IF (TARGET.EQ.'D') i = 2                                          
      BINDING = 1./(1.-EXP(-MIN(20.,7.7*(1./X+.93828**2/Q2-1.))))       
      IF (i.EQ.1) BINDING = 1.                                          
                                                                        
      !OMEGA9 MODEL FIRST:                                              
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                             
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                             
           Y    = 1.-XP                                                 
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+      
     >            B(i,9)*Y**7                                           
           F2B  = X/XPP*BINDING*POLY                                    
          !-----------------------------------------------------------  
          !V(k) is the derivative of F_2 with respect to parameter k.   
           V(1) = -F2B/XPP/(Q2/X+B(i,2))                                
           V(2) =  F2B/XPP/(Q2/X+B(i,2))**2*(Q2+B(i,1))                 
           POL1 =  3.*B(i,5)*Y**2+4.*B(i,6)*Y**3+5.*B(i,7)*Y**4+        
     >             6.*B(i,8)*Y**5+7.*B(i,9)*Y**6                        
           V(3) = -F2B*POL1/POLY/(Q2/X+B(i,4))                          
           V(4) =  F2B*POL1/POLY/(Q2/X+B(i,4))**2*(Q2+B(i,3))           
           DO 10 j = 5,9                                                
10         V(j) =  F2B/POLY*Y**(j-2)                                    
           STB = 0.                                                     
           DO 11 j = 1,9                                                
           DO 11 k = 1,9                                                
11         STB = STB + L(i,j,k)*V(j)*V(k)                               
           STB = SQRT(STB)*BINDING                                      
                                                                        
      !LAMBDA12 MODEL NEXT:                                             
           Y    = 1.-X                                                  
           q    = LOG(Q2)                                               
           qth  = .2+3.2*X                                              
           dq   = q-qth                                                 
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                  
     >            C(i,4)*Y**6+C(i,5)*Y**7                               
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                   
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq       
           IF (q.GT.qth) QUAD = 0.                                      
           SCB  = (1.+LIN+QUAD)                                         
           F2L  = F2th*SCB*BINDING                                      
          !-----------------------------------------------------------  
          !Z(k) is the derivative of F_2 with respect to parameter k.   
           DO 20 j = 1,5                                                
20         Z(j) = SCB*Y**(j+2)                                          
           Z(6) = 0.                                                    
           IF (q.LT.qth) Z(6) = F2th*dq**2                              
           Z(7) = Z(6)*X                                                
           Z(8) = Z(6)*X**2                                             
           DO 21 j = 9,12                                               
21         Z(j) = F2th*X**(j-9)*dq                                      
           STL = 0.                                                     
           DO 22 j = 1,12                                               
           DO 22 k = 1,12                                               
22         STL = STL + M(i,j,k)*Z(j)*Z(k)                               
           STL = SQRT(STL)*BINDING                                      
                                                                        
          !U(k) is the derivative of slope with respect to parameter k. 
           SLOPE= F2th*LIN/dq*BINDING                                   
           DO 30 j = 1,5                                                
30         U(j) = LIN/dq*Y**(j+2)                                       
           DO 31 j = 6,8                                                
31         U(j) = 0.                                                    
           DO 32 j = 9,12                                               
32         U(j) = Z(j)/dq                                               
           DSLOPE = 0.                                                  
           DO 33 j = 1,12                                               
           DO 33 k = 1,12                                               
33         DSLOPE = DSLOPE + M(i,j,k)*U(j)*U(k)                         
           DSLOPE = SQRT(DSLOPE)                                        
      !---------------------------------------------------------------- 
                                                                        
      F2 = 0.                                                           
      ST = 0.                                                           
      IF (MODEL.EQ. 9) THEN                                             
           F2 = F2B                                                     
           ST = STB                                                     
      ELSEIF (MODEL.EQ.12) THEN                                         
           F2 = F2L                                                     
           ST = STL                                                     
      ELSE                                                              
           WRITE(*,'('' F1990: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')    
      ENDIF                                                             
      SY = ABS(F2B-F2L)                                                 
                                                                        
      GOODFIT = .TRUE.                                                  
      !The following cuts define the region of applicability of F1990.  
      !In order they are:                                               
      !     [radiative corrections convergence criteria in Q2] .and.    
      !     [radiative corrections convergence criteria in x]  .and.    
      !     [stay out of resonance region, W2.ge.3.0]          .and.    
      !     [limitation imposed by maximum beam energy].                
      IF ((Q2.LT..566).OR.(X.LT..062).OR.(X.LT.Q2/(2.*.93828*21.))      
     >   .OR.(X.GT.1./((3.-.93828**2)/Q2+1.)))     THEN                 
                                                                        
C         WRITE(*,'('' WARNING[F1990]: OUTSIDE RECOMMENDED RANGE.'')')  
          GOODFIT=.FALSE.                                               
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               

C====================================================================== 

!---------------------------------------------------------------------
! -*-Mode: Fortran; compile-command: "f77 -o f2nmc.o -c f2nmc.f"; -*-
      Subroutine F2NMC_new(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the NMC parametrisation   
! in CERN_PPE/95-138  Sept 4, 1995  Proton and Deuteron F2 Structure Functions 
! in deep inelastic muon scattering   P. Arneodo et. al.               
!  Published in Phys.Lett.B364:107-115,1995 
!   e-Print Archive: hep-ph/9509406 
!                   
      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi ,FNP_NMC
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2NMC_NEW_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       F2 = F2 * FNP_NMC(X4,QSQ4)
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF

      RETURN
      END

      SUBROUTINE F2NMC_NEW_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)
     >  /-0.02778,  2.926,   1.0362, -1.840, 8.123, -13.074, 6.215, 
     >   -0.04858,  2.863,    .8367, -2.532, 9.145, -12.504, 5.473/ 
      real*8 b(4,2)/ 0.285,   -2.694,   0.0188,  0.0274,  
     >              -0.008,   -2.227,   0.0551,  0.0570/
      real*8 c(4,2)/-1.413,    9.366, -37.79,   47.10, 
     >              -1.509,    8.553, -31.20,   39.98/

!lower limits
      real*8 al(7,2)
     >  /-0.01705,  2.851,   0.8213, -1.156, 6.836, -11.681, 5.645, 
     >   -0.02732,  2.676,    .3966, -0.608, 4.946,  -7.994, 3.686/ 
      real*8 bl(4,2)/ 0.325,   -2.767,   0.0148,  0.0226,  
     >                0.141,   -2.464,   0.0299,  0.0396/
      real*8 cl(4,2)/-1.542,   10.549, -40.81,   49.12, 
     >               -2.128,   14.378, -47.76,   53.63/

!upper limits
      real*8 au(7,2)
     >  /-0.05711,  2.887,   0.9980, -1.758, 7.890, -12.696, 5.992, 
     >    -0.04715,  2.814,    .7286, -2.151, 8.662, -12.258, 5.452/ 
      real*8 bu(4,2)/ 0.247,   -2.611,   0.0243,  0.0307,  
     >               -0.048,  -2.114,   0.0672,  0.0677/
      real*8 cu(4,2)/-1.348,    8.548, -35.01,   44.43, 
     >               -1.517,    9.515, -34.94,   44.42/
      

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0. !new on 3/15/00  SER
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28
      
      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               


      SUBROUTINE R1998(X,Q2,R,DR,GOODFIT)
     +  
                                                                        
!----------------------------------------------------------------       
! X      : Bjorken x                                                    
! Q2     : Q squared in (GeV/c)**2                                      
! R      :                                                              
! DR     : Absolute error on R                                          
! GOODFIT:  = .TRUE. if the X,Q2 are within the range of the fit.       
!-------------------------------------------------------------------    
! Model for R, based on a fit to world R measurements. Fit performed by 
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details   
! see Reference.                                                        
!                                                                       
! Three models are used, each model has three free parameters.  The     
! functional forms of the models are phenomenological and somewhat      
! contrived.  Each model fits the data very well, and the average of    
! the fits is returned.  The standard deviation of the fit values is    
! used to estimate the systematic uncertainty due to model dependence.  
!                                                                       
! Statistical uncertainties due to fluctuations in measured values have 
! have been studied extensively.  A parametrization of the statistical  
! uncertainty of R1990 is presented in FUNCTION DR1990.                 
!                                                                       
!                                                                       
! Each model fits very well.  As each model has its own strong points   
! and drawbacks, R1998 returns the average of the models.  The          
! chisquare for each fit (237 points with 6 parameters) are:  
!          ALL DATA  #PTS=237         |     X<  0.07 #PTS= 28
! FIT  #PARAM CHISQ  CHISQ/DF PROB(%) |  CHISQ  CHISQ/DF   PROB(%)
! R1998         217.4   0.94    73.1        28.7   1.06    37.5
! R1998a   6    219.4   0.95    69.8        28.8   1.07    37.1
! R1998b   6    217.7   0.94    72.6        27.8   1.03    42.2
! R1998c   6    221.9   0.96    65.5        30.9   1.15    27.4                

!                                                                       
! This subroutine returns reasonable values for R for all x and for all 
! Q2 greater than or equal to .3 GeV.                                   
!                                                                       
! The uncertainty in R originates in three sources:                     
!                                                                       
!     D1 = uncertainty in R due to statistical fluctuations of the data 
!          as reflected in the error of the fit. 
!          It is parameterized in FUNCTION DR1990, for details see     
!          Reference.                                                   
!                                                                       
!     D2 = uncertainty in R due to possible model dependence, approxi-  
!          mated by the variance between the models.                    
!                                                                       
!     D3 = uncertainty in R due to possible epsilon dependent errors    
!          in the radiative corrections, taken to be +/- .025.  See     
!          theses (mine or Dasu's) for details.  This is copied from R1990                       
!                                                                       
! and the total error is returned by the program:                       
!                                                                       
!     DR = is the total uncertainty in R, DR = sqrt(D1+2D) 7
!          DR is my best estimate of how well we have measured R.  At   
!          high Q2, where R is small, DR is typically larger than R.  If
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,  
!          you might wish to assume DR = 0 at very high Q2.             
!                                                                       
! NOTE:    In many applications, for example the extraction of F2 from  
!          measured cross section, you do not want the full error in R  
!          given by DR.  Rather, you will want to use only the D1 and D2
!          contributions, and the D3 contribution from radiative        
!          corrections propogates complexely into F2.  For more informa-
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or     
!          for explicite detail, see Reference.                         
!                                                                       
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!***** modified 10/28/01 pyb to give R at  max(0.5,Q**2)

      IMPLICIT NONE                                                     
      REAL QP2,FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1998,X,Q2
      REAL Q2_SAVE

!!** Note, in S. Rock version, this is 0.2
      REAL Q2_LIMIT/.5/
      REAL A(6) /4.8520E-02,  5.4704E-01,  2.0621E+00,
     >          -3.8036E-01,  5.0896E-01, -2.8548E-02/
      REAL B(6) /4.8051E-02,  6.1130E-01, -3.5081E-01, 
     >          -4.6076E-01,  7.1697E-01, -3.1726E-02/
      REAL C(6) /5.7654E-02,  4.6441E-01,  1.8288E+00,
     >           1.2371E+01, -4.3104E+01,  4.1741E+01/


      LOGICAL GOODFIT                                                   
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      Q2_SAVE = Q2
      IF(Q2.LT.Q2_LIMIT) THEN   ! added 10/28/01 to match Blok et al
       Q2 = Q2_LIMIT
      ENDIF
                                 
      FAC = 1.+12.*Q2/(Q2+1.)*(.125**2 /(.125**2+X**2))
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!      


!Model A
      QP2 = (A(2)/(Q2**4+A(3)**4)**.25) * (1.+A(4)*X +A(5)*X**2) !1.07   .030
      R_A   = A(1)*RLOG +QP2*X**A(6)                                           
!Model B
      QP2 = (1.+B(4)*X+B(5)*X**2)*(B(2)/Q2 +B(3)/(Q2**2+.09)) !1.06   .042
      R_B   =  B(1)*RLOG+QP2*X**B(6)   
!Model C
      Q2thr =C(4)*(X) +c(5)*x**2 +c(6)*X**3  
      QP2 =  C(2)/SQRT((Q2-Q2thr)**2+C(3)**2) 
      R_C   =  C(1)*RLOG+QP2    

      R     = (R_A+R_B+R_C)/3.                                          

      D1    = DR1998(X,Q2)                                              

      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)               
      D3    = .023*(1.+.5*R)                                            
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                       


      DR    = SQRT(D1**2+D2**2+D3**2)                                   
                                                                        
      GOODFIT = .TRUE.                                                  
      IF ((X.LT.0.02).OR.(Q2.LT.0.3)) GOODFIT = .FALSE.                                   

! Restore Q2
      Q2 = Q2_SAVE

!*** commented out 10/28/01
!      IF(Q2.LT.Q2_LIMIT) THEN   ! added 11/15 to avoid low Q2 problems caused by RLOG
!       R = Q2/Q2_LIMIT * R
!      ENDIF
      

      RETURN                                                            
      END                                                               

C====================================================================== 

      SUBROUTINE D2MODEL_IOANA(QSQ,WSQ,W1,W2)
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
*****************************************************************************
        IMPLICIT NONE

        REAL*4  WSQ, QSQ
        REAL*4  R_NRES, SIG_RES(2), SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, 
     >          DIPOLE, K, NU,TAU, PI, AM, ALPHA, W1, 
     >          W2,CONV,sig_roper(2)
	integer i
	real*4  xval(37)
	logical goroper/.false./
*
        COMMON/ioana/xval
*       
*
* N.I. ...
*

        DATA PI/3.14159265/, AM/.93828/, ALPHA/0.00729735/,
     >          CONV/0.0025767/

*
! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.17) return
        IF(WSQ.LT.1.17.OR.WSQ.GT.5.OR.
     >    QSQ.GT.10.0) THEN
           do i=1,2	
              SIG_NRES(i) = 0.0
              SIG_RES1(i) = 0.0
              SIG_RES2(i) = 0.0
              SIG_ROPER(i)= 0.0 
           enddo
!          WRITE(*,*)'H2MODEL_IOANA called outside of kinematic range'
          RETURN
        ENDIF

!  returns transverse cross sections in units of
! microbarns/(dipole FF)**2
        CALL i_d2_model(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     &                  SIG_ROPER,goroper,xval)
*        write(*,*)'1',SIG_NRES
        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.
        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
*        write(*,*)'1',SIG_NRES(1),DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1)+sig_roper(1))*
     &                 DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES(1)
        SIGL = R_NRES*SIG_NRES(1)
*	write(*,*) 'I am here'

!        write(33,*)SIG_NRES(1),SIG_RES(1),R_NRES,DIPOLE,sigt
! 33     format(F15.3,2x,F15.3,2x,F15.3,2x,F15.3,2x,F15.3)      

! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)
*        write(*,*)'from program',w1,w2 

!*** PYB DIVIDE BY 2 TO GET PER NUCLEON
        W1=W1/2.
        W2=W2/2.

        RETURN
        END


CCC  Version 031606  -  Author:  M.E. Christy                        CCC
C*** changed to 7/7/06 version (see below)
CCC  Subroutine to get Transvese and Longitudinal eP cross sections  CCC 
CCC  from fits to L/T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


      SUBROUTINE CHRISTY31606(W2,Q2,F1,R)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval1(50),xvall(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1,fl,r
      integer i,npts,sf
 
      mp = .9382727
      mp2 = mp*mp   
      pi = 3.141593
      alpha = 1./137.036

      data xval1 / 
c Iteration of July 26, 2006
     & 0.12292E+01,0.15000E+01,0.15136E+01,0.16930E+01,0.16167E+01,
     & 0.14115E+01,0.14189E+00,0.23400E+00,0.78397E-01,0.94916E-01,
     & 0.16490E+00,0.23595E+00,0.86273E+02,0.10226E+02,0.36352E+00,
     & 0.38246E+01,0.13157E+02,0.69880E+00,0.93290E-01,0.25211E+01,
     & 0.17480E+02,0.45287E+01,0.25665E+00,0.30034E+01,0.93534E+01,
     & 0.10069E+02,0.27692E+00,0.29360E+01,0.23039E+02,0.49782E+01,
     & 0.19785E+00,0.48402E+01,0.15297E+01,0.23360E+02,0.26987E+00,
     & 0.25097E+01,0.27634E+03,0.90974E-01,0.12214E+01,0.13983E+00,
     & -.29210E+03,0.75702E+00,0.26522E+01,-.37007E+00,-.43681E-03,
     & 0.82147E-01,0.20200E+01,0.59499E+00,0.34074E+02,0.93053E-01 /


      data xvalL/
C  Simona - July 28, 2006
     & 0.12438E+01,0.14944E+01,0.12778E+01,0.17327E+01,0.16700E+01,
     & 0.13921E+01,0.70781E-01,0.45000E-01,0.53394E-01,0.75084E-01,
     & 0.60000E-01,0.45810E+00,0.28904E+02,0.73691E+04,0.20492E+05,
     & 0.00000E+00,0.51975E+01,0.49889E+05,0.64123E+05,0.00000E+00,
     & 0.13818E+05,0.10000E+05,0.34644E+04,0.79746E+05,0.41526E+04,
     & 0.38928E+06,0.93410E+05,0.43281E+04,0.00000E+00,0.92972E+04,
     & 0.87089E+04,0.46171E+05,0.17525E+02,0.70074E+04,0.75226E+04,
     & 0.79992E+00,0.70392E+03,0.13460E+01,0.00000E+00,0.34373E+01,
     & 0.19070E+01,0.22171E+00,0.95006E+03,0.74309E+01,0.13754E+00,
     & 0.15949E+01,0.00000E+00,0.00000E+00,0.66140E+00,0.38369E+00 /

      F1=0.
      R=0.
      if(w2.lt.1.155) return

      xb = q2 / ( w2 + q2 - mp2)
      if(xb.le.0.0) return

      call resmod316(1,w2,q2,xval1,sigT)
      call resmod316(2,w2,q2,xvalL,sigL)

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      FL = sigL*2.*xb*(w2-mp2)/8./pi/pi/alpha/0.3894e3
      R = FL/ (2.D0 * XB * F1)

   
      end
C====================================================================== 


      SUBROUTINE RESMOD316(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,2),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,2),x0(7),dip,dip2
      REAL*8 sig_res,sig_4L,sigtemp,slope,q2low,dq2,t,xpr
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2
      dq2 = 0.05

      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (mp + mpi)
      wr = wdif/w

      br(1,1) = 1.0     !!! single pion branching ratios
      br(2,1) = 0.5
      br(3,1) = 0.65
      br(4,1) = 0.65
      br(5,1) = 0.4
      br(6,1) = 0.65
      br(7,1) = 0.6

      if(sf.EQ.2) then 
        br(6,1) = xval(48)
        br(2,1) = xval(49)
      endif 

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  ? 4th resonance region

      do i=1,7
        x0(i) = 0.165
      enddo
      x0(4) = 0.6

      do i=1,7
        br(i,2) = 1.-br(i,1)
      enddo
    

      if(sf.EQ.1) then
        q2low = 0.00
      else
        q2low = 0.1
      endif 

      if(q2.LT.q2low) then
        lowq2 = .true.
        lmax = 2
      endif

c      write(6,*) q2,lmax

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = q2low
        elseif(l.EQ.2.AND.lowq2) then
          q2 = q2low + dq2
        endif

        dip = 1./(1.+q2/0.71)**2             !!!  Dipole parameterization  !!!
        dip2 = dip*dip

        xb = q2/(q2+w2-mp2)
        xpr = 1.+(w2-(mp+mpi)**2)/(q2+xval(50))
        xpr = 1./xpr
c        t = log(log((q2+xval(50))/0.330**2)/log(xval(50)/0.330**2))

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6              !!!  Read in resonance masses     !!!
          num = num + 1
          mass(i) = xval(i)
        enddo
        do i=1,6              !!!  Read in resonance widths     !!!
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

        if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
          mass(7) = xval(41)
          intwidth(7) = xval(42)
          width(7) = intwidth(7)
        else
          mass(7) = xval(47)
          intwidth(7) = xval(48)
          width(7) = intwidth(7) 
        endif

        do i=1,7
          kr(i) = (mass(i)**2-mp2)/2./mp
          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

          pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)

          pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+1.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))**ang(i)

          if(i.EQ.2) then
            pwid(i,2) =  intwidth(2)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
          endif 

          pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

          pgam(i) = intwidth(i)*pgam(i)

          width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)

        enddo
 

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo

          if(sf.EQ.1) then

            if(i.eq.6) height(i) = rescoef(i,1)/
     &        (1.+ q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)


             height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &          (1.+q2/0.71)**rescoef(i,4)

          else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif

        enddo

CCC    End resonance Q^2 dependence calculations   CCC

     
        do i=1,3               !!!  Non-Res coefficients  !!!
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo


        if(sf.EQ.2) then      !!!  4th resonance region  !!!
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
        else
          height(7) = xval(49)*dip2 
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC

        sig_res = 0.0

        do i=1,7
          sigr(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)
          sigr(i) = sigr(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
          sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
          sig_res = sig_res + sigr(i)   
        enddo


CCC    Finish resonances / start non-res background calculation   CCC

 
        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
          enddo

          sig_nr = sig_nr*xpr


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xpr)**2.*q2/(1.+nr_coef(i,3)*q2)
     &       /(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif


        sig = sig_res + sig_nr

          
        if(L.EQ.1) sigtemp = sig  

      enddo
       

CCC   Now extrapolate sig_L linearly to zero for Q^2 less than q2min   CCC

      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/dq2
          sig = sigtemp + slope*(q2low-q2)
        else
          slope = sig/q2low
          sig = sig - slope*(q2low-q2)     
        endif
      endif


 1000  format(9f12.5)

      RETURN 
      END 


CCC  Version 050521  -  Author:  M.E. Christy   modified by S. Malace CCC
CCC  Subroutine to get Transvese eD cross sections  CCC 
CCC  from fits to T cross sections.  The subroutine resmod.f is    CCC
CCC  required.  Units are in ub/Sr/Gev.                              CCC


C S. Malace 7/22/06 to get D2 

      SUBROUTINE rescssim(W2,Q2,F1)

      IMPLICIT NONE

      real*8 w2,q2,fn,fnerr,xval11(50),xval12(50),w2inc,temp(4)
      real*8 mp,mp2,pi,alpha,xb,sigt,sigl,f1
      integer i,npts,sf
 
c following is iteration 3 starting with Bodek. /rc/simonaD2SimFit6.out
C  Simona D2 fit to (simona data) + (ioana-delta)
      data xval11/
     & 0.118651E+01,0.15247E+01,0.14354E+01,0.213435E+01,0.17131E+01,
     & 0.29094E+01,0.30152E+00,0.14436E+00,0.400969E+00,0.55492E+00,
     & 0.29446E+00,0.97992E+00,0.460104E+04,0.17992E+04,0.25114E+04,
     & -.37399E+03,0.24733E+04,0.75449E+03,0.8939804E+05,0.3523602E+05,
     & 0.30091E+02,0.665424E+05,0.790445E+05,-.42655E+04,0.99998E+02,
     & 0.9217486E+06,0.2988304E+06,0.869219E+05,0.512361E+05,0.9782E+03,
     & 0.13294E+02,0.38278E+01,0.10849E+04,0.100891E+05,0.4726618E+05,
     & -.46058E+04,0.36094E+03,0.17463E+01,0.130750E+01,0.670575E-01,
     & 0.969567E+05,0.361745E+02,0.335147E+01,-.86421E+00,0.21919E-02,
     & 0.28124E+00,0.135992E+01,0.65969E+00,0.16259E+03,0.26693E+01 /


      data xval12/


C  Simona D2 fit to ioana data
     & 0.12271E+01,0.155191E+01,0.148304E+01,0.169691E+01,0.17799E+01,
     & 0.19778E+01,0.125836E+00,0.25999E+00,0.267779E+00,0.208961E+00,
     & 0.277916E+00,0.35815E+00,0.3517854E+05,0.342156E+04,0.173279E+05,
     & 0.665358E+04,0.235974E+04,0.851821E+02,0.717906E+05,0.821708E+04,
     & 0.738373E+03,0.245138E+04,0.8417701E+05,0.1622773E+05,0.8442E+01,
     & 0.872436E+06,0.4964696E+06,0.57930E+05,0.8911091E+06,0.12752E+04,
     & 0.10143E+04,0.319242E+04,0.831498E+04,0.932843E+04,0.3098729E+05,
     & 0.38488944E+07,0.38934E+03,0.17783E+00,0.151086E+01,0.9301E-01,
     & -.889461E+04,0.13923E+02,0.454747E+01,-.65654E-01,-.43802E-02,
     & 0.939828E+03,0.13554E+01,0.18729E+00,0.38683E+02,0.52966E-01 /

       if(q2.lt.4.5.and.w2.lt.2.and.q2.gt.1.7) then

         call resmodsim(1,w2,q2,xval12,sigt)

       else

         call resmodsim(1,w2,q2,xval11,sigt)

      endif  

      mp = .9382727
      mp2 = mp*mp
      pi = 3.141593
      alpha = 1./137.036

      F1 = sigT*(w2-mp2)/8./pi/pi/alpha/0.3894e3

      return
      end



CCC  Version 061105  -  Author:  M.E. Christy                        CCC
CCC  Fit form is empirical.  Please do not try to interpret physics  CCC
CCC  from it.  This routine needs the parameter files f1parms.dat    CCC
CCC  and flparms.dat.  Units are ub/Sr/Gev.                          CCC

C THIS IS FOR D2 from S. MALACE 7/22/07
             
      SUBROUTINE RESMODSIM(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),wdif,wdif2,wr,sig_nr,sig_4,q2temp
      REAL*8 mpi,meta,intwidth(6),k,kcm,kcmr(6),ppicm,ppi2cm,petacm
      REAL*8 ppicmr(6),ppi2cmr(6),petacmr(6),epicmr(6),epi2cmr(6)
      REAL*8 eetacmr(6),epicm,epi2cm,eetacm,br_21_1,br_21_2
      REAL*8 sig_res,sig_4L,sigtemp,slope
      INTEGER i,j,l,lmax,num,sf
      LOGICAL lowq2

      lowq2 = .false.
      lmax = 1
      q2temp = q2

      mp = 0.9382727
      mpi = 0.136
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif = w - (0.937272 + 0.136)
      wr = wdif/w


c      if(q2.LT.0.3.AND.sf.EQ.1) then
      if(q2.LT.0.15) then
        lowq2 = .true.
        lmax = 2
      endif

      do l=1,lmax

               
        if(l.EQ.1.AND.lowq2) then
          q2 = 0.15
        elseif(l.EQ.2.AND.lowq2) then
          q2 = 0.25
        endif

        xb = q2/(q2+w2-mp2)
        xth(1) = (q2 + xval(50))/(w2-mp2-0.136+q2)


CCC  Calculate kinematics needed for threshold Relativistic B-W  CCC

        k = (w2 - mp2)/2./mp
        kcm = (w2-mp2)/2./w

        epicm = (W2 + mpi**2 -mp2 )/2./w
        ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
        epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
        ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
        eetacm = (W2 + meta*meta -mp2 )/2./w
        petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

        num = 0

        do i=1,6
          num = num + 1
          mass(i) = xval(i)

          kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
          epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
          ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
          epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
          ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
          eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
          petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

       enddo
 
        do i=1,6
          num = num + 1
          intwidth(i) = xval(num)
          width(i) = intwidth(i)
        enddo

           
c      write(6,*) "1:  ",num

        do i=1,6
          do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
          enddo
          if(sf.EQ.1) then

            height(i) = rescoef(i,1)/
     &        (1.+q2/rescoef(i,2))**(rescoef(i,3)+rescoef(i,4)*q2)

            if(i.EQ.1) height(i) = 3.0*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.2) height(i) = 1.4*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 
            if(i.EQ.5) height(i) = 0.3*
     &         (1.+rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)/
     &         (1.+q2/rescoef(i,3))**rescoef(i,4))/(1.+q2/0.71)**2. 

          else
c            height(i) = rescoef(i,1)*
c     &            (1.+q2/rescoef(i,2))**(-1.*rescoef(i,3))
            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,4)*q2)/
     &            (1.+q2/rescoef(i,2))**rescoef(i,3)     !!!  Galster Shape  !!!             

          endif
          if(height(i).LT.0) height(i) = 0. 

        enddo
     

        do i=1,3
          do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
          enddo
        enddo

        if(sf.EQ.2) then      !!!  Put in Roper  !!!
          mass(7) = xval(41)
          width(7) = xval(42)
          height(7) = xval(43)*q2/(1.+xval(44)*q2)/
     &            (1.+q2/xval(45))**xval(46)     !!!  Galster Shape  !!!
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        else
          mass(7) = xval(47)
          width(7) = xval(48)
          height(7) = xval(49)/(1.+q2/0.71)**3.    
          sig_4L   = height(7)*width(7)/((W-mass(7))**2. 
     &               + 0.25*width(7)*width(7))  
        endif


CC   Calculate Breit-Wigners for the 3 resonance regions   CC


        sig_32 = width(5)/((W-mass(5))**2. 
     &               + 0.25*width(5)*width(5))
        sig_4   = width(6)/((W-mass(6))**2. 
     &               + 0.25*width(6)*width(6))

        br_21_1 = 0.5
        br_21_2 = 0.5
        if(sf.EQ.2) then
          br_21_1 = xval(48)
          br_21_2 = 1.- br_21_1
        endif

        width(1)=intwidth(1)*ppicm/ppicmr(1)
        width(2)=intwidth(2)*(br_21_1*ppicm/ppicmr(2)
     &            +br_21_2*petacm/petacmr(2))
        width(3)=intwidth(3)*(0.5*ppicm/ppicmr(3)+0.5*ppi2cm/ppi2cmr(3))
        width(4)=intwidth(4)*
     &                     (0.65*ppicm/ppicmr(4)+0.35*ppi2cm/ppi2cmr(4))

c      write(6,*) ppicm,ppicmr(3),petacm,petacmr(3),intwidth(3)

        sig_del = ppicm/kcm/((W2 - mass(1)**2.)**2. 
     &              + (mass(1)*width(1))**2.)
        sig_21 =  (0.5*ppicm+0.5*petacm)/kcm/
     &           ((W2 - mass(2)**2.)**2. + (mass(2)*width(2))**2.)
        sig_22 =  (0.5*ppicm+0.5*ppi2cm)/2./kcm/
     &           ((W2 - mass(3)**2.)**2. + (mass(3)*width(3))**2.)
        sig_31 =  (0.65*ppicm+0.35*ppi2cm)/2./kcm/
     &           ((W2 - mass(4)**2.)**2. + (mass(4)*width(4))**2.)
        if(sf.EQ.2) then
          width(5)=intwidth(5)*
     &     (xval(47)*petacm/petacmr(5)+(1.-xval(5))*ppi2cm/ppi2cmr(5))

          sig_32 =  (xval(47)*petacm+(1.-xval(47))*ppi2cm)/2./kcm/
     &           ((W2 - mass(5)**2.)**2. + (mass(5)*width(5))**2.)

        endif
        

        sig_del = height(1)*sig_del
        sig_21 = height(2)*sig_21
        sig_22 = height(3)*sig_22
        sig_31 = height(4)*sig_31
        sig_32 = height(5)*sig_32
        sig_4   = height(6)*sig_4

        sig_nr = 0.

        if(sf.EQ.1) then

          do i=1,2  
            sig_nr = sig_nr +(nr_coef(i,1)*(wdif)**(float(2*i+1)/2.))
     &       /(q2+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2*q2)
          enddo

          sig_nr = sig_nr*xb


        elseif(sf.EQ.2) then
          do i=1,1 
            sig_nr = sig_nr +(nr_coef(i,1)*wdif**(float(i)/2.)*
     &       (1.-xth(1))**2.+nr_coef(i,3)*wdif**(float(3*i-1)/2)
     &       *(1.-xth(1))**3.)/(1.+ q2/nr_coef(i,2))**nr_coef(i,4)/w2

          enddo                           
        endif

        sig_res = sig_del + sig_21 + sig_22 + sig_31 + sig_32 + sig_4
 
        sig_res = sig_res + sig_4L

        if(sf.EQ.2) then
c          sig_res = sig_res + sig_4L
          sig_nr = sig_nr*q2/(1.+xval(49)*q2)
        endif

        sig = sig_res + sig_nr

c        sig = sig_res  

        if(w2.LE.1.16.OR.sig.LT.0) sig = 0.d0
          
        if(L.EQ.1) sigtemp = sig  

      enddo
       
      q2 = q2temp

      if(lowq2) then
        if(sf.EQ.1) then
          slope = (sigtemp - sig)/0.1
          sig = sigtemp + slope*(0.15-q2)
        else
          slope = sig/0.15
          sig = sig - slope*(0.15-q2)     
        endif
      endif
c      if(lowq2) write(6,*) q2, sig,sigtemp,slope

c      if(sf.eq.1.AND.q2.GT.5) write(6,1000) sig,sig_res,sig_nr 

 1000  format(9f12.5)

      RETURN 
      END 

!-----------------------------------------------------------------------
      SUBROUTINE Q_E_VANORDEN(Q2_ELAS_GEV,E_GEV,SUPPRESSION)
!---------------------------------------------------------------------------
C   This program compute the quasi-elastic cross section
C   based on the Van Orden calculation using the fermi gas model
!  input energy now in GeV
! It returns the Suppression factor for the quasi-elastic peak.
!-----------------------------------------------------------------------------

      IMPLICIT NONE
      REAL E_GEV,Z,A,KF_GEV,SUPPRESSION
      REAL SUM,AP,FMT,FACTOR,N,MT,TH,SINSQ,FMOTT,WMOTT,TANSQ,EP_ELAS,
     >  EF,RLR,RTR,SRL,RATIO,W,WBAR,QV,Q2BAR,F2,XETA,ZETA,T1,T2,E1,
     >  E,SL,ST,RS,RSSIG,XSECT,CROSS,QT,C,Q2_ELAS,E2,W0,W1,Q2,EINC,KF,
     >  KF3,Z1,A1,QV2,WSQ,MT4PI,MP2,Q2_ELAS_GEV
      INTEGER IOPT,I0,I00,I,IOPT1
      REAL MUP2/7.784/,MUN2/3.65/
      REAL EBAR/1./,WDEL/2./ !$$ are these OK?
      REAL  ALPH/.0072972/, PI/3.1415927/,AMASS/931.5/,MP/939./
      LOGICAL FIRST/.TRUE./
!-----------------------------------------------------------------------------------------
! In this notation, W=nu
!---------------------------------------------------------------------------------------     
      EINC = 1000.*E_GEV
      Q2_ELAS = 1.E6 * Q2_ELAS_GEV

C$$      OPEN(7,FILE='Q_E_VANORDEN2')
      SUM=0.

      IF(IOPT.EQ.1)THEN
 1     EP_ELAS = EINC - Q2_ELAS/(2.*MP)
       IF(EP_ELAS.GT.0) THEN
        SINSQ = Q2_ELAS/(4.*EINC*EP_ELAS)
        TH = 2.*ASIN(SQRT(SINSQ))
       ELSE
        EINC = EINC +5.
        GO TO 1
       ENDIF

       FMOTT=ALPH*ALPH*COS(TH/2.)**2/4./SINSQ/SINSQ*197.3**2*1.E-26
       WMOTT=FMOTT/EINC**2
       TANSQ=TAN(TH/2.)**2
       QT = sqrt(Q2_ELAS**2/(4.*MP2) + Q2_ELAS)
       IF(QT.GT. 2.*KF) THEN
        SUPPRESSION=1.
        RETURN
       ENDIF
       W0=MAX(EINC-(EP_ELAS+2.*KF),2.)
       W1=MAX(EINC-(EP_ELAS-2.*KF),0.)
C$$       WRITE(7,'(/''E_INC,THET,WDEL,FERMI_MOM='',1P9E10.2)')
C$$     >   EINC,THETA,WDEL,KF
C$$       WRITE(7,'(''EBAR,Z,A,Q2_ELAS,W0,W1='',1P10E9.2)')
C$$     >  EBAR,Z,A,Q2_ELAS,w0,w1
C$$       WRITE(7,'(/''        NU          SIGMA(nb)'')')
      ENDIF
      
      RLR=0.
      RTR=0.
      SRL=0.
C      QV=0.
      RATIO=1.
C
      W=W0
      I00=W0/WDEL
      I0=W1/WDEL
      DO 17 I=I00,I0
       WBAR=W-EBAR
       IF(WBAR.LE.0.)GO TO 15
       IF(IOPT.EQ.1)Q2=4.*EINC*(EINC-W)*SINSQ
c       WSQ = W**2
       QV2 = Q2+WSQ
       QV=SQRT(QV2)
       Q2BAR=QV2-WBAR**2
       E1=SQRT(QV2*MP2/Q2BAR+QV2/4.)-WBAR/2.
       IF(EF.LT.E1)GO TO 18     ! do not calculate 
       RATIO=Q2/(Q2+WSQ)
       F2=1./(1.+Q2/855.**2)**4
       XETA=Q2/(4.*MP2)
       ZETA=1.+XETA
!       T1=F2*Q2/2.*(((1.+2.79*XETA)/ZETA+(1.79/ZETA))**2+N/Z*3.65)
!         T1 = 2Mp**2 * DIPOLE *Tau*(MuP2 + N/Z * MuN2) = Tau*(Gmp**2 + Gmn**2 )
!        = 2Mp*Tau*(F1+(Mu-1)F2)**2 where
!          F1= DIPOLE*(1+TAU*Mu)/(1+Tau)    F2= DIPOLE/(1+Tau)
!       T2=2.*MP22*(((1.+2.79*XETA)/ZETA)**2+XETA*((1.79/ZETA)**2+
!     >  N/Z*1.91**2))*F2  !$$$*** I think the neutron term should be divided by ZETA
!         T2=2Mp**2 *DIPOLE* (Gep +Tau*Gmp)/(1+Tau)  + neutron
!           =2Mp**2 * (F1**2 +Tau*(MuP-1)F2**2)
     
! Below is Steve's Redoing
       T1 =F2*Q2/2.*(MUP2 +N/Z*MUN2)
       T2 =2.*MP2*F2*
     >   ( (1.+MUP2*XETA)/ZETA +N/Z* (0.+MUN2*XETA)/ZETA)
       E2=EF-WBAR
       E=E1
       IF(E2.GT.E1)E=E2

       RLR=(.75*Z/(KF3*QV))*(T2*((EF**3-E**3)/3.+W*(EF**2-E**2)/2.
     >  +WSQ*(EF-E)/4.)/MP2-QV2*T1*(EF-E)/Q2)

       RTR=(.75*Z/(KF3*QV))*(2.*T1*(EF-E)+T2*(Q2BAR*(EF**3-E**3)/
     > (3.*QV2)+Q2BAR*WBAR*(EF**2-E**2)/(2.*QV2)-
     > (Q2BAR**2/(4.*QV2)
     > +MP2)*(EF-E))/MP2)
15     CONTINUE
       SL=RLR*MT4PI
       ST=RTR*MT4PI

       RS=RATIO*RATIO*SL+(0.5*RATIO+TANSQ)*ST
       RSSIG=WMOTT*RS
       XSECT=FACTOR*WMOTT*RS
       IF(SL.EQ.0.AND.ST.EQ.0.)GO TO 18
C       SRL=SRL+RLR
       IF(IOPT.EQ.1)THEN
        SUM = SUM + XSECT* WDEL
C$$        WRITE(7,35)W,XSECT
       ELSE
C$$        WRITE(7,36)W,RLR,RTR
       ENDIF
35     FORMAT(1X,1P1E15.2,1X,1P1E15.2)
36     FORMAT(1X,F6.2,2E15.8)
 18    W=W+WDEL
 17   CONTINUE

      F2=1./(1.+Q2_ELAS/855.**2)**4
      XETA=Q2_ELAS/(4.*MP2)
      ZETA=1.+XETA
 
      CROSS =WMOTT* F2* 1.E33 *
     >         ( Z*( (1.+MUP2*XETA)/ZETA +2*TANSQ*MUP2 *XETA)+
     >           N*( (0.+MUN2*XETA)/ZETA +2*TANSQ*MUN2 *XETA))
C$$   WRITE(7,'(''SIGMA: SUM='',1PE9.2,'';  PEAK='',1PE9.2)')SUM,CROSS
      SUPPRESSION = SUM/CROSS
! Tsai
      QT = sqrt(Q2_ELAS**2/(2.*938.28)**2 + Q2_ELAS)
      IF(QT.GT. 2.*KF) THEN
       C=1.
      ELSE
       C = .75*QT/KF *(1.-.083333*(QT/KF)**2)
      ENDIF
C$$      WRITE(7,'('' SUPPRESSION: TSAI='',F5.3,'':   Meziani='',F5.3)')
C$$     >   C,sum/cross
      RETURN


      ENTRY  Q_E_VANORDEN_INIT(Z1,A1,KF_GEV,IOPT1)  
       Z = Z1
       A = A1
       KF= 1000.*KF_GEV
       KF3 = KF**3
       IOPT=IOPT1
       AP=ALPH/PI
       FMT=A*AMASS
       FACTOR= 4.0*PI/FMT *1.E33
       N=A-Z
       MP2 = MP**2
       MT=931.5*(Z+N)
       EF=SQRT(KF**2+MP2)
       MT4PI = MT/4./PI
      RETURN
      END
                          
      SUBROUTINE F2PGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(proton) for E665                            C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2PHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2PHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2PHIW(QQ,NU)

c.....low W parametrization for proton
      LOWWF2 = F2LOWW(1,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END


!---------------------------------------------------------------------
      SUBROUTINE F2DGLO(DQ2,DX,DF2,DR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     A V Kotwal 6jan94  (kotwal@fnal.gov)                  C
C                                                           C
C ***  Trial F2(deuteron) for E665                          C
C ***  Constructed from various parametrizations of data    C
C                                                           C
C    arguments are double precision                         C
C inputs:                                                   C
C    DQ2 = Q2 in GeV^2                                      C
C    DX  = Xbj                                              C
C outputs:                                                  C
C    DF2 = F2 structure function                            C
C    DR  = R  (sigma_l/sigma_t)                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      DOUBLE PRECISION DQ2,DX,DF2,DR,DFX
      REAL QQ,R,ERR,F2LOWW,LOWWF2,NU,W,HIWF2,F2DHIW,FW,NEWF2,XX
      DOUBLE PRECISION PM,PM2,TWOPM
      PARAMETER (PM=0.93828D0)  ! proton mass
      PARAMETER (PM2=PM*PM,TWOPM=2.0D0*PM)
      EXTERNAL F2LOWW,F2DHIW
C
      NU = SNGL(DQ2/(DX*TWOPM))
      QQ = SNGL(DQ2)

      DFX = 1.0D0/DX - 1.0D0
      W  = SNGL ( DSQRT ( DFX*DQ2 + PM2 ) )

c.....smooth function of W switching at W=5
      IF (W.GT.6.0) THEN
        FW = 0.0
      ELSEIF (W.LT.4.0) THEN
        FW = 1.0
      ELSE
        FW = 1.0/(EXP(((W-5.0)/0.4)) + 1.0)
      ENDIF

c.....nmc+DOLA F2 for high W
      HIWF2 = F2DHIW(QQ,NU)

c.....low W parametrization for deuteron
      LOWWF2 = F2LOWW(2,QQ,NU)

c.....connect at W=5
      NEWF2 = LOWWF2*FW + (1.0-FW)*HIWF2
      DF2 = DBLE(NEWF2)

      XX = SNGL(DX)
      CALL R1990F(XX,QQ,R,ERR)
      DR=DBLE(R)
      RETURN
      END



      Subroutine F2SMC98(T,x4,qsq4,F2,err_lo,err_hi) 
! This subroutine returns a value for F2 from the SMC parametrisation
! in SMC9x/xx 3/27/98
! in Tables 12 and 13.


      IMPLICIT NONE                                                          
      real x4,qsq4,f2,err_lo,err_hi,errd_lo,errd_hi,FNP_NMC,F2D
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron
                                ! 3 = neutron

      IF(T.LE.2) THEN
       CALL F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)
      ELSE  !neutron
       CALL F2SMC98_DO(1,x4,qsq4,F2,err_lo,err_hi)  !proton
       CALL F2SMC98_DO(2,x4,qsq4,F2D,errd_lo,errd_hi)  !deuteron
       F2 = 2.*F2D- F2
       ERR_LO = ERR_LO * FNP_NMC(X4,QSQ4)
       ERR_HI = ERR_HI * FNP_NMC(X4,QSQ4)
      ENDIF
      RETURN
      END

      SUBROUTINE F2SMC98_DO(T,x4,qsq4,F2,err_lo,err_hi)                                         
      IMPLICIT NONE            
      real x4,qsq4,f2,f2l,f2u,err_lo,err_hi
      real*8 a(7,2)/
     >  -0.24997, 2.3963, 0.22896, 0.08498, 3.8608, -7.4143, 3.4342,  !p      
     >  -0.28151, 1.0115, 0.08415,-0.72973, 2.8647, -2.5328, 0.47477 /
      real*8 b(4,2)/ 0.11411, -2.2356, 0.03115, 0.02135,    !p
     >               0.20040, -2.5154, 0.02599, 0.01859/ 

      real*8 c(4,2)/ -1.4517, 8.4745, -34.379, 45.888,    !p
     >               -1.3569, 7.8938, -29.117, 37.657 /
              
!lower limits
      real*8 al(7,2)/
     >  -0.25196, 2.4297, 0.21913,  0.21630, 3.4645, -6.9887, 3.2771, !p  
     >  -0.28178, 1.1694, 0.09973, -0.85884, 3.4541, -3.3995, 0.86034/
      real*8 bl(4,2)/ 0.13074, -2.2465, 0.02995, 0.02039, 
     >                0.20865, -2.5475, 0.02429, 0.01760/  
      real*8 cl(4,2)/ -1.4715, 8.9108, -35.714, 47.338,
     >                -1.3513, 8.3602, -31.710, 41.106/

!upper limits
      real*8 au(7,2)/
     > -0.24810, 2.3632, 0.23643, -0.03241, 4.2268, -7.8120, 3.5822,    !p
     > -0.28047, 0.8217, 0.06904, -0.60191, 2.2618, -1.6507, 0.08909 /
      real*8 bu(4,2)/0.09734, -2.2254, 0.03239, 0.02233,
     >               0.18711, -2.4711, 0.02802, 0.01973  /
      real*8 cu(4,2)/ -1.4361, 8.1084, -33.306, 44.717,
     >                -1.3762, 7.6113, -27.267, 35.100  /

      real Lam/0.25/                                                    
      real*8 Qsqo/20.0/                                                   
      real*8 log_term                                                     
      real*8 A_x,B_x,C_x,dl,lam2
      real*8 x,qsq,z,z2,z3,z4,f28                                  
      integer T                 ! 1 = Proton                            
                                ! 2 = Deuteron                          
      logical first/.true./


                  
      if(first) then
        dl =dlog(qsqo/lam**2)  
        lam2 =lam**2
        first=.false.
      endif

      x=x4
      z=1.-x
      z2=z*z
      z3=z*Z2
      z4 =z3*z

      qsq=qsq4                                                      
      A_x=x**a(1,t)*z**a(2,t)*(a(3,t)+a(4,t)*z+a(5,t)*z2  
     >    +a(6,t)*z3+a(7,t)*z4)                             
      B_x=b(1,t)+b(2,t)*x+b(3,t)/(x+b(4,t))                             
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4                  
      Log_term = dlog(qsq/lam2)/dl                     
      if(Log_term.lt.0) Log_term=0.     !new on 3/15/00  SER      
      F28=A_x*(log_term)**B_x*(1+C_x/qsq)
      f2 = f28

      A_x=x**au(1,t)*z**au(2,t)*(au(3,t)+au(4,t)*z+au(5,t)*z2  
     >    +au(6,t)*z3+au(7,t)*z4)                             
      B_x=bu(1,t)+bu(2,t)*x+bu(3,t)/(x+bu(4,t))                             
      C_x=cu(1,t)*x+cu(2,t)*x**2+cu(3,t)*x**3+cu(4,t)*x**4     
      f2u   =A_x*(log_term)**B_x*(1+C_x/qsq)

      A_x=x**al(1,t)*z**al(2,t)*(al(3,t)+al(4,t)*z+al(5,t)*z2  
     >    +al(6,t)*z3+al(7,t)*z4)                             
      B_x=bl(1,t)+bl(2,t)*x+bl(3,t)/(x+bl(4,t))                             
      C_x=cl(1,t)*x+cl(2,t)*x**2+cl(3,t)*x**3+cl(4,t)*x**4     
      f2l   =A_x*(log_term)**B_x*(1+C_x/qsq)

      err_hi = f2u-f2
      err_lo = f2-f2l

      return                                                            
      end                                                               


c========================================================================
C     ===

      SUBROUTINE i_d2_model(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper,XVAL)
********************************************************************************
*
* This subroutine calculates model cross sections for deuterium in resonance region. 
* Cross section is returned in nanobarns. 
* For all the resonance regions we use nonrelativistc
* Breit-Wigner shapes (L.Stuart). 
* This subroutine is based on h2model.f from SLAC
* 
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
********************************************************************************
      IMPLICIT NONE
*     
      REAL*4 XVAL(37)
*     
      logical goroper
      INTEGER I, J
*
* modified by I.N. 04/09/97
*

       
      REAL*4 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2)
      REAL*4 KRCM, EPIRCM, PPIRCM   
      REAL*4 W, KCM, K, EPICM, PPICM, WDIF 
 
      REAL*4 PI, ALPHA, AM, MPPI, MPI, XR  

      integer jj
      real*4 kr,meta
      real*4 mrho,mn,am2
      real*4 gamma_gamma,gamma_pi,denom_pi
      real*4 sigma_delta,sigma_second,sigma_third
*
      real*4 mass(3),width(3),qdep(3)
      real*4 coeff(3,0:3),nres_coeff(0:4,3),temp
* Define vectors : MASS(3), WIDTH(3), which will contain the masses and widths of
* the three resonances: 1=DELTA, 2=second res reg, 3=third res reg
*
      DATA PI/3.14159265/
      DATA ALPHA/0.00729735/
      DATA AM/.93828/
      DATA MPI/0.13957/
      DATA MPPI/1.07783/ 
      data Meta/0.54745/ 
      data Mrho/0.7681/
*     DATA XR/0.18/
      data Mn/0.938/
      Am2 = Am*Am
*     
*     The fit parameters are called  XVAL 
*     Assign them to masses, widths, etc.
*     
      xr=XVAL(1)
      do i=1,3
        mass(i)  = XVAL(i+1)
*     write(*,*)xval(i)
      enddo
      do i=1,3
        width(i)  = XVAL(i+4)
*     write(*,*)xval(i+3)
      enddo
      do i=1,3
        qdep(i)   = XVAL(i+7)
*     write(*,*)xval(i+6)
      enddo     
      k=11
      do i=1,3
        do j=0,3
          coeff(i,j) = XVAL(k)
*     write(*,*)xval(k)
          k=k+1
        enddo
      enddo
      do i=0,4
        do j=1,3
          nres_coeff(i,j) = XVAL(k)
***   !!!            write(*,*)'3',xval(k)
          k=k+1
        enddo
      enddo
*     write(*,*)'Done assigning parameters'
*     
*     write(*,*) 'w2=', w2
      W = SQRT(W2) 
      temp=w-mppi
      WDIF = MAX(real(0.0001), real(temp))
****************************************************************
*     This part is not used for deuterium.
****************************************************************
*     K equivalent real photon energy needed to excite resonance W
*     K is in Lab frame, KCM is in CM frame
*     
      K = (W2 - Am2)/2./Am
      KCM = (W2 - Am2)/2./W
      KR = (Mass(1)*Mass(1)-Am2)/2./Am
      KRCM = (Mass(1)*Mass(1) - Am2)/2./Mass(1)
*     Decay mode for Delta: N,pi:
********
      EPICM = 0.5*( W2 + Mpi*Mpi - Am2 )/W
      
      PPICM = SQRT(MAX(real(0.0),(EPICM*EPICM - MPI*MPI)))
      EPIRCM = 0.5*(Mass(1)*Mass(1) + Mpi*Mpi - Am2 )/Mass(1)   
      PPIRCM = SQRT(MAX(real(0.0),(EPIRCM*EPIRCM - MPI*MPI)))
      
********
*     Now calculate the partial widths:
*     gamma_pi(i)
*     write(*,*)'mass(1),w2 ',mass(1),' ' ,w2
*     write(*,*)'width(1) ',width(1)
*     write(*,*)'PPICM, PPIRCM ',PPICM,' ',PPIRCM
*     write(*,*)'XR ' ,xr
      gamma_pi = width(1)*(PPICM/PPIRCM)**3*(PPIRCM*PPIRCM+XR*XR)/
     >  (PPICM*PPICM+XR*XR)
*     
*     write(*,*)'iuhu 2'
      gamma_gamma = width(1)*(KCM/KRCM)**2*(KRCM*KRCM+XR*XR)/
     >  (KCM*KCM+XR*XR)
      
      denom_pi = (W2 - Mass(1)*Mass(1))**2 + (Mass(1)*gamma_pi)**2
*     write(*,*)'iuhu 3'
***********************************************************************
***********************************************************************
*     Now calculate cross sections for Delta:
      sigma_delta = width(1)/((W-Mass(1)*(1. + Q2*qdep(1)))**2 
     >  + 0.25*width(1)*width(1))
*     For the second and third resonance regions:
      sigma_second = width(2)/((W-Mass(2)*(1. + Q2*qdep(2)))**2 
     >  + 0.25*width(2)*width(2))
      
*     write(*,*)'width(2)',' ',width(2)
*     write(*,*)'mass(2)',' ',mass(2)
      
*     write(*,*)'sigma_second',' ',sigma_second
      
      sigma_third  = width(3)/((W-Mass(3)*(1. + Q2*qdep(3)))**2
     >  + 0.25*width(3)*width(3))
      
*     write(*,*)'sigma_third',' ',sigma_third
      
*     write(*,*)'width(3)',' ',width(3)
*     write(*,*)'mass(3)',' ',mass(3)
      
*     write(*,*)'iuhu 5'
*     
*     Put in the Q2 dependence of AH(Q2):
*     AH(Q2) = coeff(i,j)*Q2**j
*     
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(1) = 0.0
      do j=0,3  
        sig_res1(1) = sig_res1(1)+ sigma_delta*coeff(1,j)*Q2**j
        sig_res2(1) = sig_res2(1)+
     >    sigma_second*coeff(2,j)*Q2**j+
     >    sigma_third *coeff(3,j)*Q2**j
      enddo   
*     write(*,*)'sigma_d,sigma_23',' ',sig_res1(1),' ',sig_res2(1)
      do j=0,4
        do jj=1,3
          sig_nres(1) = sig_nres(1)+
     >      nres_coeff(j,jj)*q2**j*sqrt(wdif**(2*jj-1))
*     write(*,*)'2',sig_nres(1),nres_coeff(j,jj),j,wdif,jj
        enddo
      enddo     
*     
      RETURN
      END



      SUBROUTINE INEFT(QQ,W,W1,W2)
C Modified 6feb87 by lww to accept target information passed through    
C common block /targt/.                                                 
                                                                        
C This program takes the old slac structure function model (Atwood,     
C Bodek, et.al.) and outputs values for W1 and W2 at given kinematics. 
! As of 11/3/95 this version is per NEUCLEON   ! Steve Rock

! amuM is atomic number, ie. 1. 2.xxx etc.
! 11/05 changed to FITEMC_N (with neutron excess) instead of FITEMC
! 8/06 changed back to FITEMC and put in neutron excess from
! diff of d and p fits pyb peb
c this doesnt work: put in CONSTNAT n/p ratio for test
      Implicit None 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      integer iz,ia
      REAL avgN, avgA, avgM, amuM               
      REAL*8 QQ,W,W1,W2,WW,V,VV,OMEGAP,SP,UNIV,BRES,SLACF2,B,fneut
      REAL*8 VW2,X,EMCFAC,UNIVD,BRESD,UNIVP,BRESP,vW2p,vW2d,vW2n
      REAL*8    C(24),CF(11),CD(24),CFD(11)                          
      REAL*8    EF(7) 
      REAL FITEMC
c      REAL FITEMC_N
                                                  
      REAL*8         PM / .93828/,PMPM/.8803/,TPM/1.876512/            
      REAL*8         R /  .18/,ALPHAX/137.0388/,THCONST/0.0174533/
      LOGICAL GOODFIT
      common/testing/prttst
      logical prttst
      DATA         EF / -0.00136693,-.00510425,-.0375986,-.0946004,     
     +                  -.122435,-.0112751,0.406435/                    
                                                                        
C FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE TERMS)    
                                                                        
      DATA   C(1) / 0.10741163E 01/,  C(2) / 0.75531124E 00/,           
     *       C(3) / 0.33506491E 01/,  C(4) / 0.17447015E 01/,           
     *       C(5) / 0.35102405E 01/,  C(6) / 0.10400040E 01/,           
     *       C(7) / 0.12299128E 01/,  C(8) / 0.10625394E 00/,           
     *       C(9) / 0.48132786E 00/,  C(10)/ 0.15101467E 01/,           
     *       C(11)/ 0.81661975E-01/,  C(12)/ 0.65587179E 00/,           
     *       C(13)/ 0.17176216E 01/,  C(14)/ 0.12551987E 00/,           
     *       C(15)/ 0.74733793E 00/,  C(16)/ 0.19538129E 01/,           
     *       C(17)/ 0.19891522E 00/,  C(18)/-0.17498537E 00/,           
     *       C(19)/ 0.96701919E-02/,  C(20)/-0.35256748E-01/,           
     *       C(21)/ 0.35185207E 01/,  C(22)/-0.59993696E 00/,           
     *       C(23)/ 0.47615828E 01/,  C(24)/ 0.41167589E 00/            
                                                                        
C CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION) OMEGAW FIT   
                                                                        
      DATA CF(1) / 0.25615498E 00/,  CF(2) / 0.21784826E 01/,           
     *     CF(3) / 0.89783738E 00/,  CF(4) /-0.67162450E 01/,           
     *     CF(5) / 0.37557472E 01/,  CF(6) / 0.16421119E 01/,           
     *     CF(7) / 0.37635747E 00/,  CF(8) / 0.93825625E 00/,           
     *     CF(9) / 0.10000000E 01/,  CF(10)/ 0.0           /,           
     *     CF(11)/ 0.50000000E 00/                                      
                                                                        
C FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONANT TERMS)   
                                                                        
      DATA  CD(1) / 0.10521935E 01/, CD(2) / 0.76111537E 00/,           
     *      CD(3) / 0.41469897E 01/, CD(4) / 0.14218146E 01/,           
     *      CD(5) / 0.37119053E 01/, CD(6) / 0.74847487E 00/,           
     *      CD(7) / 0.12399742E 01/, CD(8) / 0.12114898E 00/,           
     *      CD(9) / 0.11497852E-01/, CD(10)/ 0.14772317E 01/,           
     *      CD(11)/ 0.69579815E-02/, CD(12)/ 0.12662466E 00/,           
     *      CD(13)/ 0.15233427E 01/, CD(14)/ 0.84094736E-01/,           
     *      CD(15)/ 0.74733793E 00/, CD(16)/ 0.19538129E 01/,           
     *      CD(17)/ 0.19891522E 00/, CD(18)/-0.24480414E 00/,           
     *      CD(19)/ 0.14502846E-01/, CD(20)/-0.35256748E-01/,           
     *      CD(21)/ 0.35185207E 01/, CD(22)/-0.21261862E 00/,           
     *      CD(23)/ 0.69690531E 01/, CD(24)/ 0.40314293E 00/            
                                                                        
C CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)        
C OMEGAW FIT                                                            
                                                                        
      DATA CFD(1) / 0.47708776E 00/, CFD(2) / 0.21601918E 01/,          
     *     CFD(3) / 0.36273894E 01/, CFD(4) /-0.10470367E 02/,          
     *     CFD(5) / 0.49271691E 01/, CFD(6) / 0.15120763E 01/,          
     *     CFD(7) / 0.35114723E 00/, CFD(8) / 0.93825625E 00/,          
     *     CFD(9) / 0.10000000E 01/, CFD(10)/ 0.0           /,          
     *     CFD(11)/ 0.50000000E 00/                                     
                                                                        
C COMPUTE SOME KINEMATIC QUANTITIES                                     
                                                                        
      WW     = W**2                                                     
      V      = (WW+QQ-PMPM)/2.D0/PM                                     
      VV     = V*V                                                      
      OMEGAP = TPM*V/QQ+PMPM/QQ                                         
                                                                        
C OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION                      
      OMEGAP = DMIN1(20.0D0,OMEGAP)                                       
                                                                        
      SP = 1.0-EXP(-7.7*(OMEGAP-1.0))                                   
C pyb *** modified to get proton too when IA ne IZ*2
      IF (IA. ne. 2*IZ) THEN 
C          UNIVERSAL AND RESONANCE FIT FOR HYDROGEN                     
           UNIVP = SLACF2(W,QQ,CF)                                      
           BRESP = B(W,QQ,C)                                            
      endif
      if(amuM.ge.1.5) then
C          UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM                    
           UNIVd = SLACF2(W,QQ,CFD)/SP
           BRESd = B(W,QQ,CD)          
      ENDIF                                                             
                                                                        
C COMPUTE VW2,W2,W1                                                     
                                                                        
      if(amuM.le.1.5) then
        vW2 = UNIVp * BRESp
      else
        VW2p    = UNIVp * BRESp 
        VW2d    = UNIVd * BRESd / 2. ! per nucleon
        vW2n    = 2.*vW2d - vW2p
        fneut = float(IA - 2*IZ) / float(IA)
c       if(fneut.lt.0.) write(6,'(''ERROR, fneut='',f8.3,2i3)') 
c    >    fneut,ia,iz
c *** set fneut back to zero: assume n/p=1!
        fneut=0.
        vW2 = (1. - fneut) * vW2d + fneut * vW2n
        if(prttst) write(98,'(8f8.3)')  qq,w,vw2p,vw2d,vw2n,
     >    vw2,vw2n/vw2d,fneut
      endif

      W2     = VW2/V                                                    
      W1     = (1.0D0+VV/QQ)/(V*(1.0D0+R))*VW2                          
!      if(prttst) write(*,'(1x,''univ...='',6f10.4)') sp,univ,bres,
!     >  vw2,w2,w1

      IF (amuM.LE.2.5) RETURN                                               
      X      = QQ/2./PM/V
c Modified 11/05 to include neutron excess pyb peb
c Modified 8/06 to do neutron excess as above, not in EMC
      EMCFAC= FITEMC(REAL(X),REAL(amuM),GOODFIT)
c      EMCFAC= FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GOODFIT)
                                                                        
      W2     = W2*EMCFAC                                                
      W1     = W1*EMCFAC                                                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C-----------------------------------------------------------------------
!---------------------------------------------------------------------
      REAL FUNCTION FITEMC_N(X,A,Z,GOODFIT)
     +  
!---------------------------------------------------------------------  
! Modified FITEMC.F with Neutron excess correction and proton=1 added 8/19/98
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! Z= number of protons
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 8/19/98 **  If proton, return 1.
! 8/19/98 **  Add neutron excess correction.
! 11/05 Modified PYB to return value at x=0.7 if x>0.7 since
!       Fermi smearing in INELASTU will account for Fermi smearing rise
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL ALPHA, C,LN_C,X,A,Z ,X_U,SIG_N_P,F_IS
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points
!Term    Coeficient     Error
      REAL*8 ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                 
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /                                  

             !
!Chisq=         22.    for 30 points
!Term    Coeficient     Error 
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        

      IF(A.LT.1.5) THEN    ! Added 8/19/98
       FITEMC_N=1.
       GOODFIT=.TRUE.
       RETURN
      ENDIF                                                                                
c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c      IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC_N  =  C *A**ALPHA    !isoscaler
      SIG_N_P = 1.-0.8*X_U
      F_IS = .5*(1.+ SIG_N_P)/(Z/A +(1.-Z/A)*SIG_N_P)
      FITEMC_N = FITEMC_N/F_IS
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------
      REAL FUNCTION FNP_NMC(X,QSQ)
!---------------------------------------------------------
! NMC FIT TO NMC,SLAC,NMC data in CERN-PPE/91-167
! No Fermi Motion Corrections
!  Steve Rock 1 Feb 1996
!-------------------------------------------------------------
      IMPLICIT NONE
      REAL X,QSQ,A,B,X2,X3
    
      X2 = X*X
      X3 = X2*X
      A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
C$$      B = -.171*X2 + .277*X3
      B = -.171*X2 + .244*X3  ! replaced 10/22/97 by correct value on x3
      FNP_NMC = A *(1+X2/QSQ)*(QSQ/20.)**B
      RETURN
      END

!---------------------------------------------------------------------
      SUBROUTINE R1990F(X,QQ2,R,ERR)
C
C    Ref: L.W.Whitlow SLAC Report 357 (1990)         and
C         L.W.Whitlow et al.: PL B250 (1990) 193
C
C    for Q2 < 0.35 we extrapolate R as a constant with a rough error of 100 %
C
      REAL A(3)  / .06723, .46714, 1.89794 /,
     >     B(3)  / .06347, .57468, -.35342 /,
     >     C(3)  / .05992, .50885, 2.10807 /
C
      DATA QMAX /64./
      Q2=QQ2
      IF(Q2.LT.0.35) Q2=0.35
      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))
      RLOG  = FAC/LOG(Q2/.04)
      Q2THR = 5.*(1.-X)**5
      RA   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))
      RB   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)
      RC   = C(1)*RLOG + C(2)/SQRT((Q2-Q2THR)**2+C(3)**2)
      R     = (RA+RB+RC)/3.
      IF (Q2.GE.0.35) THEN
        Q = MIN(Q2,QMAX)
        S = .006+.03*X**2
        AA= MAX(.05,8.33*X-.66)
        XLOW  = .020+ABS(S*LOG(Q/AA))
        XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)
        D1SQUARE=XLOW**2+XHIGH**2
        D2SQUARE= ((RA-R)**2+(RB-R)**2+(RC-R)**2)/2.
        D3    = .023*(1.+.5*R)
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3
        ERR    = SQRT(D1SQUARE+D2SQUARE+D3**2)
      ELSE
        ERR = R
      ENDIF
      RETURN
      END

!---------------------------------------------------------------------

      FUNCTION f2allm(x,q2)             
!-----------------------------------------------------------------
c       allm97, NMC published measured points Q2>0.75 GeV2
c       for values Q<1 use data of E665!
c       parameterization of F2 , according to
c       H.Abramowicz and A.Levy, hep-ph/9712415
c
c       3*10-6 < x  < 0.85, W2>3GeV2
c       0.   < Q2 < 5000 GeV2, dof=0.97
c
! From Peter Bosted, Nov 00 from Antje Burrel.
!------------------------------------------------------------------
      IMPLICIT NONE
      REAL F2ALLM,X,Q2
      REAL SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      COMMON/ALLM/SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
C  POMERON
      REAL S11,S12,S13,A11,A12,A13,B11,B12,B13,M12
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      REAL S21,S22,S23,A21,A22,A23,B21,B22,B23,M22
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      REAL M02,LAM2,Q02,ALFA,XMP2
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      REAL W2,W,Z
C                                                                               
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
C
      IF(Q2.EQ.0.) THEN                                                        
       S=0.
       Z=1.           
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                                   
       AP=A11                                                            
       BP=B11                                                               
       SP=S11                                                             
       F2P=SP*XP**AP                                                
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))          
       AR=A21                                 
       BR=B21                                 
       SR=S21
       F2R=SR*XR**AR              
C                                                                               
      ELSE                                                                      
       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))                       
       Z=1.-X                                      
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)         
       BP=B11+B12*S**B13                                                  
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)         
       F2P=SP*XP**AP*Z**BP                            
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))                                          
       AR=A21+A22*S**A23                                                     
       BR=B21+B22*S**B23                                                
       SR=S21+S22*S**S23                                                     
       F2R=SR*XR**AR*Z**BR                                                   
    
C                                                                               
      ENDIF

c      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z              
c      SIGal=CIN*(F2P+F2R)                                             
c      f2allm=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)
      F2ALLM = q2/(q2+m02)*(F2P+F2R)
 

      RETURN                                                                    
      END                                                                       

!---------------------------------------------------------------------

      REAL FUNCTION FITEMC(X,A,GOODFIT)                                         
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
!    beyond that assume Fermi motion is giving the rise, and we
!    already are taking that into account with the y-smearing of
!    the inelastic
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*4 ALPHA, C,LN_C,X,A ,X_U
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c       IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------

      REAL FUNCTION  ROCK_RES9(TARG,X,Q2)
! April 22 2003 Version of deuterium fit.
! No JLAB data

      IMPLICIT NONE
      INTEGER TARG ! 1=H, 2=D
      REAL Q2
      REAL X,F2INEL
      REAL*8 WM,B03
      CHARACTER*1 TARG_STR ! D or H      
      REAL MP2/.8803/,MP/.93828/ 
! DATE =2003/ 4/22  TIME=13: 3
! MODEL# 9
! HYDROGEN
      REAL*8 CH(40)/
     >  1.074100000, 0.279494473, 4.712278073, 1.739263857, 1.431155682,
     >  0.184080601, 1.225717385, 0.110613494, 0.214285318, 1.510392164,
     >  0.105595846, 0.309575075, 1.725500000, 0.138478762,-1.452323856,
     >  1.953800000, 0.342265136,-0.044517487, 0.093869839,-0.029404837,
     >  1.199536964, 0.671722156, 2.033898985, 0.024107348, 0.452992650,
     >  0.624721200, 0.054377925, 0.099384033, 0.897445331,13.205228916,
     >  6.284002261, 0.013794386, 0.030764850,-0.007500832, 0.025566716,
     >  0.016759263, 0.966305525, 1.058950737, 0.000000000, 0.000000000/
     > 
!deuterium
      REAL*8 CD(40)/
     >   1.05220000,  0.75530000,  3.30000000,  1.70000000,  3.50000000,
     >   1.04000000,  1.22000000,  0.10000000,  0.48000000,  1.50000000,
     >   0.08000000,  0.65000000,  1.70000000,  0.12000000,  0.74000000,
     >   1.90000000,  0.19000000, -0.17000000,  0.00960000,  0.03500000,
     >   3.50000000, -0.60000000,  4.70000000,  0.41100000,  0.41100000,
     >   0.10000000,  0.41100000,  0.10000000,  0.10000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     >   0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000/
    
      IF(TARG.EQ.1) TARG_STR ='H'
      IF(TARG.EQ.2) TARG_STR ='D'
      WM = SQRT(MP2 +Q2*(1./X -1.) )

      CALL F2GLOB_FAST(X,Q2,TARG_STR,9,F2INEL)

      IF(TARG.EQ.1) THEN  ! H2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CH)   
      ELSE     ! D2
        ROCK_RES9 = F2INEL*B03(WM,DBLE(Q2),0,CD)   
      ENDIF
      RETURN
      END
!----------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                                        
      REAL FUNCTION DR1998(X,Q2)                                             
                                                                        
! Parameterizes the uncertainty in R1990 due to the statistical         
! fluctuations in the data.  Values reflect an average of the R-values  
! about a neighborhood of the specific (x,Q2) value.  That neighborhood 
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see       
! Reference.                                                            
!                                                                       
! This subroutine is accurate over all (x,Q2), not only the SLAC deep   
! inelastic range.  Where there is no data, for example in the resonance
! region, it returns a realistic uncertainty, extrapolated from the deep
! inelastic region (suitably enlarged).  We similarly estimate the      
! uncertainty at very large Q2 by extrapolating from the highest Q2     
! measurments.  For extremely large Q2, R is expected to fall to zero,  
! so the uncertainty in R should not continue to grow.  For this reason 
! DR1990 uses the value at 64 GeV for all larger Q2.                    
!                                                                       
! XHIGH accounts for the rapidly diminishing statistical accuracy for   
! x>.8, and does not contribute for smaller x.                          
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      REAL Q2,X                   

      DR1998 = .0078 -.013*X +(.070 -.39*X+.70*X**2)/(1.7+Q2)
      RETURN                                                            
      END                                                               

!------------------------------------------------------------------------

                                                                        
C=======================================================================
                                                                        
      REAL Function Fgauss(T)                                                
      IMPLICIT NONE                                                                             
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL T, RADIUS,X2,CHAR
                                                                        
      Fgauss = 0.                                                       
      Radius = 1.07*avgA**(1./3.)  
! from H. de Vries: Nuclear Charge Density Distributions from Elastic Electron
!   Scattering. in Atomic Data and Nuclear Data Tables 36, 495(1987)
! 8/9/96
      IF(IA.EQ.205)RADIUS=5.470
      IF(IA.EQ.56) RADIUS=3.729    
      If(IA.EQ.28) RADIUS=3.085
      IF(IA.EQ.27) RADIUS=3.035                                                        
      x2     = (T/.197328**2)*Radius**2                                 
      char   = (T/.197328**2)*(2.4**2)/6.                               
      if (char.lt.80) Fgauss = exp(-char)/(1.+x2/6.)                    
      Return                                                            
      End                                                               
                                                                        
C=======================================================================
                                                                        
      REAL Function Fshell(T)                                                
      IMPLICIT NONE                                                                          
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM   
      REAL T, RADIUS,X2,ALP,CHAR
                                                                     
      Fshell = 0.                                                       
      Radius = 1.07*avgA**(1./3.)     
! from H. de Vries: Nuclear Charge Density Distributions from Elastic Electron
!   Scattering. in Atomic Data and Nuclear Data Tables 36, 495(1987)
!8/9/96
      IF(IA.EQ.16) RADIUS=2.737
      IF(IA.EQ.15) RADIUS=2.611
      IF(IA.EQ.14) RADIUS=2.540
      IF(IA.EQ.12) RADIUS=2.464  
      IF(IA.EQ. 9) RADIUS=2.519
      IF(IA.EQ. 7) RADIUS=2.39 
      IF(IA.EQ. 6) RADIUS=2.56 
      IF(IA.EQ. 4) RADIUS=1.676                                       
      x2     = (T/.197328**2)*Radius**2                                 
      alp    = (iZ-2.)/3.                                               
      char   = x2*(2.+3.*alp)/(12.+30.*alp)                             
      if (char.lt.80) Fshell = exp(-char)*(1.-alp*x2/(6.+15.*alp))      
      Return                                                            
      End                                                               
                                                                        
C=======================================================================
                                                                        
      REAL FUNCTION FDEL(T)   
      IMPLICIT NONE       
      REAL SQF,T
                                        
      SQF  = SQRT(T/.197328**2)                                         
      FDEL = 1.58/SQF*(ATAN(SQF/0.93)-2.*ATAN(SQF/3.19)+ATAN(SQF/5.45)) 
      IF (FDEL.LE.0.) FDEL = 0.                                         
      RETURN                                                            
      END                                                               
                                                                        
!----------------------------------------------------------------
      REAL function FF_BESSEL ( T ,OUT_OF_RANGE)
      IMPLICIT NONE
C     Calculates PWBA form factor  for Si-28 or O-16 using 
C     Fourier-Bessel coefficients given by H. de Vries et al., 
C     Atomic Data and Nuclear Data Tables (1986).
C     
C     Note: Other nuclides can be entered, but at this time this
C     is only good for 
C     He-3, C-12, N-15, O-16, Al-27, Si-28, Fe-56, Cu-65,
!----------------------------------------------------------------
!     3/28/03 Corrected for divide 0/0 when qm =0  - Steve Rock
      
      REAL T
      LOGICAL OUT_OF_RANGE
      common /targt/ iZ, iA, avgN, aavgA, avgM, amuM
      real avgN, aavgA, avgM, amuM
      integer iZ, iA
      real*8 a(20),a16(14),a28(14),a3(12),a12(16),a15(8),a27(12),
     >  a56(17), a65(17),a196(15)
      INTEGER N,J,I
      REAL*8 R_MAX,QMU,Q,FF,QMUX,QM,QP,SINQR
      
      REAL*8 PI4/12.56637062/   ! 4*PI
      REAL*8 PI2/ 6.28318531/   ! 2*PI
C     
C     de Vries, de Vries & de Jager FB coefficients:


!3He
      data a3/0.20020E-01, 0.41934E-01, 0.36254E-01, 0.17941E-01,
     1        0.46608E-02, 0.46834E-02, 0.52042E-02, 0.38280E-02,
     2        0.25661E-02, 0.14182E-02, 0.61390E-03, 0.22929E-03/
!12C
      data a12/0.15721E-01, 0.38732E-01, 0.36808E-01, 0.14671E-01,
     1        -0.43277E-02,-0.97752E-02,-0.68908E-02,-0.27631E-02,
     2        -0.63568E-03, 0.71809E-04, 0.18441E-03, 0.75066E-04,
     3         0.51069E-04, 0.14308E-04, 0.23170E-05, 0.68465E-06/


!15N7
      data a15/0.25491E-01, 0.50618E-01, 0.29822E-01, -0.55196E-02,
     1        -0.15913E-01,-0.76184E-02, -.23992E-02, -0.47940E-03/
!16 O8
      data a16/0.20238E-01, 0.44793E-01, 0.33533E-01, 0.35030E-02,
     1          -0.12293E-01,-0.10329E-01,-0.34036E-02,-0.41627E-03,
     2          -0.94435E-03,-0.25571E-03, 0.23759E-03,-0.10603E-03,
     3           0.41480E-04, 0.00000E-03/

!27Al
      data a27/0.43418E-01, 0.60298E-01,  0.28950E-02, -0.23522E-01,
     1        -0.79791E-02, 0.23010E-02,  0.10794E-02,  0.12574E-03,
     2        -0.13021E-03, 0.56563E-04, -0.18011E-04,  0.42869E-05/

!28Si
      data a28/0.33495E-01, 0.59533E-01, 0.20979E-01,-0.16900E-01,
     1          -0.14998E-01,-0.93248E-03, 0.33266E-02, 0.59244E-03,
     2          -0.40013E-03, 0.12242E-03,-0.12994E-04,-0.92784E-05,
     3           0.72595E-05,-0.42096E-05/

!56Fe26 
      data a56/ 
     1  .42018E-1,  .62337E-1,  .23995e-3, -.32776E-1, -.79941E-2,
     2  .10844E-1,  .49123e-2, -.22144e-2, -.18146E-3,  .37261E-3,
     3 -.23296E-3,  .11494E-3, -.50596E-4,  .20652E-4, -.79428E-5,
     4  .28986E-5, -.10075E-5/         

!65Cu29 
        data a65/0.45444E-01, 0.59544E-01, -.94968E-02, -.31561e-01,
     1           0.22898E-03, 0.11189E-01, 0.37360E-02, -.64873E-03,
     2	         -.51133E-03, 0.43765E-03, -.24276E-03, 0.11507E-03, 
     3           -.49761E-04, 0.20140E-04, -.76945E-05, 0.28055E-05,
     4		 -.97411E-06 /


!196Pt78
         data a196/
     1     .50218e-1,  .53722e-1, -.35015e-1, -.34588e-1,  .23564e-1,
     2     .14340e-1, -.13270e-1, -.51212e-2,  .56088e-2,  .14890e-2,
     3    -.10928e-2,  .55662e-3, -.50557e-4, -.19708e-3,  .24016e-3/  

C	Change to units of fm**(-1)
	q = sqrt ( T ) / 0.197328
	if( q .le. 0.005 ) then
            OUT_OF_RANGE=.FALSE. ! line added 4/29/98
	    FF_BESSEL = 1.00000
	    return
	endif

        FF_BESSEL=0.
        OUT_OF_RANGE=.FALSE.
        R_max=0.
	if( iA .eq. 28 ) then
            if(q.gt.2.64) OUT_OF_RANGE=.TRUE.
	    R_max = 8.
	    n = 14
	    do i = 1,n
	        a(i) = a28(i)
	    enddo
        elseif( iA .eq. 16 ) then
            if(q.gt.2.77) OUT_OF_RANGE=.TRUE.
            R_max=8.
	    n = 13
	    do  i = 1,n
	      a(i) = a16(i)
	    enddo
        elseif( iA .eq. 3 ) then
            if(q.gt.10. ) OUT_OF_RANGE=.TRUE.
            R_max=5.
	    n = 12
	    do  i = 1,n
	      a(i) = a3(i)
	    enddo
        elseif( iA .eq. 12 ) then
            if(q.gt.4.01 ) OUT_OF_RANGE=.TRUE.
            R_max=8.
	    n = 16
	    do  i = 1,n
	      a(i) = a12(i)
	    enddo
        elseif( iA .eq. 15 ) then
            if(q.gt.3.17) OUT_OF_RANGE=.TRUE.
            R_max=7
	    n = 8
	    do  i = 1,n
	      a(i) = a15(i)
	    enddo
        elseif( iA .eq. 27 ) then
            if(q.gt.2.70) OUT_OF_RANGE=.TRUE.
            R_max=7
	    n = 12
	    do  i = 1,n
	      a(i) = a27(i)
	    enddo
        elseif( iA .eq. 56 ) then
            if(q.gt.2.22) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.51) OUT_OF_RANGE=.TRUE.
            R_max=9
	    n = 17
	    do  i = 1,n
	      a(i) = a56(i)
	    enddo

        elseif(( iA .eq. 64).or.(iA.eq.65) ) then
            if(q.gt.2.22) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.51) OUT_OF_RANGE=.TRUE.
            R_max=9
	    n = 17
	    do  i = 1,n
	      a(i) = a65(i)
	    enddo
        elseif( iA .eq. 196)  then
            if(q.gt.2.28) OUT_OF_RANGE=.TRUE.
            if(q.lt.0.34) OUT_OF_RANGE=.TRUE.
            R_max=12
	    n = 15
	    do  i = 1,n
	      a(i) = a196(i)
	    enddo
	else      
            out_of_range=.true.
	endif
        if(out_of_range.or.r_max.eq.0.) then
          ff_bessel=0.
          Return
        endif 


	qmu = 3.14159265 / R_max

        ff=0.
        sinqR = sin(q*R_max)
        do j=1,n
	 qmux = qmu * float(j)
	 qm =  q - qmux
         qp =  q + qmux
         if(abs(qm).gt.1.e-6) then
          ff=ff+ a(j)*((-1.)**j)*sinqR/(qm*qp)
         else
          ff= ff +a(j)*R_max**2/(PI2*j)
         endif
        enddo
        if((q*R_max.gt.1.E-20).and.(ff.lt.1.E20)) then
         ff_bessel =PI4/FLOAT(IZ)/q *ff
        else
         ff_bessel=0.
        endif 
	Return
	End

                                                                       

      REAL FUNCTION F2LOWW(ITARGE,Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*   A V Kotwal 30dec93                                       *
*                                                            *
*Return F2 for 1.1<W<5 over a large range of Q2 down to Q2=0 *
*                                                            *
* Based on parametrizations of SLAC,DESY and Daresbury data. *
* Calls two other routines for the resonance region 1.1<W<2  *
* and the inelastic region 2<W<5                             *
* Assume proton and deuteron same for resonance region       *
*                                                            *
* ITARG = 1 for H2, 2 for D2                                 *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
*                                                            *
*     PGS Feb 20, 94                                         *
* Proton and neutron are not the same, see Close chapter 7.2 *
* Modify deuteron resonances accordingly.                    *
* Simple scaling of the D resonance parametrization          *
* according to data (flat line fit in 8 points) by 0.88      *
* Data from J. Franz et al Z.Phys.C 96 P. and F. 10,         *
* 105-116(1981)                                              *
* fit a flat line in the H/D plots p. 110                    *
**************************************************************
      REAL M2,TWOM,F2,Q2,NU,W2,W,WMIN,WINTER,WMAX,NUMAX,WFUNC
      PARAMETER (M2=0.8803,TWOM=1.87654) ! protonmass**2, 2*protonmass
      PARAMETER (WMIN=1.12,WINTER=1.95,WMAX=5.0)
      INTEGER ITARGE,ITARG
      REAL F2LONU,F2PRES,LONUF2,PRESF2,WDIFF
      EXTERNAL F2LONU,F2PRES
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      ITARG = ITARGE

      W2 = M2 + TWOM*NU - Q2
      IF (W2.LT.M2)  W2=M2
      W  = SQRT( W2 )

      IF (NU.LT.0.001) NU = 0.001
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (ITARG.GT.2) ITARG = 2
      IF (ITARG.LT.1) ITARG = 1

      IF (W.LT.WMIN) THEN
c.......force to zero smoothly as W->M
        PRESF2 = F2PRES(Q2,WMIN)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        F2 = PRESF2 * ((W2-M2)/(WMIN**2-M2))**2
      ELSEIF (W.LE.WMAX) THEN
c.....inelastic region, joined smoothly to resonance region
        WDIFF = (W-WINTER)/0.03
        IF (WDIFF.GT.10.0) THEN
          WFUNC = 0.0
        ELSEIF (WDIFF.LT.-10.0) THEN
          WFUNC = 1.0
        ELSE
          WFUNC = 1.0/(1.0+EXP(WDIFF))
        ENDIF
        PRESF2 = F2PRES(Q2,W)
        IF(ITARG.EQ.2)PRESF2 = PRESF2*0.88
        LONUF2 = F2LONU(ITARG,Q2,NU)
        F2=WFUNC*PRESF2+(1.0-WFUNC)*LONUF2
      ELSE
C.......force to zero at high W
        NUMAX = (WMAX**2-M2+Q2)/TWOM
        LONUF2 = F2LONU(ITARG,Q2,NUMAX)
c        F2 = LONUF2 * EXP(10.0*(WMAX-W))
        F2 = LONUF2
      ENDIF

      IF (F2.LT.0.0) F2 = 0.0

      F2LOWW = F2

      RETURN
      END

      REAL FUNCTION F2DHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(deuteron) for W>5                                *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2DDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2DI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2DI15,F2DDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2DI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2DDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2DHIW = F2

      RETURN
      END
      REAL FUNCTION F2PHIW(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) for W>5                                  *
* obtained by combining the model of Donnachie and Landshoff *
* valid at small Q2 with the NMC parametrization of their    *
* data at higher Q2.                                         *
* The two functions are merged at an intermediate Q2 where   *
* both functions are fits to NMC data                        *
*                                                            *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,TWOM,F2NMC,F2DL,F2,Q2MIN,Q2MAX,F2PDL,Q2MERG
      REAL Q2EXT,NUEXT,Q2NMC,Q2DOLA,XNMC,Q2FUNC
      DOUBLE PRECISION DX,DQ2,DF2NMC,F2HI15
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      PARAMETER (Q2MIN=0.1,Q2MAX=10.0)
      PARAMETER (Q2MERG=3.0)
      EXTERNAL F2HI15,F2PDL

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

c.....do not evaluate NMC parametrization below Q2MIN
      Q2NMC=Q2
      IF (Q2NMC.LT.Q2MIN) Q2NMC=Q2MIN
      XNMC = Q2NMC/(TWOM*NU)
      IF (XNMC.GT.0.99) XNMC=0.99
      IF (XNMC.LT.0.0001) XNMC=0.0001

      DX = DBLE(XNMC)
      DQ2 = DBLE(Q2NMC)
      DF2NMC = F2HI15(DX,DQ2)
      F2NMC = SNGL(DF2NMC)

c.....do not evaluate DOLA above Q2MAX
      Q2DOLA=Q2
      IF (Q2DOLA.GT.Q2MAX) Q2DOLA = Q2MAX
      F2DL = F2PDL(Q2DOLA,NU)

c.....now merge DOLA with NMC
      IF (Q2.GT.(Q2MERG+3.0)) THEN
        Q2FUNC = 0.0
      ELSE
        Q2FUNC = 1.0/(1.0+EXP((Q2-Q2MERG)/0.2))
      ENDIF
      F2 = F2DL*Q2FUNC + F2NMC*(1.0-Q2FUNC)

      IF (F2.LT.0.0) F2 = 0.0
      F2PHIW = F2

      RETURN
      END


       REAL FUNCTION F2PRES(Q2EXT,WEXT)
       IMPLICIT NONE
***********************************************************
* returns F2(proton) in the resonance region 1.1<W<2.0    *
* paremetrization of the cross-section is obtained from   *
* F.W.Brasse et al,                                       *
* ' Parametrization of the Q2 Dependence of the           *
* Gamma_v-P Total Cross-Sections in the Resonance Region' *
* Nuclear Physics B110, 413 (1976)                        *
*                                                         *
* for high energy muon scattering, the virtual photon     *
* polarization is close to 1 for low W, low Q2 scattering *
* such as resonance excitation.                           *
* Hence the parametrization for epsilon>0.9 is used from  *
* this reference. The average epsilon for this data is    *
* 0.957 . This is used together with R=sigmaL/sigmaT to   *
* convert the cross-section to F2                         *
*                                                         *
* A V Kotwal 28dec93                                      *
* inputs: Q2EXT = Q2 in GeV^2                             *
*         WEXT  = W  in GeV                               *
***********************************************************
       REAL M2,TWOM,F2P,Q2,NU,W2,W,Q,Q0,WLO,WHI,LNQOQ0,M,XX
       INTEGER NWBIN,IWBIN,N
       PARAMETER (NWBIN=56)
       REAL WMIN,WINTER,WMAX,DW1,DW2
       PARAMETER (WMIN=1.11,WINTER=1.77,WMAX=1.99)
       PARAMETER (DW1=0.015,DW2=0.02)
       PARAMETER (M=0.93828)    ! proton mass
       PARAMETER (M2=M*M,TWOM=2.0*M)
       REAL EPSILON,R,SIGMA,PI2AL4,Y,Y1,Y2,GD2,STRUW2,ERR
       PARAMETER (EPSILON=0.957)
       PARAMETER (PI2AL4=112.175)  ! 4.pi**2.alpha_em (microbarn.GeV**2)
       REAL Q2EXT,WEXT
       REAL A(NWBIN),B(NWBIN),C(NWBIN)
       DATA A /
     &  5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178
     & ,6.125,5.999,5.769,5.622,5.431,5.288,5.175,5.131
     & ,5.003,5.065,5.045,5.078,5.145,5.156,5.234,5.298
     & ,5.371,5.457,5.543,5.519,5.465,5.384,5.341,5.328
     & ,5.275,5.296,5.330,5.375,5.428,5.478,5.443,5.390
     & ,5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158
     & ,5.159,5.178,5.182,5.195,5.160,5.195,5.163,5.172 /
       DATA B /
     &  0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878
     & ,1.887,1.927,2.041,2.089,2.148,2.205,2.344,2.324
     & ,2.535,2.464,2.564,2.610,2.609,2.678,2.771,2.890
     & ,2.982,3.157,3.188,3.315,3.375,3.450,3.477,3.471
     & ,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519
     & ,4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322
     & ,5.546,5.623,5.775,5.894,6.138,6.151,6.301,6.542 /
       DATA C /
     &   0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080
     & ,-0.065,-0.056,-0.065,-0.056,-0.043,-0.034,-0.054,-0.018
     & ,-0.046,-0.015,-0.029,-0.048,-0.032,-0.046,-0.084,-0.115
     & ,-0.105,-0.159,-0.164,-0.181,-0.203,-0.220,-0.245,-0.264
     & ,-0.239,-0.302,-0.299,-0.318,-0.388,-0.393,-0.466,-0.588
     & ,-0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798
     & ,-1.048,-0.980,-1.021,-1.092,-1.313,-1.341,-1.266,-1.473 /

      Q2 = Q2EXT
      W  = WEXT
      IF (Q2.LT.0.0) Q2 = 0.0
      IF (W.LT.M) W = M
      W2 = W*W
      NU = (W2+Q2-M2)/TWOM

      IF (NU.LT.0.001) NU = 0.001

* Q = three momentum transfer to the hadronic system in the lab frame
      Q = ABS(SQRT( Q2 + NU*NU ))
      IF (Q.LT.0.01) Q = 0.01

* Q0 is the value of Q at q2=0 for the same W
      Q0 = (W2 - M2)/TWOM
      IF (Q0.LT.0.01) Q0 = 0.01

c.....find the bin of W, and bin edges
      IF (W.GE.WMAX) THEN
        IWBIN = NWBIN
        WLO   = WMAX
        WHI   = WMAX
      ELSEIF (W.GT.WINTER) THEN
        N     = INT((W-WINTER)/DW2)
        IWBIN = N + 45
        WLO   = FLOAT(N)*DW2 + WINTER
        WHI   = WLO + DW2
      ELSEIF (W.GT.WMIN) THEN
        N     = INT((W-WMIN)/DW1)
        IWBIN = N + 1
        WLO   = FLOAT(N)*DW1 + WMIN
        WHI   = WLO + DW1
      ELSE
        IWBIN = 1
        WLO   = WMIN
        WHI   = WMIN
      ENDIF

      LNQOQ0 = ALOG(Q/Q0)

c.....eqn 3
c.....parameter d is fixed d=3 in section 3.1
      Y1 = A(IWBIN) + B(IWBIN)*LNQOQ0 + C(IWBIN)*(ABS(LNQOQ0))**3
      IF (Y1.LT.0.0) Y1 = 0.0

      IF (IWBIN.LT.NWBIN) THEN
       Y2=A(IWBIN+1)+B(IWBIN+1)*LNQOQ0+C(IWBIN+1)*(ABS(LNQOQ0))**3
c.....linear interpolation inside bin
       Y = Y1 + (W-WLO)*(Y2-Y1)/(WHI-WLO)
      ELSE
       Y = Y1
      ENDIF
      IF (Y.LT.0.0) Y = 0.0

* The parametrization returns log(sigma/GD**2) where GD is the
* nucleon dipole form factor, and sigma is the total virtual photoabsorption
* cross-section
      GD2   = 1.0 / (1.0 + Q2/0.71)**4
      SIGMA = GD2*EXP(Y)

c.....R=sigmaL/sigmaT
c.....use 1990 SLAC analysis result
      XX = Q2/(TWOM*NU)
      IF (XX.GT.0.99) XX=0.99
      IF (XX.LT.0.0) XX=0.0
      CALL R1990F(XX,Q2,R,ERR)

c.....'An introduction to Quarks and Leptons', F.E.Close
c......eqn 9.44,9.46
c......also eqn 1 from reference
      STRUW2 = SIGMA*(1.0+R)*Q2 / ((1.0+EPSILON*R)*(Q2+NU*NU))
c......the Hand convention for the virtual photon flux has been followed
c......(reference 3 in this reference)
      STRUW2 = STRUW2*Q0/PI2AL4
      F2P    = NU*STRUW2
      IF (F2P.LT.0.0) F2P = 0.0

9999  CONTINUE
      F2PRES = F2P

      RETURN
      END


      DOUBLE PRECISION FUNCTION F2HI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** hydrogen BCDMS-like function
      DATA DPEMC /-0.1011D0,2.562D0,0.4121D0,-0.518D0,5.957D0,
     &           -10.197D0,4.695D0,
     &             0.364D0,-2.764D0,0.015D0,0.0186D0,-1.179D0,
     &             8.24D0,-36.36D0,47.76D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2HI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)
C
      RETURN
      END
 
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION SLACF2(WM,QSQ,CF)                                        
                                                                        
C UNIVERSAL FUNCTION FOR ATWOOD'S FIT                                   

      Implicit none                                      
      REAL*8    WM,QSQ,CF(11)                                               
      REAL*8    PM2/1.876512/, PMSQ/.8803/, PHTR/.61993/
      REAL*8    V,OMEGA,XX,XPX,OMEGAW,ARG

                                                                        
C OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING                             
                                                                        
      V      = (WM**2+QSQ-PMSQ)/PM2                                     
      OMEGA  = 2.*CF(8)*V/QSQ                                           
      XX     = 1./OMEGA                                                 
      XPX    = CF(9)+CF(10)*(XX-CF(11))**2                              
      OMEGAW = (2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))                         
      ARG    = 1.-1./OMEGAW                                             
                                                                        
      SLACF2 = OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+                    
     >         CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)                  
      SLACF2 = SLACF2*XPX                                               
                                                                        
      RETURN                                                            
      END                                                               
      
C-----------------------------------------------------------------------
      
      
      REAL*8 FUNCTION B03(WM,QSQ,NC2,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(80),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,BRES
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.8803/, PM2/1.876512/, PM/.93828/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      INTEGER   NC2 !offset for coeficients (0 for h2)
      DATA      LSPIN/1,2,3,2/                                       

C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(NC2+22)+C(NC2+23)*(X-C(NC2+24))**2                                 
      PIEMSQ = (C(NC2+1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(NC2+1)) B1 = C(NC2+2)                                         
      EB1 = C(NC2+3)*(WM-C(NC2+1))                                              
      IF( (EB1.LE.25.).AND.(B1.NE.0.)) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(NC2+4)) B2 = (1.-C(NC2+2))                                    
      EB2 = C(NC2+5)*(WSQ-C(NC2+4)**2)                                          
      IF( (EB2.LE.25.0).AND.(B2.NE.0.)) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(NC2+2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(NC2+INDEX)                                            
           IF (I.EQ.1) RAM=C(NC2+INDEX)+C(NC2+18)*SQRT(QSQ)! +C(NC2+30)*QSQ**2
     >       + C(NC2+25)/(QSQ+C(NC2+26))
           RMA    = C(NC2+INDEX+1)
           IF (I.EQ.2) RAM =RAM + C(NC2+27)/(C(NC2+28)+QSQ)
           IF (I.EQ.3) THEN
              RMA=RMA*(1.D0 +C(NC2+20)/(1.D0+C(NC2+21)*QSQ))
              RAM =RAM + C(NC2+19)/(C(NC2+29)+QSQ) 
           ENDIF
           IF (I.EQ.4) RAM =RAM + C(NC2+30)/(C(NC2+31)+QSQ)
           RWD    = C(NC2+INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
                                                                        
C FORM VW2/F2                                                           
                                                                        
      B03 = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
                                                                        
      RETURN                                                            
      END             



! -*-Mode: Fortran; compile-command: "f77 -o f2glob.o -c f2glob.f"; -*-
!File E.13. F1990.FORTRN.                                               
!Reference:  L.W.Whitlow, SLAC-Report-357,                              
!            Ph.D. Thesis, Stanford University,                         
!            March 1990.                                                
!For details see file HELP.DOCUMENT.                                    
                                                                        
!Program contains 145 lines of Fortran code, of 72 characters each, with
!no subroutines.  Program requires File E.14 as input.                  
                                                                        
                                                                        
      SUBROUTINE F2GLOB_FAST(X,Q2,Target,MODEL,F2)
                                                                        
! Returns F2 and related quantities from the either the LAMBDA12 model  
! (MODEL=12), or the OMEGA9 model (MODEL=9), both of 27Jan90.           
!                                                                       
! F2 Deuterium is for average nucleon. i.e. approximately (F2p +F2n)/2  
!                                                                       
! Further, program returns uncertainty in F2 based on both statistics   
! (ST) and systematic effects due to our choice of models (SY).         
! Also, program calculates the slope d[F2]/d[logQ2] plus the statistical
! uncertainty in this slope.                                            
!                                                                       
! Best model is LAMBDA12.  SY is estimated to be the difference between 
! the two models.                                                       
!                                                                       
! Errors due to overall normalization are not included in program output
! and they are:  +/- 2.1% for hydrogen and +/- 1.7% for deuterium.      
! Systematic errors due to radiative corrections are shown in Reference 
! to be very kinematic independent, and are everywhere <.5%, and thus,  
! are ignored by this routine (also see documentation to dFRC in file   
! HELP.DOCUMENT).                                                       
!                                                                       
! Coefficients and correlation matrix elements are from File            
! E.13 F1990.MATRICES, dated 27Jan90.                                   
                                                                        
      IMPLICIT NONE                                                     
      LOGICAL FIRST/.TRUE./                                     
      REAL   X, Q2, F2,ST
      REAL   XPP,XP,Y,POLY,Q,QTH,DQ,F2TH,QUAD,SCB      
      REAL   BINDING                                               
      INTEGER MODEL, I, J, K                                            
      CHARACTER*1 TARGET                                                
      REAL*8    B(2,9),L(2,9,9)                      ! OMEGA9 variables 
      REAL*8    C(2,12),M(2,12,12)                   ! LAMBDA12 variable
      REAL*8 LIN                                       
                                                                        
                                                                        
! Model #9  27 Jan 90.                                                  
      DATA (B(1,J),J=1,9)/                        !  HYDROGEN Ci        
     >   0.7338659870D0,  11.0245522588D0,   2.6185804129D0,            
     >   4.0956321483D0,   0.1206495422D0,   1.9714128709D0,            
     >   3.8893348719D0, -14.0507358314D0,   8.8080576075D0/            
      DATA ((L(1,J,K),K=1,9),J=1,9)/              !  HYDROGEN MijDiDj   
     >  0.0006676790D0,   0.0088218048D0,  -0.0007305188D0,             
     > -0.0015980319D0,   0.0000814499D0,   0.0022889591D0,             
     > -0.0153597481D0,   0.0257681937D0,  -0.0129827203D0,             
     >  0.0088218048D0,   0.4084284036D0,   0.0479735629D0,             
     >  0.0472083864D0,   0.0007306896D0,  -0.0267770531D0,             
     >  0.0663676188D0,  -0.1319505427D0,   0.1028644511D0,             
     > -0.0007305188D0,   0.0479735629D0,   0.0141871362D0,             
     >  0.0188269696D0,  -0.0000772884D0,  -0.0209539831D0,             
     >  0.1024234116D0,  -0.1688799776D0,   0.0910043198D0,             
     > -0.0015980319D0,   0.0472083864D0,   0.0188269696D0,             
     >  0.0264316633D0,  -0.0001541384D0,  -0.0321703747D0,             
     >  0.1590906780D0,  -0.2577418883D0,   0.1356424745D0,             
     >  0.0000814499D0,   0.0007306896D0,  -0.0000772884D0,             
     > -0.0001541384D0,   0.0021536048D0,  -0.0190110257D0,             
     >  0.0585567801D0,  -0.0758507669D0,   0.0352107941D0,             
     >  0.0022889591D0,  -0.0267770531D0,  -0.0209539831D0,             
     > -0.0321703747D0,  -0.0190110257D0,   0.2220310596D0,             
     > -0.7858318126D0,   1.0974127015D0,  -0.5309260823D0,             
     > -0.0153597481D0,   0.0663676188D0,   0.1024234116D0,             
     >  0.1590906780D0,   0.0585567801D0,  -0.7858318126D0,             
     >  2.9565217889D0,  -4.2563361422D0,   2.0922424569D0,             
     >  0.0257681937D0,  -0.1319505427D0,  -0.1688799776D0,             
     > -0.2577418883D0,  -0.0758507669D0,   1.0974127015D0,             
     > -4.2563361422D0,   6.2376383315D0,  -3.1028661049D0,             
     > -0.0129827203D0,   0.1028644511D0,   0.0910043198D0,             
     >  0.1356424745D0,   0.0352107941D0,  -0.5309260823D0,             
     >  2.0922424569D0,  -3.1028661049D0,   1.5586723492D0/             
                                                                        
      DATA (B(2,J),J=1,9) /                       !  Deuterium Ci       
     >  0.6087459014D0,   8.4283440045D0,   1.8643042857D0,             
     >  3.1298831009D0,   0.1952690820D0,   0.8207482504D0,             
     >  3.2808011387D0,  -8.2972794804D0,   4.4892920417D0/             
                                                                        
      DATA ((L(2,J,K),K=1,9),J=1,9)/              !  Deuterium MijDiDj  
     >  0.0004823134D0,   0.0055128707D0,  -0.0003158223D0,             
     > -0.0008664550D0,   0.0000058824D0,   0.0013253049D0,             
     > -0.0072791640D0,   0.0109300741D0,  -0.0049461930D0,             
     >  0.0055128707D0,   0.2107333442D0,   0.0259720298D0,             
     >  0.0248189032D0,   0.0007144468D0,  -0.0145424906D0,             
     >  0.0405570442D0,  -0.0721227448D0,   0.0486265355D0,             
     > -0.0003158223D0,   0.0259720298D0,   0.0068492388D0,             
     >  0.0088813426D0,   0.0001809208D0,  -0.0091545289D0,             
     >  0.0388897684D0,  -0.0588631696D0,   0.0295266467D0,             
     > -0.0008664550D0,   0.0248189032D0,   0.0088813426D0,             
     >  0.0124007760D0,   0.0002241085D0,  -0.0138537368D0,             
     >  0.0599295961D0,  -0.0889074149D0,   0.0432637631D0,             
     >  0.0000058824D0,   0.0007144468D0,   0.0001809208D0,             
     >  0.0002241085D0,   0.0010114008D0,  -0.0090339302D0,             
     >  0.0277972497D0,  -0.0356355323D0,   0.0162553516D0,             
     >  0.0013253049D0,  -0.0145424906D0,  -0.0091545289D0,             
     > -0.0138537368D0,  -0.0090339302D0,   0.0957852750D0,             
     > -0.3188133729D0,   0.4239206981D0,  -0.1961729663D0,             
     > -0.0072791640D0,   0.0405570442D0,   0.0388897684D0,             
     >  0.0599295961D0,   0.0277972497D0,  -0.3188133729D0,             
     >  1.1017824091D0,  -1.4925639539D0,   0.6968068686D0,             
     >  0.0109300741D0,  -0.0721227448D0,  -0.0588631696D0,             
     > -0.0889074149D0,  -0.0356355323D0,   0.4239206981D0,             
     > -1.4925639539D0,   2.0479986415D0,  -0.9652124406D0,             
     > -0.0049461930D0,   0.0486265355D0,   0.0295266467D0,             
     >  0.0432637631D0,   0.0162553516D0,  -0.1961729663D0,             
     >  0.6968068686D0,  -0.9652124406D0,   0.4591313874D0/             
                                                                        
!    MODEL #12:    27Jan90.                                             
      DATA (C(1,J),J=1,12)/                !     HYDROGEN Ci            
     >  1.4168453160D0,  -0.1076464631D0,   1.4864087376D0,             
     > -5.9785594887D0,   3.5240257602D0,  -0.0106079410D0,             
     > -0.6190282831D0,   1.3852434724D0,   0.2695209475D0,             
     > -2.1790402676D0,   4.7223977551D0,  -4.3633393929D0/             
      DATA ((M(1,J,K),K=1,12),J=1,12)/     !     HYDROGEN MijDiDj       
     >  0.0014961921D0,  -0.0114525491D0,   0.0302843702D0,             
     > -0.0334635318D0,   0.0132208899D0,   0.0000371728D0,             
     > -0.0004173300D0,   0.0007986253D0,   0.0000132630D0,             
     > -0.0000712621D0,  -0.0001056593D0,   0.0004288772D0,             
     > -0.0114525491D0,   0.0967765603D0,  -0.2740561190D0,             
     >  0.3184559770D0,  -0.1307971364D0,  -0.0011246012D0,             
     >  0.0095305519D0,  -0.0155069847D0,  -0.0010495929D0,             
     >  0.0090797755D0,  -0.0200963251D0,   0.0116773587D0,             
     >  0.0302843702D0,  -0.2740561190D0,   0.8159191015D0,             
     > -0.9844443599D0,   0.4163716693D0,   0.0049245087D0,             
     > -0.0379185977D0,   0.0567662659D0,   0.0051689160D0,             
     > -0.0439817571D0,   0.0995835938D0,  -0.0638367188D0,             
     > -0.0334635318D0,   0.3184559770D0,  -0.9844443599D0,             
     >  1.2221697276D0,  -0.5286404057D0,  -0.0072551971D0,             
     >  0.0521650844D0,  -0.0735924860D0,  -0.0082081518D0,             
     >  0.0683850387D0,  -0.1551044074D0,   0.1026791211D0,             
     >  0.0132208899D0,  -0.1307971364D0,   0.4163716693D0,             
     > -0.5286404057D0,   0.2327515525D0,   0.0034606631D0,             
     > -0.0235526467D0,   0.0317074158D0,   0.0041807175D0,             
     > -0.0342135427D0,   0.0775630764D0,  -0.0522714782D0,             
     >  0.0000371728D0,  -0.0011246012D0,   0.0049245087D0,             
     > -0.0072551971D0,   0.0034606631D0,   0.0006331410D0,             
     > -0.0035750486D0,   0.0043493144D0,   0.0005207326D0,             
     > -0.0035419381D0,   0.0068329087D0,  -0.0038428417D0,             
     > -0.0004173300D0,   0.0095305519D0,  -0.0379185977D0,             
     >  0.0521650844D0,  -0.0235526467D0,  -0.0035750486D0,             
     >  0.0234071623D0,  -0.0312734982D0,  -0.0029088270D0,             
     >  0.0220336426D0,  -0.0446325428D0,   0.0252355182D0,             
     >  0.0007986253D0,  -0.0155069847D0,   0.0567662659D0,             
     > -0.0735924860D0,   0.0317074158D0,   0.0043493144D0,             
     > -0.0312734982D0,   0.0455043874D0,   0.0034940236D0,             
     > -0.0283748709D0,   0.0601210472D0,  -0.0342674110D0,             
     >  0.0000132630D0,  -0.0010495929D0,   0.0051689160D0,             
     > -0.0082081518D0,   0.0041807175D0,   0.0005207326D0,             
     > -0.0029088270D0,   0.0034940236D0,   0.0007624603D0,             
     > -0.0058108049D0,   0.0129263887D0,  -0.0087097278D0,             
     > -0.0000712621D0,   0.0090797755D0,  -0.0439817571D0,             
     >  0.0683850387D0,  -0.0342135427D0,  -0.0035419381D0,             
     >  0.0220336426D0,  -0.0283748709D0,  -0.0058108049D0,             
     >  0.0487297250D0,  -0.1154000355D0,   0.0812897233D0,             
     > -0.0001056593D0,  -0.0200963251D0,   0.0995835938D0,             
     > -0.1551044074D0,   0.0775630764D0,   0.0068329087D0,             
     > -0.0446325428D0,   0.0601210472D0,   0.0129263887D0,             
     > -0.1154000355D0,   0.2885784358D0,  -0.2128276155D0,             
     >  0.0004288772D0,   0.0116773587D0,  -0.0638367188D0,             
     >  0.1026791211D0,  -0.0522714782D0,  -0.0038428417D0,             
     >  0.0252355182D0,  -0.0342674110D0,  -0.0087097278D0,             
     >  0.0812897233D0,  -0.2128276155D0,   0.1642123699D0/             
                                                                        
      DATA (C(2,J),J=1,12)/                !     DEUTERIUM Ci           
     >  0.9483220437D0,  -0.1153382195D0,   1.8614034534D0,             
     > -4.7333791157D0,   2.3483754563D0,  -0.0651156444D0,             
     > -0.2243092198D0,   1.0850340284D0,   0.2125643792D0,             
     > -1.6872146840D0,   3.4085883231D0,  -3.2545701111D0/             
      DATA ((M(2,J,K),K=1,12),J=1,12)/     !     DEUTERIUM MijDiDj      
     >  0.0007144431D0,  -0.0055332437D0,   0.0148345485D0,             
     > -0.0166543296D0,   0.0066913067D0,  -0.0000063353D0,             
     > -0.0000313908D0,   0.0001476921D0,  -0.0000519937D0,             
     >  0.0004518877D0,  -0.0011993941D0,   0.0010410232D0,             
     > -0.0055332437D0,   0.0464241060D0,  -0.1316100281D0,             
     >  0.1539289430D0,  -0.0638038463D0,  -0.0004724619D0,             
     >  0.0037853638D0,  -0.0060936945D0,  -0.0000911765D0,             
     >  0.0007345446D0,  -0.0009520769D0,  -0.0006845386D0,             
     >  0.0148345485D0,  -0.1316100281D0,   0.3889562708D0,             
     > -0.4695521254D0,   0.1995114383D0,   0.0024459109D0,             
     > -0.0172286634D0,   0.0247997549D0,   0.0014795243D0,             
     > -0.0120036957D0,   0.0259982755D0,  -0.0151245338D0,             
     > -0.0166543296D0,   0.1539289430D0,  -0.4695521254D0,             
     >  0.5810365405D0,  -0.2518031702D0,  -0.0037864499D0,             
     >  0.0248168137D0,  -0.0334709127D0,  -0.0029341015D0,             
     >  0.0235187814D0,  -0.0525907667D0,   0.0342155275D0,             
     >  0.0066913067D0,  -0.0638038463D0,   0.1995114383D0,             
     > -0.2518031702D0,   0.1108885979D0,   0.0018333157D0,             
     > -0.0113882800D0,   0.0146376694D0,   0.0016469653D0,             
     > -0.0130947155D0,   0.0297048474D0,  -0.0201812916D0,             
     > -0.0000063353D0,  -0.0004724619D0,   0.0024459109D0,             
     > -0.0037864499D0,   0.0018333157D0,   0.0005976780D0,             
     > -0.0033294157D0,   0.0040280997D0,   0.0004270733D0,             
     > -0.0027573603D0,   0.0049156906D0,  -0.0024136903D0,             
     > -0.0000313908D0,   0.0037853638D0,  -0.0172286634D0,             
     >  0.0248168137D0,  -0.0113882800D0,  -0.0033294157D0,             
     >  0.0207148104D0,  -0.0268964589D0,  -0.0023283682D0,             
     >  0.0162308979D0,  -0.0297645179D0,   0.0142701075D0,             
     >  0.0001476921D0,  -0.0060936945D0,   0.0247997549D0,             
     > -0.0334709127D0,   0.0146376694D0,   0.0040280997D0,             
     > -0.0268964589D0,   0.0372995011D0,   0.0027664597D0,             
     > -0.0203157999D0,   0.0385356275D0,  -0.0183131702D0,             
     > -0.0000519937D0,  -0.0000911765D0,   0.0014795243D0,             
     > -0.0029341015D0,   0.0016469653D0,   0.0004270733D0,             
     > -0.0023283682D0,   0.0027664597D0,   0.0005581515D0,             
     > -0.0041387256D0,   0.0089984380D0,  -0.0059280886D0,             
     >  0.0004518877D0,   0.0007345446D0,  -0.0120036957D0,             
     >  0.0235187814D0,  -0.0130947155D0,  -0.0027573603D0,             
     >  0.0162308979D0,  -0.0203157999D0,  -0.0041387256D0,             
     >  0.0334835563D0,  -0.0777433187D0,   0.0540437564D0,             
     > -0.0011993941D0,  -0.0009520769D0,   0.0259982755D0,             
     > -0.0525907667D0,   0.0297048474D0,   0.0049156906D0,             
     > -0.0297645179D0,   0.0385356275D0,   0.0089984380D0,             
     > -0.0777433187D0,   0.1924237194D0,  -0.1418467794D0,             
     >  0.0010410232D0,  -0.0006845386D0,  -0.0151245338D0,             
     >  0.0342155275D0,  -0.0201812916D0,  -0.0024136903D0,             
     >  0.0142701075D0,  -0.0183131702D0,  -0.0059280886D0,             
     >  0.0540437564D0,  -0.1418467794D0,   0.1109342554/               
      !---------------------------------------------------------------- 
      !---------------------------------------------------------------- 
                                                                        
                                                                        
                                                                        
      i = 1                                                             
      IF (TARGET.EQ.'D') i = 2                                          
      BINDING = 1./(1.-EXP(-MIN(20.,7.7*(1./X+.93828**2/Q2-1.))))       
      IF (i.EQ.1) BINDING = 1.                                          
                                                                        
      !OMEGA9 MODEL FIRST:
      IF(MODEL.EQ.9) THEN                                        
           XPP  = (Q2+B(i,1))/(Q2/X+B(i,2))                             
           XP   = (Q2+B(i,3))/(Q2/X+B(i,4))                             
           Y    = 1.-XP                                                 
           POLY = B(i,5)*Y**3+B(i,6)*Y**4+B(i,7)*Y**5+B(i,8)*Y**6+      
     >            B(i,9)*Y**7                                           
           F2  = X/XPP*BINDING*POLY                                    

      ELSEIF(MODEL.EQ.12) THEN

      !LAMBDA12 MODEL NEXT:                                             
           Y    = 1.-X                                                  
           q    = LOG(Q2)                                               
           qth  = .2+3.2*X                                              
           dq   = q-qth                                                 
           F2th = C(i,1)*Y**3+C(i,2)*Y**4+C(i,3)*Y**5+                  
     >            C(i,4)*Y**6+C(i,5)*Y**7                               
           QUAD = (C(i,6)+C(i,7)*X+C(i,8)*X**2)*dq**2                   
           LIN  = (C(i,9)+C(i,10)*X+C(i,11)*X**2+C(i,12)*X**3)*dq       
           IF (q.GT.qth) QUAD = 0.                                      
           SCB  = (1.+LIN+QUAD)                                         
           F2  = F2th*SCB*BINDING                                      
      ELSE                                                              
           WRITE(*,'('' F2GLOB: OOPS! MODEL.NE.9.AND.MODEL.NE.12'')')  
           F2=0.
           ST=-1.  
           RETURN
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               



       REAL FUNCTION F2LONU(ITARG,Q2EXT,NUEXT)
       IMPLICIT NONE
***************************************************************
* parametrizations for F2 at low nu    A V Kotwal 27dec93     *
* Original reference is F.W.Brasse et al, NP B39, 421 (1972)  *
* SLAC and DESY electroproduction and photoproduction         *
* data has been used in the fit. This fit can be used in the  *
* range 2<W<5 over a large Q2 range down to Q2=0              *
*                                                             *
* This parametrization is used in J. Franz et al,             *
* 'Inclusive Electron Scattering in the Low Q2 region'        *
* Z. Physics C 10 105-116 (1981)                              *
* in the range W>2, nu<6.5, 0.08<Q2<1.0                       *
*                                                             *
* ITARG=1 for H2, 2 for D2                                    *
* Q2EXT = Q2 in GeV^2                                         *
* NUEXT = Nu in GeV                                           *
***************************************************************
       REAL OMEGAW,NU,Q2,OMEGA,YOW,XPRIME,F2P,F2N,F2D,F2
       INTEGER ITARG
       REAL MW2,A2,B3,B4,B5,B6,B7
       PARAMETER (MW2=1.43,A2=0.42,B3=0.933,B4=-1.494,B5=9.021)
       PARAMETER (B6=-14.50,B7=6.453)
       REAL Q2EXT,NUEXT

       Q2 = Q2EXT
       NU = NUEXT
       IF (NU.LT.0.0) NU = 0.001
       IF (Q2.LT.0.0) Q2 = 0.001

       OMEGAW = (1.87654*NU+MW2)/(Q2+A2)
       OMEGA  =  1.87654*NU/Q2
       IF (OMEGA.LT.1.0) OMEGA = 1.0

       YOW    = 1.0 - 1.0/OMEGAW
       XPRIME = Q2/(1.87654*NU+0.8803)
       IF (XPRIME.GT.0.9999) XPRIME=0.9999

       F2P = B3+B4*YOW+B5*YOW**2+B6*YOW**3+B7*YOW**4
       F2P = OMEGAW*F2P*YOW**3/OMEGA
       IF (F2P.LT.0.0) F2P=0.0

       F2N = F2P*(1.0-XPRIME)
       IF (F2N.LT.0.0) F2N=0.0

       F2D = (F2N+F2P)/2.0

       IF (ITARG.EQ.1) THEN
         F2=F2P
       ELSEIF (ITARG.EQ.2) THEN
         F2=F2D
       ELSE
         F2=0.0
       ENDIF

       IF (F2.LE.0.0) F2 = 0.0

9999   CONTINUE
       F2LONU=F2
       RETURN
       END

      DOUBLE PRECISION FUNCTION F2DI15(X,Q2)
      IMPLICIT NONE
      DOUBLE PRECISION DPEMC,X,Q2,Z,AAA,BBB,CCC,ALAMB,Q2ZERO
      PARAMETER (ALAMB=0.250D0, Q2ZERO=20.D0)
      DIMENSION DPEMC(15)
*** deuterium BCDMS-like function
      DATA DPEMC /-0.0996D0,2.489D0,0.4684D0,-1.924D0,8.159D0,
     &            -10.893D0,4.535D0,
     &              0.252D0,-2.713D0,0.0254D0,0.0299D0,-1.221D0,
     &              7.5D0,-30.49D0,40.23D0/
C
      Z = 1. - X
      AAA = X**DPEMC(1)*Z**DPEMC(2)*(DPEMC(3)+DPEMC(4)*z+DPEMC(5)*Z**2
     +                        + DPEMC(6)*Z**3 + DPEMC(7)*Z**4 )
      BBB = dpemc(8)+dpemc(9)*X+dpemc(10)/(X+dpemc(11) )
      CCC = X*(dpemc(12)+dpemc(13)*X+dpemc(14)*X**2+dpemc(15)*X**3 )
      F2DI15 = AAA * ( (DLOG(Q2)    -DLOG(ALAMB**2))
     +           /(DLOG(Q2ZERO)-DLOG(ALAMB**2)) )**BBB *(1.D0+CCC/Q2)
C
      RETURN
      END
      REAL FUNCTION F2DDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(deuteron)                                        *
* Uses n/p parametrization of E665 and model of Donnachie and*
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
*                                                            *
* currently sets n/p to 1                                    *
*                                                            *
* 21may94 AVK                                                *
* set n/p = 0.94                                             *
* private communication of E665 measurement at low x and Q2  *
* from P. G. Spentzouris                                     *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)
      F2 = F2*(1.0+0.94)/2.0

      IF (F2.LT.0.0) F2 = 0.0
      F2DDL = F2

      RETURN
      END
      REAL FUNCTION F2PDL(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 6jan94                                      *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
* units of Q2 and Nu are GeV^2, GeV                          *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2EXT,NUEXT,F2,F2DOLA
      EXTERNAL F2DOLA

      F2 = F2DOLA(Q2EXT,NUEXT)

      IF (F2.LT.0.0) F2 = 0.0
      F2PDL = F2

      RETURN
      END
 
      REAL FUNCTION F2DOLA(Q2EXT,NUEXT)
      IMPLICIT NONE
**************************************************************
*     A V Kotwal 31dec93                                     *
*                                                            *
* Return F2(proton) according to the model of Donnachie and  *
* Landshoff, 'Proton Structure Function at Small Q2'         *
* May 1993 preprint M/C-TH 93/11 , DAMTP 93-23               *
*                                                            *
* The form of the function is a sum of two Regge powers, and *
* is fit to photoproduction data over a wide range of energy *
* and to NMC data for Q2<10 Gev2                             *
* Only the simple function represented by Eqns 4a and 4b are *
* used here, the refinements discussed in the paper are      *
* omitted. F2(neutron) is also omitted.                      *
*                                                            *
* units of Q2ext and Nuext are GeV^2, GeV                    *
* Q2EXT = Q2 in GeV^2                                        *
* NUEXT = Nu in GeV                                          *
**************************************************************
      REAL Q2,NU,A,SMALLA,B,SMALLB,POWER1,POWER2,TWOM,F2
      PARAMETER (A=0.324,B=0.098,SMALLA=0.562,SMALLB=0.01113)
      PARAMETER (POWER1=0.0808,POWER2=-0.4525)
      PARAMETER (TWOM=1.87654)  ! twice proton mass
      REAL TERM1,TERM2,TWOMNU
      REAL Q2EXT,NUEXT

      Q2 = Q2EXT
      NU = NUEXT
      IF (NU.LT.1.0) NU = 1.0
      IF (Q2.LT.0.0) Q2 = 0.0

      TWOMNU = TWOM*NU

      TERM1= A * TWOMNU**POWER1 / (Q2+SMALLA)**(POWER1+1.0)

      TERM2= B * TWOMNU**POWER2 / (Q2+SMALLB)**(POWER2+1.0)

      F2 = Q2 * (TERM1 + TERM2)

      IF (F2.LT.0.0) F2 = 0.0
      F2DOLA = F2

      RETURN
      END
