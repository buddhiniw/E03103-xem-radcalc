      subroutine RLT(flag, ebeam, eprime, angled, R)
      
!--------------------------------------------------------
!     this subroutine returns the value of R = sigma_L/sigma_T
!     chosen from one of three models.
!     
!     flag=0 -> R_xem = 0.32/Q^2
!     flag=1 -> R_1990: L.W. Whitlow et al., Phys. Let. B 250 193 (1990)
!     flag=2 -> R_1998: hep-ex/9808028
      include "math_physics.inc"

      IMPLICIT NONE                                                     
      
      REAL*8 ebeam, eprime, angled, R
      REAL*8 Q2, XBJ
      REAL*8 THETA

! R1990 parameters:
      REAL A90(3)  / 0.06723, 0.46714,  1.89794 /,
     >     B90(3)  / 0.06347, 0.57468, -0.35342 /,
     >     C90(3)  / 0.05992, 0.50885,  2.10807 /
C     
! R1998 parameters:
      REAL QP2,FAC,RLOG,Q2THR,R_A,R_B,R_C, D1,D2,D3,DR,DR1998
      REAL Q2_SAVE
      
!     !** Note, in S. Rock version, this is 0.2
      REAL Q2_LIMIT/.5/
      REAL A98(6) /4.8520E-02,  5.4704E-01,  2.0621E+00,
     >     -3.8036E-01,  5.0896E-01, -2.8548E-02/
      REAL B98(6) /4.8051E-02,  6.1130E-01, -3.5081E-01, 
     >     -4.6076E-01,  7.1697E-01, -3.1726E-02/
      REAL C98(6) /5.7654E-02,  4.6441E-01,  1.8288E+00,
     >     1.2371E+01, -4.3104E+01,  4.1741E+01/
      
      INTEGER FLAG   
      LOGICAL GOODFIT                                                   
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      DATA QMAX /64./

      THETA = ANGLED*D_R
      Q2  = 4*ebeam*eprime*sin(THETA/2.0)*sin(THETA/2.0)
      XBJ = Q2/2/(ebeam-eprime)/M_P

C choose your R!     
      if(flag.eq.0) then
         R = 0.32/Q2
         
!---------------------------------------------------------
         
      else if(flag.eq.1) then
         
!     R1990
C     
C     Ref: L.W.Whitlow SLAC Report 357 (1990)         and
C     L.W.Whitlow et al.: PL B250 (1990) 193
C     
C     for Q2 < 0.35 we extrapolate R as a constant with a rough error of 100 %
C     

         IF(Q2.LT.0.35) Q2=0.35
         
         FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(XBJ**2+.125**2))
         RLOG  = FAC/LOG(Q2/.04)
         Q2THR = 5.*(1.-XBJ)**5
         R_A   = A90(1)*RLOG + A90(2)/SQRT(SQRT(Q2**4+A90(3)**4))
         R_B   = B90(1)*RLOG + B90(2)/Q2 + B90(3)/(Q2**2+.3**2)
         R_C   = C90(1)*RLOG + C90(2)/SQRT((Q2-Q2THR)**2+C90(3)**2)
         R     = (R_A+R_B+R_C)/3.

c$$$         IF (Q2.GE.0.35) THEN
c$$$            Q = MIN(Q2,QMAX)
c$$$            S = .006+.03*XBJ**2
c$$$            AA= MAX(.05,8.33*XBJ-.66)
c$$$            XLOW  = .020+ABS(S*LOG(Q/AA))
c$$$            XHIGH = .1*MAX(.1,XBJ)**20/(.86**20+MAX(.1,XBJ)**20)
c$$$            D1SQUARE=XLOW**2+XHIGH**2
c$$$            D2SQUARE= ((RA-R)**2+(RB-R)**2+(RC-R)**2)/2.
c$$$            D3    = .023*(1.+.5*R)
c$$$            IF (Q2.LT.1.OR.XBJ.LT..1) D3 = 1.5*D3
c$$$            ERR    = SQRT(D1SQUARE+D2SQUARE+D3**2)
c$$$         ELSE
c$$$            ERR = R
         
!---------------------------------------------------------------------
      else
!     R1998
         
!----------------------------------------------------------------
!     XBJ      : Bjorken x                                                    
!     Q2     : Q squared in (GeV/c)**2                                      
!     R      :                                                              
!     DR     : Absolute error on R                                          
!     GOODFIT:  = .TRUE. if the XBJ,Q2 are within the range of the fit.       
!-------------------------------------------------------------------
!     Model for R, based on a fit to world R measurements. Fit performed by 
!     program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details   
!     see Reference.                                                        
!     
!     Three models are used, each model has three free parameters.  The     
!     functional forms of the models are phenomenological and somewhat      
!     contrived.  Each model fits the data very well, and the average of    
!     the fits is returned.  The standard deviation of the fit values is    
!     used to estimate the systematic uncertainty due to model dependence.  
!     
!     Statistical uncertainties due to fluctuations in measured values have 
!     have been studied extensively.  A parametrization of the statistical  
!     uncertainty of R1990 is presented in FUNCTION DR1990.                 
!     
!     
!     Each model fits very well.  As each model has its own strong points   
!     and drawbacks, R1998 returns the average of the models.  The          
!     chisquare for each fit (237 points with 6 parameters) are:  
!     ALL DATA  #PTS=237         |     XBJ<  0.07 #PTS= 28
!     FIT  #PARAM CHISQ  CHISQ/DF PROB(%) |  CHISQ  CHISQ/DF   PROB(%)
!     R1998         217.4   0.94    73.1        28.7   1.06    37.5
!     R1998a   6    219.4   0.95    69.8        28.8   1.07    37.1
!     R1998b   6    217.7   0.94    72.6        27.8   1.03    42.2
!     R1998c   6    221.9   0.96    65.5        30.9   1.15    27.4                
         
!     
!     This subroutine returns reasonable values for R for all x and for all 
!     Q2 greater than or equal to .3 GeV.                                   
!     
!     The uncertainty in R originates in three sources:                     
!     
!     D1 = uncertainty in R due to statistical fluctuations of the data 
!     as reflected in the error of the fit. 
!     It is parameterized in FUNCTION DR1990, for details see     
!     Reference.                                                   
!     
!     D2 = uncertainty in R due to possible model dependence, approxi-  
!     mated by the variance between the models.                    
!     
!     D3 = uncertainty in R due to possible epsilon dependent errors    
!     in the radiative corrections, taken to be +/- .025.  See     
!     theses (mine or Dasu's) for details.  This is copied from R1990                       
!     
!     and the total error is returned by the program:                       
!     
!     DR = is the total uncertainty in R, DR = sqrt(D1+2D) 7
!     DR is my best estimate of how well we have measured R.  At   
!     high Q2, where R is small, DR is typically larger than R.  If
!     you have faith in QCD, then, since R1990 = Rqcd at high Q2,  
!     you might wish to assume DR = 0 at very high Q2.             
!     
!     NOTE:    In many applications, for example the extraction of F2 from  
!     measured cross section, you do not want the full error in R  
!     given by DR.  Rather, you will want to use only the D1 and D2
!     contributions, and the D3 contribution from radiative        
!     corrections propogates complexely into F2.  For more informa-
!     tion, see the documentation to dFRC in HELP.DOCUMENT, or     
!     for explicit detail, see Reference.                         
!     
         
         IF(Q2.LT.Q2_LIMIT) THEN ! added 10/28/01 to match Blok et al
            Q2 = Q2_LIMIT
         ENDIF
         
         FAC = 1.+12.*Q2/(Q2+1.)*(.125**2 /(.125**2+XBJ**2))
         RLOG  = FAC/LOG(Q2/.04) !   <--- we use natural logarithms only!      
         
         
!     Model A
         QP2 = (A98(2)/(Q2**4+A98(3)**4)**.25)*(1.+A98(4)*XBJ +A98(5)*XBJ**2) !1.07   .030
         R_A   = A98(1)*RLOG +QP2*XBJ**A98(6)                                           
!     Model B
         QP2 = (1.+B98(4)*XBJ+B98(5)*XBJ**2)*(B98(2)/Q2 +B98(3)/(Q2**2+.09)) !1.06   .042
         R_B   =  B98(1)*RLOG+QP2*XBJ**B98(6)   
!     Model C
         Q2thr =C98(4)*(XBJ) +C98(5)*XBJ**2 +C98(6)*XBJ**3  
         QP2 =  C98(2)/SQRT((Q2-Q2thr)**2+C98(3)**2) 
         R_C   =  C98(1)*RLOG+QP2    
         
         R     = (R_A+R_B+R_C)/3.                                          
         
!         D1    = DR1998(XBJ,Q2)                                              
!         
!         D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)               
!         D3    = .023*(1.+.5*R)                                            
!         IF (Q2.LT.1.OR.XBJ.LT..1) D3 = 1.5*D3                       
         
         
!         DR    = SQRT(D1**2+D2**2+D3**2)                                   
         
!         GOODFIT = .TRUE.                                                  
!         IF ((XBJ.LT.0.02).OR.(Q2.LT.0.3)) GOODFIT = .FALSE.                                   
      ENDIF
      
      RETURN
      END
      
