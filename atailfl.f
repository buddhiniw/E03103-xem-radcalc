      real*8 function atailfl_139(es,ep,snsq,tb,ta,zp,zn,a,iflg)  
C     Calculate radiative tail from elastic scattering. Adapted from
C     atailfl_e139. Most changes are cosmetic to improve formatting.
C     NEVER use iflg=1. This part is broken. I leave it here in case 
C     we ever want to make this work again. 
C     DJG Feb. 2007

      implicit none

      real*8 es,ep,snsq,tb,ta,zp,zn,a
      integer iflg

      real*8 funxl_139,elastcl_139,phis,phip,btr,peaks,peakp,coeff,
     >       vac_term_s,vac_term_p,facp,facs,part1,part2,part3,abrem,
     >       atail,cost,quadmo,tm,q2,xi,omp,oms,q2s,q2p,vs,vp,eta,dlz,
     >       b, vertex_term_s, vertex_term_p,tot
      integer nlvl
      real*4 brems_139,delvac
      real*8 ddilog
      external funxl_139
      common/tal1/esx,epx,svec,pvec,sp,usq,u0,uvec,cosp,coss,tmsq,fac1  
      real*8  esx,epx,svec,pvec,sp,usq,u0,uvec,cosp,coss,tmsq,fac1  
      common/tal2/zpx,znx,ax,iflgx
      real*8 zpx,znx,ax
      integer iflgx
      real*8 etail,qtail,qe_smear_peak,qe_delta_exact,qe_delta_peak, 
     >       target_radiation,qe_soft                                           
      real*8 emsq/.26113d-6/,pi/3.1415926535d0/                           
      real*8 alpha/7.2973515d-3/,two_alph_pi/.00464564d0/                                            



c------------------------------------------------------------------------
c     calculates the radiative tail of elastic or quasi-elastic         
c     scattering using the exact expression (tsai)                      
c                                                                       
c     given es = incident energy in gev                                 
c           ep = scattered energy in gev                                
c           snsq = sin(theta/2)**2 of the scattering angle                                
c           tb = radiator "before" in radiation lengths                 
c           ta = radiator "after"  in radiation lengths                 
c           zp = atomic number (no. of protons)                         
c           zn = average number of neutrons                             
c           a= atomic weight (carbon-12 scale) h=1.00797                
c           iflg = 0 for elastic, 1 for quasi-elastic                   
c                                                                       
c     answer is in pb/(gev-sr)                                          
! version of atailfl, originally from e139/e140 code, 
! adopted for use in external 7/25/96
!------------------------------------------------------------------------------
      atailfl_139=0.d0                                                      
      if (ep.le.0.d0) return                                            
c                                                                       
c     start filling the common block for the integrand                  
      epx=ep                                                            
      esx=es                                                            
      zpx=zp                                                            
      znx=zn                                                            
      ax=a                                                              
      iflgx=iflg                                                        
c                                                                       
c     calculate target dependent parameters                             
c                                                                       
c *** note factors of 1440 & 183 have been changed to 1194 & 184.15     
c     to agree with tsai, rev. mod. phys. 46, 4 (1975), on p. 828       
c                                              25 june 83 afs           
c                                                                       
      dlz=dlog(184.15d0/zp**(1.d0/3.d0))                                
      eta=dlog(1194.d0/zp**(2.d0/3.d0))/dlz                             
      b=4.d0/3.d0*(1.d0+(zp+1.d0)/(9.d0*(zp+eta)*dlz))                  
      xi=.11d0*(tb+ta)/((zp+eta)*dlz)                                   
      tm=.938256d0                                                      
      if (iflg.eq.0) tm=a*tm/1.00797d0                                  
      tmsq=tm**2                                                        
c                                                                       
c     calculate kinematics                                              
      q2=4.d0*es*ep*snsq                                                
      oms=es-ep/(1.d0-2.d0*ep*snsq/tm)                                  
      omp=es/(1.d0+2.d0*es*snsq/tm)-ep                                  
      if (omp.le.0.) return                                             
      q2s=4.d0*(es-oms)*ep*snsq                                         
      q2p=4.d0*es*(ep+omp)*snsq                                         
      vs=oms/es                                                         
      vp=omp/(ep+omp)                                                   
c                                                                       
c     calculate qed terms and mult. photon normalization                
      tot=tb+ta                                                         
c second order term in b*tot added 5/28/87 5:07 pm- steve rock          
      fac1=(1.d0+0.5772d0*b*tot- .66*(b*tot)**2)                        
     >   -alpha/(2.d0*pi)*dlog(es/ep)**2                                
cxxx dropped this term temporarily to see if can get code to work
c     *+alpha/pi*(pi**2/6.d0-dilog(1.d0-snsq))                           
     *+alpha/pi*(pi**2/6.d0-ddilog(1.d0-snsq))                           
c     facs=fac1+2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2s/emsq))  
c     facp=fac1+2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2p/emsq))  
                                                                        
                                                                        
                                                                        
c*******8/19/87**  steve rock ***************************************   
c use the bardin calculation of vertex terms instead of only electron   
c  and tsai's vertex term.                                              
c     fac2=2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(-qsq/emsq))      
                                                                        
c  vacuum terms of bardin including 3 leptons and quarks                
      vac_term_s= 2. * delvac(q2s) !delvac was implicitly real*8 before
      vac_term_p= 2. * delvac(q2p)                                
                                                                        
c  vertex term of tsai (eq. 2.11)                                       
      vertex_term_s = two_alph_pi *(-1.d0 +.75d0 *dlog(q2s/emsq) )      
      vertex_term_p = two_alph_pi *(-1.d0 +.75d0 *dlog(q2p/emsq) )      
c7/29/96      fac2 = vac_term + vertex_term                                     
c     facs=fac1+2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2s/emsq))  
c     facp=fac1+2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2p/emsq))  
      facs = fac1 + vac_term_s + vertex_term_s                          
      facp = fac1 + vac_term_p + vertex_term_p                          
                                                                        
c*******************************************************************    
                                                                        
                                                                        
c     equivalent radiator for internal                                  
      btr=alpha/pi*(dlog(q2/emsq)-1.d0)                                 
c                                                                       
c     bremsstrahlung spectrum                                           
c$$   phis=1.d0-vs+0.75d0*vs**2                                         
c$$   phip=1.d0-vp+0.75d0*vp**2                                         
                                                                        
c installed 6/4/87 to get better value of bremstrung spectrum using     
c  tsai's rev mod physics article eq (3.83) steve                       
c     this is the old-e139.e140 version of brems which is different than
c     the one traditionally in external.  ser-7/25/96
      phis= brems_139(vs,zp)                                                
      phip= brems_139(vp,zp)                                                
c                                                                       
c     s-peak                                                            
      peaks=(tm+2.d0*(es-oms)*snsq)/(tm-2.d0*ep*snsq)*facs*             
     * elastcl_139(es-oms,snsq,zp,zn,a,iflg)*(b*tb/oms*phis                 
     * +xi/(2.d0*dmax1(oms,10.d0*xi)**2))                               
c                                                                       
c     p-peak                                                            
      peakp=facp*                                                       
     * elastcl_139(es,snsq,zp,zn,a,iflg)*(b*ta/omp*phip                     
     * +xi/(2.d0*dmax1(omp,10.d0*xi)**2))                               
c                                                                       
c     calculate quantities independent of photon angle                  
      svec=dsqrt(es**2-emsq)                                            
      pvec=dsqrt(ep**2-emsq)                                            
      cost=1.d0-2.d0*snsq                                               
      sp=es*ep-svec*pvec*cost                                           
      usq=2.d0*emsq+tmsq+2.d0*tm*(es-ep)-2.d0*sp                        
      u0=es+tm-ep                                                       
      uvec=dsqrt(u0**2-usq)                                             
      cosp=(svec*cost-pvec)/uvec                                        
      coss=(svec-pvec*cost)/uvec                                        
c                                                                       
c     integrate over cos(thetak) in three pieces                        
c     making certain that the s and p peaks are found                   
                                                                        
      part1=quadmo(funxl_139,-1.d0,cosp,1.d-4,nlvl,.false.)                         
      part2=quadmo(funxl_139,cosp,coss,1.d-4,nlvl,.false.)                          
      part3=quadmo(funxl_139,coss,1.d0,1.d-4,nlvl,.false.)                          
                                                                        
c$$   part1=dcadre(funxl,-1.d0,cosp,0.d0,1.d-4,error,ier)               
c$$   if((ier.ne.0).and.(ier.ne.65).and.(ier.ne.66))                    
c$$  > write(6,'('' ** dcadre integration:part1 in atailfl; ier='',i4)')
c$$  >  ier                                                             
c$$   part2=dcadre(funxl,cosp,coss,0.d0,1.d-4,error,ier)                
c$$   if((ier.ne.0).and.(ier.ne.65).and.(ier.ne.66))                    
c$$  > write(6,'('' ** dcadre integration:part2 in atailfl; ier='',i4)')
c$$  >  ier                                                             
c$$   part3=dcadre(funxl,coss,1.d0,0.d0,1.d-4,error,ier)                
c$$   if((ier.ne.0).and.(ier.ne.65).and.(ier.ne.66))                    
c$$  > write(6,'('' ** dcadre integration:part3 in atailfl; ier='',i4)')
c$$  >  ier                                                             
                                                                        
                                                                        
      coeff=389.44d6*alpha**3*tm/pi*ep/es                               
      atail=coeff*(part1+part2+part3)                                   
c                                                                       
c     put it all together                                               
      abrem=peaks+peakp                                                 
      qe_soft=vs**(b*tb+btr)*vp**(b*ta+btr)                             
                                                                        
c     target_radiation=tarcor(es,ep,snsq,zp,a,iflg)                     
      target_radiation= 1.                                              
                                                                        
      atailfl_139=(atail*target_radiation+abrem)*qe_soft 


      return                                                            
      end                                                               

*---------------------------------------------------------------------------------*
!------------------------------------------------------------------------
                                   
c $member=funx, date=75061922, user=kce                                 
      double precision function funxl_139(cosk,modelhack) 
      implicit none                            

      real*8 cosk
      common/tal1/es,ep,svec,pvec,sp,usq,u0,uvec,cosp,coss,tmsq,fac1 
      real*8  es,ep,svec,pvec,sp,usq,u0,uvec,cosp,coss,tmsq,fac1 
      common/tal2/zp,zn,at,iflg                                         
      real*8 zp,zn,at
      integer iflg
      real*8 om,qsq,fac2,vac_term,vertex_term,fac,fw1,fw2,a,ap,bp,gnu
      real*8 x,y,term1,term21,term22,term23,term24,term25,term26,term2
      real*8 emsq/.26113d-6/,pi/3.1415926535d0/                           
      real*8 alpha/7.2973515d-3/
      real*8 two_alph_pi/.00464564d0/  
c      real*4 qsq_4,w1,w2,ff,esom_4,delvac
      real*8 qsq_4,w1,w2,ff,esom_4,delvac
      real*8 mp,tau,gep,gmp,gen,gmn
      logical modelhack ! this is just a dummy variable bacause of JRA's modified 
                        ! "quadmo" routine

      parameter (mp=0.938272)

                                                                        
c -----------------------------------------------------------------     
c 8/28/87: single preciession log - steve rock                          
c----------------------------------------------------------------- 
      om=.5d0*(usq-tmsq)/(u0-uvec*cosk)                                 
      qsq=2.d0*emsq-2.d0*sp-2.d0*om*(es-ep-uvec*cosk)  
      qsq_4 = -qsq
      esom_4 = es-om
      tau = qsq_4/4./mp**2
      if(iflg.eq.0)call nform(14.0,qsq_4,gep,gen,gmp,gmn) !use Arrington for proton
      W1 = tau*gmp**2                                              
      W2 = (gep**2+W1)/(1.+tau)   
cdjg      if(iflg.eq.0)call nuc_form_factor(qsq_4,w1,w2,ff) 
cdjg      if(iflg.eq.1) call qe_peak(qsq_4,esom_4,w1,w2)   
c$$      call astrul_139(-qsq,zp,zn,at,w1_8,w2_8,iflg)

                                                                        
c*******8/13/87**  steve rock ***************************************   
c use the bardin calculation of vertex terms instead of only electron   
c  and tsai's vertex term.                                              
c     fac2=2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(-qsq/emsq))      
c                                                                       
c  vacuum terms of bardin including 3 leptons and quarks                
      vac_term = 2. * delvac(-qsq)                                 
c  vertex term of tsai (eq. 2.11)                                       
      vertex_term =two_alph_pi *(-1.+.75 *alog(-sngl(qsq)/sngl(emsq)) ) 
      fac2 = vac_term + vertex_term                                     
c*******************************************************************    
                                                                        
                                                                        
      fac=fac1+fac2                                                     
      fw1=fac*w1                                                        
      fw2=fac*w2                                                        
      a=om*(ep-pvec*cosp*cosk)                                          
      ap=om*(es-svec*coss*cosk)                                         
                                                                        
c      tsai on page 40 of slac-pub-898 (1971) says that an uncertainty  
c      of zero / zero occurs when a = a' in his equation a.24, the      
c      integrand of which is calculated by this function.  to avoid this
c                                                                       
      funxl_139= 0.d0                                                       
      if (a.eq.ap) return                                               
c                                                                       
      bp=-om*pvec*dsqrt((1.d0-cosp**2)*(1.d0-cosk**2))                  
      gnu=1.d0/(ap-a)                                                   
      x=dsqrt(a**2-bp**2)                                               
      y=dsqrt(ap**2-bp**2)                                              
      term1=(a/x**3+ap/y**3)*emsq*(2.d0*emsq+qsq)+4.d0                  
     * +4.d0*gnu*(1.d0/x-1.d0/y)*sp*(sp-2.d0*emsq)                      
     * +(1.d0/x-1.d0/y)*(2.d0*sp+2.d0*emsq-qsq)                         
      term21=2.d0*es*(ep+om)+.5d0*qsq                                   
      term22=2.d0*ep*(es-om)+.5d0*qsq                                   
      term23=2.d0*es*ep-sp+om*(es-ep)                                   
      term24=2.d0*(es*ep+es*om+ep**2)+.5d0*qsq-sp-emsq                  
      term25=2.d0*(es*ep-ep*om+es**2)+.5d0*qsq-sp-emsq                  
      term26=emsq*(sp-om**2)+sp*term23                                  
      term2=-a*emsq/x**3*term21-ap*emsq/y**3*term22                     
     * -2.d0+2.d0*gnu*(1.d0/x-1.d0/y)*term26+1.d0/x*term24              
     * -1.d0/y*term25                                                   
      funxl_139=(fw2*term2+fw1*term1)*om/(qsq**2*(u0-uvec*cosk)) 
      if(funxl_139.gt.1.d15) then
       bp= bp*(1.+1.d-20           )
      endif
      return                                                            
      end                              

!---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ELASTCL_139(E,SNSQ,ZP,ZN,A,IFLG)
      IMPLICIT NONE
      REAL*8 E,SNSQ,ZP,ZN,A
      INTEGER IFLG            
      REAL*8 TM,RECOIL,EP,SIGM
      REAL QSQ_4,W1,W2,FF,E_4
      REAL*8 ALPHA/.00729735D0/                                           
      real*8 tau,gep,gen,gmp,gmn
                                                                       
C     GIVEN INCIDENT ELECTRON ENERGY E IN GEV                           
C           AND SNSQ OF SCATTERING ANGLE                                
C           ZP IS ATOMIC NUMBER (NO. OF PROTONS)                        
C           ZN IS AVERAGE NUMBER OF NEUTRONS                            
C           A IS ATOMIC WEIGHT (CARBON-12 SCALE) H=1.00797              
C           IFLG = 0 FOR ELASTIC, 1 FOR QUASI-ELASTIC                   
C     RETURNS THE CROSS SECTION FOR ELASTIC OR QUASI ELASTIC SCATTERING 
C     IN PB/SR (1.E-36 CM**2/SR)                                        

      TM=.938272D0                                                      
      IF (IFLG.EQ.0) TM=TM*A/1.00797D0                                  
      RECOIL=1.D0+2.D0*E/TM*SNSQ                                        
      EP=E/RECOIL                                                       
      QSQ_4=4.D0*E*EP*SNSQ                                                
      E_4 = E
      SIGM=(19732.D0*ALPHA/(2.D0*E*SNSQ))**2*(1.D0-SNSQ)/RECOIL         
C$$      CALL ASTRUL_139(QSQ,ZP,ZN,A,W1,W2,IFLG)     
cdjg      IF(IFLG.EQ.0)CALL NUC_FORM_FACTOR(QSQ_4,W1,W2,FF)
      tau = qsq_4/4./tm**2
      if(iflg.eq.0)call nform(14.0,qsq_4,gep,gen,gmp,gmn)
      W1 = tau*gmp**2                                              
      W2 = (gep**2+W1)/(1.+tau)   
cjd      IF(IFLG.EQ.1) CALL QE_PEAK(QSQ_4,E_4,W1,W2)   

      ELASTCL_139=SIGM*(W2+2.D0*W1*SNSQ/(1.D0-SNSQ))                        
      RETURN                                                            
      END                                                               

C=======================================================================
                                                                        
      REAL*8 FUNCTION DELVAC(T)                                                
                                                                        
C VACUUM POLARIZATION FOR LEPTONS AND HADRONS, includes terms for       
C electrons, muons, taus, and quarks.  The expression for small mass    
C reduces to Tsai's expression (2.10) in the slac pub.                  
      IMPLICIT NONE
      REAL*8 T, AMF2(9), COL(9), RTMF2,RT,ALT,AI34,SUML,A,B,C,SUMH
      INTEGER IC
      REAL*8  EPST / 1.E-3 /                                         
      REAL*8     ALPH_PI /.0023229 /                                    
      DATA       AMF2 / .26112E-6,.011164,3.1827,.0064,.0064,2.25,      
     +                  .09,1024.,20.25 /                               
      DATA       COL  / 3*1.,6*3. /                                     
                                                                        
C For newer treatment of vacuum correction, comment out the following   
C three lines.  These lines are exactly (2.10) of Tsai's Slac Pub '74:  
                                                                        
C     EM     = .000511                                                  
C     DELVAC = ALPH_PI*(-5./9.+ALOG(T/EM**2)/3.)                        
C     RETURN                                                            
                                                                        
      SUML = 0.      
                                                   
      DO 121 IC=1,3                                                     
           RTMF2 = T/AMF2(IC)                                           
           IF (RTMF2.LT.EPST) THEN                                      
                AI34 = RTMF2*(2./15.-1./70.*RTMF2)                      
           ELSE                                                         
                RT   = SQRT(1.+4./RTMF2)                                
                ALT  = RT*ALOG(RTMF2*(1.+RT)**2/4.)                     
                AI34 = -5./9.+4./3./RTMF2+(1./3.-2./3./RTMF2)*ALT       
           ENDIF                                                        
           SUML = SUML+COL(IC)*AI34                                     
121   CONTINUE                                                          
                                                                        
C HADRONIC VACUUM POLARIZATION TAKEN FROM BURKHARD, TASSO NOTE 192, 1982
                                                                        
      IF (T.LT.1.) THEN                                                 
           A = -1.345E-9                                                
           B = -2.302E-3                                                
           C =  4.091                                                   
      ELSE IF (T.LE.64.) THEN                                           
           A = -1.512E-3                                                
           B = -2.822E-3                                                
           C =  1.218                                                   
      ELSE                                                              
           A = -1.1344E-3                                               
           B = -3.0680E-3                                               
           C =  9.9992E-1                                               
      ENDIF                                                             
                                                                        
      SUMH   = -(A+B*ALOG(1.+C*T))                                      
      DELVAC = SUML*ALPH_PI+SUMH                                        
                                                                        
      RETURN                                                            
      END                                            

!-----------------------------------------------------------------------
      REAL FUNCTION BREMS_139 (Y,Z)                                              
                                                                        
C=======================================================================
C                                                                       
C  Returns value of Bremstrung probability multiplied by k (photon energ
C     normalized to 1 at Y=0                                            
C                                                                       
C  Formula from Tsai; Rev of Mod Physics 46,(1974)815,  equation 3.83   
C      Complete screening                                               
C                                                                       
C  Y = k/E : fraction of beam energy carried by photon                  
C  Z =     : atomic number of target                                    
C                                                                       
C      6/4/87 Steve Rock                                                
C=======================================================================
      IMPLICIT NONE
                                                                        
      REAL*8 Z,Y                                                        
      REAL LRAD,LRADP,FC,ZZ                                          
      REAL ALPHA/7.2973515E-3/                                          
      LOGICAL FIRST /.TRUE./                                                                  
                                                                        
C TSAI page 486 Table B.2  and eq. 3.67, 3.68                           
      IF(FIRST) THEN                                                                  
       IF(Z.EQ.1)THEN                                                    
        LRAD = 5.31                                                      
        LRADP= 6.144                                                     
       ELSEIF(Z.EQ.2) THEN                                               
        LRAD = 4.79                                                      
        LRADP=5.62                                                       
       ELSEIF(Z.EQ.3) THEN                                               
        LRAD = 4.74                                                      
        LRADP= 5.805                                                     
       ELSEIF(Z.EQ.4) THEN                                               
        LRAD = 4.71                                                      
        LRADP= 5.924                                                     
       ELSE                                                              
        LRAD = ALOG(184.15) - DLOG(Z)/3.                                 
        LRADP= ALOG(1194.) -2.*DLOG(Z)/3.                                
       ENDIF                                                             
                                                                        
C Get coulomb correction from Tsai eq. 3.3                              
       ZZ=  (Z*ALPHA)**2
       FC = 1.202 * ZZ - 1.0369 * ZZ**2 +                  
     >  1.008 *ZZ**3/(1.+ ZZ)
       FIRST=.FALSE.
      ENDIF                                
                                                                        
                                                                        
C  Tsai eq.(3.83) normalized to be 1 at Y=0                             
      BREMS_139 = (1 -Y +.75 *Y**2) +                                       
     >         (1.-Y)/12. *(Z+1.)/(Z *(LRAD -FC) + LRADP)               
                                                                        
      RETURN                                                            
      END                                    

!==============================================================
      SUBROUTINE NUC_FORM_FACTOR(QSQ,W1,W2,FF)
!----------------------------------------------------------------------------
! Get Nuclear Form Factor from various models
!-----------------------------------------------------------

      IMPLICIT NONE
      REAL QSQ,W1,W2,FF
      COMMON    /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                       
      INTEGER IZ,IA       
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL                                              
      REAL AVGN,AVGA,AVGM,AMUM
      REAL A_ULMAR,B_ULMAR,w2_old,PM,TAU,FDEL,FSHELL,FGAUSS,FF_BESSEL
      REAL*8 GE,GM,GEP,GEN,GMP,GMN
      LOGICAL OUT_OF_RANGE
      PARAMETER (PM    = 0.93828)  

      IG=14
      TAU = QSQ/4./PM**2 
      write(6,*) 'Q2',qsq
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                                   

      write(6,*) 'squee',GEP,GMP
c      IF (iA.EQ.1) THEN                                                 
           W1 = TAU*GMP**2                                              
           W2 = (GEP**2+W1)/(1.+TAU)                                    
c      ELSEIF (iA.EQ.2) THEN  
c        IF((IDUT.GE.11).AND.(IDUT.LE.14).AND.(QSQ.LE.3.5))THEN   !Tjon fits
c           !Ulmar d2 elastic mode                                       
c           CALL DEUT_U1(IDUT-10,QSQ,A_ULMAR,B_ULMAR)                     
c           W1 = B_ULMAR/2.  ! sigma = sig_mot(A + B*tan...)
c           W2 = A_ULMAR 
cc        ELSEIF(IDUT.EQ.1) THEN ! Linda Stuart's Model installed 5/30/96
c           CALL FFD(DBLE(QSQ),GE,GM)   
c           TAU = QSQ/4./avgM**2  
c           W1 = TAU*GM**2                                              
c           W2 = (GE**2+W1)/(1.+TAU) 
c        ELSE   ! old  elastic deuterium from original code   
c           FF  = FDEL(QSQ)                                              
c           W1  = FF**2*TAU*.6667*(GMN+GMP)**2                           
c           W2  = W1+(FF*(avgN*(GEN+TAU*GMN)+GEP+TAU*GMP)/(1.+TAU))**2   
c        ENDIF
c      ELSEIF(iA.GE.3) THEN
c       W1=0.
c       OUT_OF_RANGE =.FALSE.
c       IF(NUC_MODEL.EQ.1) THEN
c         FF = FF_BESSEL(QSQ,OUT_OF_RANGE) !only for some Nuclei,feor limited Q2
c         W2 = (iZ*FF)**2    
c       ENDIF
c       IF(OUT_OF_RANGE.OR.(NUC_MODEL.EQ.0)) THEN !use if FF_BESSEL out of range
c        IF (iA.EQ.3) THEN  !3HE  ! added 5/30/96  SER
c           CALL FFHE3(DBLE(QSQ),GE,GM)
c           TAU = QSQ/4./avgM**2  
c           W1 = TAU*GM**2                                              
c           W2 = (GE**2+W1)/(1.+TAU)  
c           W2_old  = (iz*FSHELL(QSQ) )**2
c        ELSEIF (iA.LE.20) THEN
c           FF  = FSHELL(QSQ)
c           W2  = (iZ*FF)**2                                              
c        ELSE     !ia >20
c           FF  = FGAUSS(QSQ) 
c           W2  = (iZ*FF)**2                                             
c        ENDIF
c       ENDIF                  
c      ENDIF   ! iA>+3
      RETURN
      END
