	subroutine f_to_sig(e1,e2,theta,a,z,m1,aux,sig,y)
C+______________________________________________________________________________
!
! Subroutine to compute sigma (nb/sr/MeV) from Y-scaling model. This subroutine
!   is derived from the subroutine FTOSIG by R. McKeown. It uses the same
!   basic method, but calculates the scaling variable Y from a different
!   formula, in which the separation energy is used to calculate the mass
!   of the recoiling (A-1) system. F_TO_SIG also has an argument list
!   different from FTOSIG, and hence is NOT interchangeable with FTOSIG.
!
! INPUT ARGUMENTS:
!
!   E1:		(R*8) - Incident energy in GeV
!   E2:		(R*8) - Outgoing energy in GeV
!   THETA:	(R*8) - Scattering angle in degrees
!   A:		(R*8) - 'A' of target nucleus
!   Z:		(R*8) - 'Z' of target nucleus
!   M1:		(R*8) - Mass of target nucleus in GeV/c^2.
!   MR:		(R*8) - Mass of recoiling nucleon in GeV/c^2.
!   ES:		(R*8) - Separation energy in GeV.
!
! OUTPUT ARGUMENTS:
!
!   SIG:	(R*8) - Calculated cross section
!   Y:		(R*8) - Scaling variable Y in 1/GeV.
!
!   J. Arrington, 3/5/98
!
!   AUX(1) = RESOL parameter for DEEPSIG smearing (not used).
!   AUX(2) = ES (Separation energy)
!   AUX(3) = F(0)  parameter for F(y) function
!   AUX(4) = BigB  parameter for F(y)
!   AUX(5) = a     parameter for F(y)
!   AUX(6) = b     parameter for F(y)
!   AUX(7) = alpha parameter for F(y)
!
!   D. Potterveld, 9/16/86 - AUX DEFINITIONS FROM PREVIOUS F(y) MODEL>
!
!  ! AUX(1) = MR (Recoil mass)
!  ! AUX(2) = ES (Separation energy)
!  ! AUX(3) = Y0 parameter for F(y) function
!  ! AUX(4) = B  parameter for F(y)
!  ! AUX(5) = C  parameter for F(y)
!  ! AUX(6) = D  parameter for F(y)
!  ! AUX(7) = RESOL parameter for DEEPSIG smearing.
!
! D. Potterveld, 11-24-86
C-______________________________________________________________________________

	implicit none

C Math/physics constants.

	include 'include.inc'
        real*8 pi,d_r,mp,mp_sq
	parameter (pi = 3.141592654)
	parameter (d_r = pi/180.)
	parameter (mp = .93827231)		!Proton mass in GeV/c^2.
	parameter (mp_sq = mp*mp)

C arguments.

	real*8 e1,e2,theta,a,z,m1,mr,es,sig,y
	real*8 aux(7)

C Local variables.

	real*8 th,nu,q3,m2,tau,dwdy,sig_mott
	real*8 fact
	real*8 f1p,f2p,f1n,f2n,WMA,kcm,kmin,Einit,Efinal
	real*8 prefac,qbar2,taubar,sig_p,sig_n
	real*8 ge2,gm2,w20,w10,pf,pp2,g2,f2,w2,w1

C Declare arguments for call to Hoehler(nform).

	real*8 q4_sq,gep,gmp,gen,gmn,temp,corfact

C Declare function data types.

	logical y_calc_ok
	real*8  fy, x_local, my_frac,x_high,x_low
	save

C ============================ Executable Code =================================

!!	mr = aux(1)	  		!Copy (R*8) args to int. storage
	mr = mp
	es = aux(2)

C Compute q^2 and other items.
                            
	th = theta * d_r			!Angle in radians
  	q4_sq = 4.*e1*e2*sin(th/2.)**2 		!4-mom. transf. squared.
	nu = e1-e2				!Energy loss
	q3 = sqrt(q4_sq + nu**2)		!3-momentum transfer
	m2 = m1 - mr + es			!Mass of A-1 system.
      	tau = q4_sq/(4.*mp_sq)

	x_local=q4_sq/(2*mr*nu)


C Mott cross section.

	sig_mott = cos(th/2.)/(274.072*e1*sin(th/2.)**2.)
	sig_mott = (.019733*sig_mott)**2

!- rewritten starting here--------------------------------------

        call y_calc(e1,e2,theta,m1,mr,es,y,y_calc_ok)
	if (.not.y_calc_ok) then
	   sig = 0.
	   y = 1.0E20
	   goto 900
	endif

	if(a.gt.1.0) then
	   call sig_bar_df_dble(e1, e2, theta, y, 0.0, sig_p, sig_n)
	   sig_p=sig_p/1000.
	   sig_n=sig_n/1000.

!	   write(18,*) y, sig_p, sig_n, theta, e1, e2
	else
          call nform (12.0,q4_sq,gep,gen,gmp,gmn)
	   prefac = sig_mott*e2/e1
	   sig_p = (gep**2+tau*gmp**2)/(1.0+tau) + 2.0*tau*gmp**2*tan(th/2.)**2
	   sig_p = prefac*sig_p
	   sig_n = 0.
	endif
C Calculate scaling variable.

	dwdy = q3/sqrt(mp**2+q3**2+y**2+2.0*q3*y)


	if(a.gt.1.0) then
	   fact = (Z*sig_p+(A-Z)*sig_n)/dwdy
	else
	   fact = Z*sig_p/dwdy
	endif
	sig = fact*fy(y,aux(3),aux(4),aux(5),aux(6),aux(7), A)*1.D6 !new f(y) fit
	if(x_local.gt.3.0) then
	  x_local=3.0
	endif

	if(a.eq.2) then
	  x_low=1.1
	  x_high=1.2
	endif

	if(a.gt.2) then
	  if((a.eq.63).or.(a.eq.4).or.(a.eq.197)) then
	    x_low=1.2
	    x_high=1.4
	  else
	    x_low=1.4
	    x_high=1.6
	  endif
	endif
	
!	sig=sig/1.08944    !normalization from integration
	if((x_local.ge.(x_low))) then !.and.(A.eq.2)) then
!	if(y.le.(-0.2)) then
	  corfact=(aa*exp(bb*x_local)+cc*x_local**6+dd*x_local**4)
	  if(x_local.lt.(x_high)) then
	    my_frac=(x_local-x_low)/(x_high-x_low)
	    
	    corfact=my_frac*corfact+(1.-my_frac)
	  endif
!	  corfact= aa*exp(bb*x_local)+cc*x_local**6+dd*x_local**4
!!	  if(corfact.lt.0) then
!	    sig=1.0*sig
!!	  else
	    sig=sig*corfact
!!	    write(*,*) 'corfact is ', corfact, y
!!	  endif
	endif


900	continue 

	return
	end
