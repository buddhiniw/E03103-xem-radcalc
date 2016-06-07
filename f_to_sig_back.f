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

	real*8 q4_sq,gep,gmp,gen,gmn,temp

C Declare function data types.

	logical y_calc_ok
	real*8  fy
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




C Mott cross section.

	sig_mott = cos(th/2.)/(274.072*e1*sin(th/2.)**2.)
	sig_mott = (.019733*sig_mott)**2


! Model 12 (G.K.) for form factors - maybe add Kelly later?

	call nform (12.0,q4_sq,gep,gen,gmp,gmn)
c	temp = gmp

cdgC original junk
cdg
cdg	dwdy = q3/mp
cdg	ge2 = z*gep**2 + (a-z)*gen**2
cdg	gm2 = z*gmp**2 + (a-z)*gmn**2
cdg
cdgC Calculate W2 and W1 for nucleus at rest.
cdg
cdg	w20 = (gm2*tau+ge2)/(1.+tau)
cdg	w10 = tau*gm2
cdg
cdg        call y_calc(e1,e2,theta,m1,mr,es,y,y_calc_ok)
cdg	if (.not.y_calc_ok) then
cdg	   sig = 0.
cdg	   y = 1.0E20
cdg	   goto 900
cdg	endif
cdg
cdg	pf = (1. - exp(-a/8.))*.22 + .04
cdg	pp2 = .4*pf**2
cdg
cdgC Calculate effective W2, W1. (Convection currents?)
cdg
cdg	g2 = sqrt(y**2 + pp2 + (m2)**2)
cdg	f2 = ((m1-g2-nu*y/q3)**2+(1.-nu**2/q3**2)*pp2/2.)/mp_sq
cdg
cdg	w2 = f2*w20/dwdy
cdg	w1 = (w10 + w20*pp2/(2.*mp_sq))/dwdy
cdg
cdg
cdg	fact = sig_mott*(w2 + 2*w1*tan(th/2.)**2)

CDJG OK - so here, I'm trying to reproduce the formalism used in Ciofi et al
CDJG PRC 43, p. 1155. Using this formalism, along with the parameters from their
CDJG analysis in nucl-th/9812078. Note that all this gives the "asymptotic" scaling
CDJG function (i.e. at infinite Q2). At larger -y, this increasingly disagrees
CDJG with finite Q2 results.

	f2p = (gmp-gep)/(1.0+tau)
	f1p = (tau*gmp+gep)/(1.0+tau)

	f2n = (gmn-gen)/(1.0+tau)
	f1n = (tau*gmn+gen)/(1.0+tau)

	if(a.gt.1.0) then
	   WMA = sqrt((nu+M1)**2-q3**2) ! invarient mass of gamma+MA system
	   if((WMA**2-m2**2-mr**2)**2 -4.0*m2**2*mr**2.gt.0.0) then
	      kcm = sqrt( (WMA**2-m2**2-mr**2)**2 -4.0*m2**2*mr**2 )/2.0/WMA
	   else
	      kcm=0.0
	   endif

	   kmin = (nu+M1)*abs( kcm-q3/(nu+m1)*sqrt(m2**2+kcm**2) )/WMA
	   
	   Einit = sqrt(mr**2+kmin**2)
	   Efinal = sqrt(mr**2+(kmin+q3)**2)

	   qbar2 = q3**2 - (Einit-Efinal)**2
	   taubar = qbar2/4.0/mr**2
	   
	   prefac = sig_mott/Einit/Efinal

	   sig_p = (q4_sq/q3**2)**2*( (Einit+Efinal)**2/4.0*(f1p**2+taubar*f2p**2) 
	1	                                            -q3**2/4.0*(f1p+f2p)**2 ) 
	2	+(tan(th/2.)**2 + q4_sq/2./q3**2)*qbar2/2.0*(f1p+f2p)**2
	   
	   sig_p = prefac*sig_p
	   
	   sig_n = (q4_sq/q3**2)**2*( (Einit+Efinal)**2/4.0*(f1n**2+taubar*f2n**2) 
	1	                                            -q3**2/4.0*(f1n+f2n)**2 ) 
	2	+(tan(th/2.)**2 + q4_sq/2./q3**2)*qbar2/2.0*(f1n+f2n)**2
	   sig_n = prefac*sig_n
	else
	   prefac = sig_mott*e2/e1
	   sig_p = (gep**2+tau*gmp**2)/(1.0+tau) + 2.0*tau*gmp**2*tan(th/2.)**2
	   sig_p = prefac*sig_p
	   sig_n = 0.
	endif
C Calculate scaling variable.

        call y_calc(e1,e2,theta,m1,mr,es,y,y_calc_ok)
	if (.not.y_calc_ok) then
	   sig = 0.
	   y = 1.0E20
	   goto 900
	endif

	dwdy = q3/sqrt(mp**2+q3**2+y**2+2.0*q3*y)

	if(a.gt.1.0) then
	   fact = (Z*sig_p+(A-Z)*sig_n)/dwdy
	else
	   fact = Z*sig_p/dwdy
	endif
	sig = fact*fy(y,aux(3),aux(4),aux(5),aux(6),aux(7))*1.D6 !new f(y) fit

	write(19,*) y, sig_p, sig_n,theta, e1,e2

900	continue 

	return
	end
