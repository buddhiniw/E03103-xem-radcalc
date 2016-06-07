	real*8 function radiate_calc(ep_arg,modelhack)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Calculates the radiated cross section for quasi-elastic scattering for NE3,
!   using a peaking approximation. This subroutine is based on the subroutine
!   called RADIAT (and its supporting code), from the AMU group at SLAC. This
!   version has been re-coded to eliminate redundant calls to DILOG, and to
!   make the integrand evaluation during the EP integration more efficient,
!   in order to speed execution. The argument scheme is quite different from
!   the original, so that this is not a direct replacement for RADIAT. Also,
!   we use here a model for the Unradiated cross section developed here at
!   Caltech.
!
! AUTHOR:
!
!   D. Potterveld, Kellogg Lab, Caltech. 11/18/86
!
! MODIFICATION HISTORY:
!
!   1/15/87 - Changed all real variables to type R*8. (DHP)
!   1/05/88 - Changed code dealing with splines to use new IMSL library
!             routines. (DHP)
!
! USAGE:
!
!   1) First initialize the calculation for a given target at a given angle,
!      beam energy, and MODEL_FIX parameters via:
!
!	CALL RADIATE_INIT(E_BM,TH_SCAT,T_B,T_A,Z_P,Z_N,M_NUC,AUX,DE_E)
!
!   2) Calculate the radiated cross section via:
!
!	XSEC = RADIATE(EP_ARG)
!
!      This step may be repeated for as many ep's as you wish. Note that if
!      the model correction parameters used by MODEL_FIX are changed (this
!      happens when iterating the radiative corrections), YOU MUST CALL
!      RADIATE_INIT AGAIN, so that it can recalculate the new values of the
!      model for use in the Ep integration.
!
! ARGUMENTS (All R*8):
!
!   E_BM:	Beam energy in GeV.
!   TH_SCAT:	Scattering angle in degrees.
!   T_B:	Radiation lengths BEFORE scatter.
!   T_A:	Radiation lengths AFTER  scatter.
!   Z_P:	Number of protons in target nucleus.
!   Z_N:	Number of neutrons in target nucleus.
!   AUX:	Vector of seven elements:
!			aux(1) = RESOL parameter for DEEPSIG smearing (not used).
!			aux(2) = Separation energy in GeV of target nucleus.
!			aux(3) = F(0)  parameter for F(y) function
!			aux(4) = BigB  parameter for F(y)
!			aux(5) = a     parameter for F(y)
!			aux(6) = b     parameter for F(y)
!			aux(7) = alpha parameter for F(y)
!
! AUX definitions for old f(y) model:
!		!AUX(1) = mass of recoil nucleon in GeV/c^2
!		!AUX(2) = separation energy in GeV
!		!AUX(3) = Y0 parameter for F(y) model
!		!AUX(4) = B  parameter for F(y)
!		!AUX(5) = C  parameter for F(y)
!		!AUX(6) = D  parameter for F(y)
!		!AUX(7) = RESOL smearing parameter for DEEPSIG.
!
!   DE_E:	DEE parameter in integration limits.
!   EP_ARG:	Scattered energy in GeV.
!
! FUNCTION RETURNS:
!
!   RADIATE (R*8) returns the radiated model cross section in nb/(MeV-ster).
!   RADIATE_INIT (I*4) always returns .TRUE. However, it can also be called
!   as a subroutine, as shown above.
!
C-______________________________________________________________________________

	implicit none

	include 'include.inc'
	include 'math_physics.inc'

C Declare arguments for call to radiate_init below.

	integer radiate_init_calc
	real*8  e_bm,th_scat,t_b,t_a,z_p,z_n,m_nuc,aux_arg(7),de_e
	real*8  ep_arg
	logical modelhack

C Declare locals.

	real*8  q2,esmin,esmax,epmin,epmax,dee
	real*8  b,tot,fac,btr,ta,tb,sig,cross,sig_dis,sig_qe
	real*8  dlz,eta
	real*8  nu,mmsq,temp,di_log
	integer ind,idum

C Common blocks for evaluating integrand functions FUNCA and FUNCB below.

	real*8  es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux(7)
	common/rad/ es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux

C Function definitions.

	external funca,funcb,quadmo,atailfl_139
	real*8  funca,funcb,quadmo,atailfl_139
	real*8  ddilog
	save					!NO AMNESIA PLEASE!!!

C ============================ Executable Code =================================

	ep = ep_arg
	nu = es - ep
	radiate_calc = 0.D00
	if (ep.lt.0.D00) return
	
C Calculate some kinematics.

	q2 = 4.d0*es*ep*snsq 
	r = (tm + 2.D00*es*snsq)/(tm - 2.d0*ep*snsq)
	esmin = tm*ep/(tm-2.d0*ep*snsq)
	epmax = tm*es/(tm+2.d0*es*snsq)
	esmax = es - r*dee 
	epmin = ep + dee 

	if (epmin.gt.epmax) return 

C Calculate mult. photon and qed terms.

	temp = log(q2/m_e/m_e)
	fac1 = 1.d0 + 0.5772d0*b*tot - alpha/(2.d0*pi)*log(es/ep)**2 +
     >	alpha/pi*(pi**2/6.d0 - di_log)
	fac = fac1 + 2.d0*alpha/pi*(-14.d0/9.d0 + 13.d0/12.d0*temp)

C Equivalent radiator.

	btr = (alpha/pi)*(temp - 1.d0)
	bta = b*ta + btr 
	btb = b*tb + btr 

c	bta = b*ta
c	btb = b*tb 

C Soft photon piece.

	sig = 0.D00
!	write(6,*) 'calling sigmodel_calc from radiate_calc', es, ep, q2
	call sigmodel_calc(es,ep,th,zp+zn,zp,tm,aux,sig_dis,sig_qe,sig,modelhack)
	mmsq = m_p*m_p + 2.*m_p*nu - q2
!!!	call model_fix(-mmsq,sig)
	cross = sig

 	radiate_calc = (r*dee/es)**btb*(dee/ep)**bta*fac*cross*
     >	(1.d0 - xi/((1.d0 - (btb + bta))*max(dee,10.D00*xi))) 

C Add in integrals. (The real CPU sink!)

	radiate_calc = radiate_calc + quadmo(funca,esmin,esmax,1.d-3,idum,modelhack) +
     >	quadmo(funcb,epmin,epmax,1.d-3,idum,modelhack)

C If doing hydrogen, include the elastic tail.
	if(zp+zn.lt.1.5) then
	   radiate_calc = radiate_calc + 1.0d-6*atailfl_139(es,ep,snsq,ta,tb,zp,zn,1.00797,0)
	endif




	return 

	entry radiate_init_calc(e_bm,th_scat,t_b,t_a,z_p,z_n,m_nuc,aux_arg,de_e)
C+______________________________________________________________________________
!
! Initialize calculation.
C-______________________________________________________________________________

C Store arguments in the common block or local storage.

	es = e_bm
	th = th_scat
	snsq = sin(th*pi/360.0D00)**2	! (Sin(th/2))^2.
	tb = t_b 
	ta = t_a
	zp = z_p
	zn = z_n
	tm = m_nuc
	do ind = 1,7
	   aux(ind) = aux_arg(ind)
	enddo
	dee = de_e

C Calculate target dependent parameters. 

	tot = tb + ta
	dlz = dlog(184.15d0/zp**(1.d0/3.d0))
	eta = dlog(1194.d0/zp**(2.d0/3.d0))/dlz
	b = 4.d0/3.d0*(1.d0+(zp+1.d0)/(9.d0*(zp+eta)*dlz))
	xi = 0.11d0*(tb+ta)/((zp+eta)*dlz)

C Calculate di_log.

	di_log = ddilog(1.D00-snsq)

	radiate_init_calc = .true.
	return

	end

C ##############################################################################

	real*8 function funca(es1,modelhack)
C+______________________________________________________________________________
!
! Calculates the integrand for es1 (R*8) integration in RADIATE.
!
C-______________________________________________________________________________

	implicit none

	include 'math_physics.inc'

	real*8  es1
	real*8  sig,sig_dis,sig_qe			!Model cross section.
	real*8  nu,mmsq			!Missing mass squared.
	real*8  q2s,oms,vs,facs,phis,cross
	logical modelhack

C Common block containing auxiliary parameters.

	common/rad/ es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux
	real*8  es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux(7)

C ============================ Executable Code =================================

C Calculate kinematics for incident energy of es1.

	q2s = 4.d0*es1*ep*snsq
	oms = es - es1
	vs = oms/es

C QED terms and mult. photon normalization.

	facs = fac1 + 2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2s/m_e/m_e))

C Bremsstrahlung spectrum.

	phis = 1.d0 - vs + 0.75d0*vs**2

C Put it together.

	nu = es1 - ep
	sig = 0.D00
!	write(*,*) 'calling sigmodel from funca ', es1, ep, q2s
	call sigmodel_calc(es1,ep,th,zp+zn,zp,tm,aux,sig_dis,sig_qe,sig,modelhack)
	mmsq = m_p*m_p + 2.D00*m_p*nu - q2s
!!!	call model_fix(-mmsq,sig)
	cross = sig

	funca = (oms/(r*ep))**bta*vs**btb*(btb/oms*phis +
     >	xi/(2.d0*dmax1(oms,10.d0*xi)**2))*facs*cross

C Check for error.

	if (oms.lt.(10.D00*xi)) then
	   type 100, es,ep,snsq,es1
	   stop 'error 1 in radiate'
	endif
	return

100	format(1x,5('#'),'  ERROR: function funca ... es,ep,snsq,es1',
     >	4(e10.4,2x))

	end

C ##############################################################################

	real*8 function funcb(ep1,modelhack)
C+______________________________________________________________________________
!
! Calculates the integrand for ep1 (R*8) integration in RADIATE.
!
C-______________________________________________________________________________

	implicit none
	include 'include.inc'
	include 'math_physics.inc'

	real*8  ep1
	real*8  q2p,omp,vp,facp,phip,cross
	logical	modelhack

C Common block containing auxiliary parameters.

	common/rad/ es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux
	real*8  es,ep,snsq,r,fac1,bta,btb,zp,zn,tm,xi,th,aux(7)

	real*8 sig,sig_dis,sig_qe,nu,q2,mmsq

C Calculate kinematics for scattered energy of ep1.

	q2p = 4.d0*es*ep1*snsq
	omp = ep1 - ep
	vp = omp/ep1

C QED terms and mult. photon normalization.

	facp = fac1 + 2.d0*alpha/pi*(-14.d0/9.d0+13.d0/12.d0*dlog(q2p/m_e/m_e))

C Bremsstrahlung spectrum.

	phip = 1.d0 - vp + 0.75d0*vp**2

C Calculate sigmodel.
!	write(*,*) 'calling sigmodel from funcb ', es, ep1, q2p

	   call sigmodel_calc(es,ep1,th,zp+zn,zp,tm,aux,sig_dis,sig_qe,sig,modelhack)
	   nu = es - ep1
	   q2 = 4.D0*es*ep1*snsq
	   mmsq = m_p*m_p + 2.*m_p*nu - q2
!!!	   call model_fix(-mmsq,sig)		!Model cross section.
	   cross = sig		                !Vector.

C Put it together.

	funcb = vp**bta*(r*omp/es)**btb*(bta/omp*phip +
     >	xi/(2.d0*dmax1(omp,10.d0*xi)**2))*facp*cross

C Check for error.

	if (omp.lt.(10.D00*xi)) then
	   type 100, es,ep,snsq,ep1
	   stop 'error 2 in radiate'
	endif
	return

100	format(1x,5('#'),'Error function funcb ... es,ep,snsq,ep1 = ',
     >	4(e10.4,2x))

	end
