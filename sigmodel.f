	subroutine sigmodel(e1,e2,th,a,z,m_tgt,aux,sig,modelhack)
C+______________________________________________________________________________
!
! Subroutine to compute model unradiated cross section from a combination
!   of Y-scaling and deep inelastic models. The DEEPSIG model can be smeared
!   with a gaussian in missing mass squared of width RESOL$ (one sigma).
!   If SMEAR$ is .TRUE., such smearing will occur. If SMEAR$ is .FALSE.,
!   smearing is suppresed. It should be noted that the smearing involves
!   computing a convolution integral, which significantly slows down
!   this subroutine!
!
! ARGUMENTS:
!
!   E1:		-	Incident energy in GeV.
!   E2:		- Scattered energy in GeV.
!   TH:		- Scattering angle in Degrees.
!   A:		- 'A' of nucleus.
!   Z:		- Number of protons in nucleus.
!   M_TGT:	- Mass of target nucleus in GeV/c^2.
!   M_REC:	- Mass of recoiling nucleon in GeV/c^2.
!   E_SEP:	- Separation energy for target nucleus in GeV/c^2.
!   SIG  :	- Calculated cross section in nb/(MeV-ster).
C-______________________________________________________________________________

        implicit none
	include 'math_physics.inc'

C Declare arguments.

	real*8		e1,e2,th,a,z,m_tgt,aux,sig
	logical		modelhack

C Declare locals.

	real*8		sig_qe,sig_dis,y,normfac
        real*8		thr,cs,sn,elastic_peak
	save

C If at or beyond elastic peak, return ZERO!

	thr = th*d_r
	cs = cos(thr/2.)
	sn = sin(thr/2.)
	elastic_peak = e1/(1.+2.*e1*sn**2/m_tgt)
	if (e2.ge.elastic_peak) then
	   sig = 0.D00
	   return
	endif

C Test section - Return Mott cross section (times A/25) when Ep is below
C the Elastic peak. Return 0. when at or above the elastic peak.

d	sig_qe = alpha*cs*0.019733/(2.*e1*sn**2)
d 	sig_qe = 1E6*A*sig_qe**2/(1.*2.*e1*sn**2/m_p)
d 	sig_qe = sig_qe/25.D00
d	sig = sig_qe
d
d	if (.true.) return		!Gymnastics make compiler happy.

! TEST FOR MODEL DEPENDENACE OF RADCOR.
!	if (modelhack) then
!	  call f_to_sig(e1,1.2*(e2-3.53)+3.53,th,a,z,m_tgt,aux,sig_qe,y)
!	  sig_qe=1.2*sig_qe	!20% narrower, taller. (was 10% WIDER,
!				! and no taller for first test.
!	
!!	  sig_qe=6
!	else
!	  call f_to_sig(e1,e2,th,a,z,m_tgt,aux,sig_qe,y)
!	endif

! TEST FOR MODEL DEPENDENACE OF RADCOR.
!	if (modelhack) then
!	  call dis_smear(e1,e2,th,sig_dis,a,z,.false.)
!
!!	  sig_dis=6
!	else
!	  call dis_smear(e1,e2,th,sig_dis,a,z,.true.)
!	endif

	call f_to_sig(e1,e2,th,a,z,m_tgt,aux,sig_qe,y)

C	call deepsig(e1,e2,th,a,z,m_tgt,aux,sig_dis)	!bodek deut. fit, gauss smearing
	if(a.gt.1.0) then
	   call dis_smear(e1,e2,th,sig_dis,a,z,.true.)
	else
	   call dis_smear(e1,e2,th,sig_dis,a,z,.false.)  !for Hydrogen, no smearing
	endif

	sig = sig_qe + sig_dis


C Add tweak to model at 15,23 degrees, but only if modelhack=.true.
C for radcor, want this hack, but NOT for tbincorr.

cdg	if (modelhack .and. y.le.1.e10) then
cdg	  if (abs(th-15).le.0.01) then
cdg            normfac=0.94509+0.32709*y +
cdg     >          0.30442*exp(-(y-0.10770)**2/2/0.06180**2) +
cdg     >          0.30261*exp(-(y+0.25926)**2/2/0.10962**2)
cdg          else if (abs(th-23).le.0.01 .and. z.eq.26) then
cdg	    normfac=exp ( 2.59026E-02 +0.25571*y -0.10551*y**2 
cdg     >		-4.0166*y**3 -5.6110*y**4 )
cdg	  else
cdg	    normfac=1.0
cdg	  endif
cdg
cdg	  sig = sig * normfac
cdg	endif

	return
	end
