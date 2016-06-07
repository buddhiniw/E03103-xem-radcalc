	subroutine target_info(theta_tar,theta,tb,ta,z,a,tgt_len,tgt_thick,aux,teff)
C+______________________________________________________________________________
!
! TARGET_INFO - Calculate material thickness before and after scattering.
!
! ARGUMENTS:
!
!   THETA_TAR:  (R*8):  Target angle - 0 deg is normal to beam, increase CCW.
!			(i.e. increasing angle faces SOS).
!   THETA:	(R*8):	Scattering angle in Degrees (Input).
!   TB:		(R*8):	Thickness BEFORE scatter, in radiation lengths.
!   TA:		(R*8):	Thickness AFTER  scatter, in radiation lengths.
!   Z:		(R*8):	'Z' of target nucleus.
!   A:		(R*8):	'A' of target nucleus.
!   tgt_len:    (R*8):  target length in cm.
!   M_TGT:	(R*8):	Mass in GeV/c^2 of target nucleus.
!   AUX:	(R*8):	Array of length 10 containing:
!			aux(1) = RESOL parameter for DEEPSIG smearing (not used).
!			aux(2) = Separation energy in GeV of target nucleus.
!			aux(3) = F(0)  parameter for F(y) function
!			aux(4) = BigB  parameter for F(y)
!			aux(5) = a     parameter for F(y)
!			aux(6) = b     parameter for F(y)
!			aux(7) = alpha parameter for F(y)
!
! AUX definitions for old f(y) model:
!		!	aux(1) = Mass of recoiling nucleon in GeV/c^2.
!		!	aux(2) = Separation energy in GeV of target nucleus.
!		!	aux(3) = Y0 parameter for F(y) function
!		!	aux(4) = B  parameter for F(y)
!		!	aux(5) = C  parameter for F(y)
!		!	aux(6) = D  parameter for F(y)
!		!	aux(7) = RESOL parameter for DEEPSIG smearing.
!
C-______________________________________________________________________________

	implicit none

C Get various constants.

	include 'math_physics.inc'
	include 'spect.inc'

C Declare arguments.

	real*8 theta,tb,ta,z,a,aux(7),theta_tar

C Local declarations.

	integer*4 i,j,itar
	logical found

	real*8 rho,tgt_len,tgt_thick,tgt_rl
	real*8 thrad,thrad_tar
	real*8 t1,t2,teff	!effective thickness for e+/e- pair production.

! Parameters for NE18 target calculations

	real*8 x0_al,rho_al,x0_cm_al
	real*8 X0_air,rho_air,x0_cm_air
	real*8 x0_mylar,rho_mylar,x0_cm_mylar
	real*8 x0_cm_kevlar
	parameter (X0_Al	= 24.01)
	parameter (rho_Al	= 2.70)
	parameter (X0_cm_Al	= X0_Al/rho_Al)
	parameter (X0_air	= 36.66)
	parameter (rho_air	= .001205)
	parameter (X0_cm_air	= X0_air/rho_air)
	parameter (X0_mylar	= 39.95)
	parameter (rho_mylar	= 1.39)
	parameter (X0_cm_mylar	= x0_mylar/rho_mylar)
	parameter (X0_cm_kevlar	= 55.2)		!Calculation from C.Keppel.

C Lookup table from which we find the radiation length and AUX parameters.
C DAVE G NOTE: Stuff for He3,Be,Cu mostly placeholders.
C I'm going to try and use the parameters from Atti, Faralli, and West 
C (nucl-th/9812078)

	real*8		lookup(10,10)/			!data vs. Z.

!-------------------------------------------------------------------------------
!        H     2H       3He     4He      Be      C       Al     Fe     Cu       Au
!-------------------------------------------------------------------------------
     >   1.,     1.,     2.,     2.,     4.,     6.,    13.,    26.,   29.,     79., !Z's
     >   1.,     2.,     3.,     4.,     9.,    12.,    27.,    56.,   63.,    197., !A's
     >61.28,  122.6,  65.27,  94.32,  65.19,   42.7,  24.01,  13.84, 12.86,   6.461, !R.L.s in g/cm^2.
     > 0.00,   0.14,   0.10,   0.16,   0.10,   0.25,   0.25,   0.25,  0.10,    0.25, !RESOL's
     > 0.00,0.00225,0.00549,0.02020,0.00928,0.01727, 0.0099,0.01060,0.00855,0.00693, !E_SEP's in GeV.
     > 0.0056418958,0.00899068,0.00525479,0.0047,0.00409,0.003628,0.00361,0.00367,0.0034,0.00335, !f(0)'s
     > 0.00,0.0003499,0.00129612,0.0017,0.00138707,0.0006,0.00060,0.00062,0.0009,0.0009, !BigB's
     > 0.01,0.00583193,0.00830,0.00430,0.0039427,0.0033794,0.00493,0.00460,0.00345,0.00334 !a's
     > 0.00,0.00489127,0.00828328,0.0087,0.00997585,0.00870,0.00600,0.00600,0.01,0.011, !b's
     > 20.0, 48.2884,   167.059,  160.0,  164.997,  170.0,  156.0,  138.0,  175.0, 175.0 !alpha's

!previous iteration
!-------------------------------------------------------------------------------
!        H     2H       3He     4He      Be      C       Al     Fe     Cu       Au
!-------------------------------------------------------------------------------
!     >   1.,     1.,     2.,     2.,     4.,     6.,    13.,    26.,   29.,     79., !Z's
!     >   1.,     2.,     3.,     4.,     9.,    12.,    27.,    56.,   63.,    197., !A's
!     >61.28,  122.6,  65.27,  94.32,  65.19,   42.7,  24.01,  13.84, 12.86,   6.461, !R.L.s in g/cm^2.
!     > 0.00,   0.14,   0.10,   0.16,   0.10,   0.25,   0.25,   0.25,  0.10,    0.25, !RESOL's
!     > 0.00,0.00225,0.00549,0.02020,0.00928,0.01727, 0.0099,0.01060,0.00855,0.00693, !E_SEP's in GeV.
!     > 0.0056418958,0.00914,0.0053,0.00445,0.0042,0.00358,0.00361,0.00367,0.00367,0.0033, !f(0)'s
!     > 0.00,0.00025,0.0008,0.0018,0.00140,0.0006,0.00060,0.00062,0.0011,0.0009, !BigB's
!     > 0.01,0.00525,0.00710,0.00460,0.0042,0.0037,0.00493,0.00460,0.0040,0.0038, !a's
!     > 0.00,0.00470,0.00800,0.00900,0.0100,0.00870,0.00600,0.00600,0.0095,0.0095, !b's
!     > 20.0, 45.0,   83.0,  167.0,  165.0,  166.0,  156.0,  138.0,  170.0, 170.0/ !alpha's
!-------------------------------------------------------------------------------
!        H     2H       3He     4He      Be      C       Al     Fe     Cu       Au
!-------------------------------------------------------------------------------
!     >   1.,     1.,     2.,     2.,     4.,     6.,    13.,    26.,   29.,     79., !Z's
!     >   1.,     2.,     3.,     4.,     9.,    12.,    27.,    56.,   63.,    197., !A's
!     >61.28,  122.6,  65.27,  94.32,  65.19,   42.7,  24.01,  13.84, 12.86,   6.461, !R.L.s in g/cm^2.
!     > 0.00,   0.14,   0.10,   0.16,   0.10,   0.25,   0.25,   0.25,  0.10,    0.25, !RESOL's
!     > 0.00,0.00225,0.00549,0.02020,0.00928,0.01727, 0.0099,0.01060,0.00855,0.00693, !E_SEP's in GeV.
!     > 0.0056418958,0.00914,0.00628,0.00445,0.00391,0.00358,0.00361,0.00367,0.00367,0.00367, !f(0)'s
!     > 0.00,0.00025,0.00033,0.00065,0.00060,0.00057,0.00060,0.00062,0.00062,0.00062, !BigB's
!     > 0.01,0.00610,0.00710,0.00680,0.00574,0.00510,0.00493,0.00460,0.00460,0.00460, !a's
!     > 0.00,0.00600,0.00600,0.00600,0.00600,0.00600,0.00600,0.00600,0.00600,0.00600, !b's
!     > 20.0, 45.0,   83.0,  167.0,  166.0,  166.0,  156.0,  138.0,  138.0, 138.0/ !alpha's

!-------------------------------------------------------------------------------
!        H     2H       3He     4He      Be      C       Al     Fe     Cu       Au
!-------------------------------------------------------------------------------
cdg     >   1.,     1.,     2.,     2.,     4.,     6.,    13.,    26.,   29.,     79., !Z's
cdg     >   1.,     2.,     3.,     4.,     9.,    12.,    26.,    56.,   63.,    197., !A's
cdg     >61.28,  122.6,  65.27,  94.32,  65.19,   42.7,  24.01,  13.84, 12.86,   6.461, !R.L.s in g/cm^2.
cdg     > 0.00,   0.14,   0.10,   0.16,   0.10,   0.25,   0.25,   0.25,  0.10,    0.25, !RESOL's
cdg     > 0.00,0.00225,   0.00,   0.00,   0.00,0.01727, 0.0099,0.01060,  0.00, 0.00693, !E_SEP's in GeV.
cdg     > 0.00,  0.010,  0.010,  0.010,  0.010, 0.0033, 0.0030, 0.0028, 0.0030, 0.0025, !f(0)'s
cdg     > 0.00,0.00101,  0.001,0.00101,0.00040,0.00040,0.00040,0.00040,0.00040,0.00040, !BigB's
cdg     > 0.00,0.00751,0.00751,0.00751,0.00751,0.00388,0.00388,0.00388,0.00388,0.00388, !a's
cdg     > 0.00,0.00963,0.00963,0.00963, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, !b's
cdg     > 0.00,     45.,  140.,   140.,   140.,   140.,   140.,   140.,    140.,   140./ !alpha's

!!-------------------------------------------------------------------------------
!!        H       2H       He        C       Al       Fe       Au
!!-------------------------------------------------------------------------------
!     >   1.,      1.,      2.,      6.,     13.,     26.,     79., !Z's
!     >   1.,      2.,      4.,     12.,     26.,     56.,    197., !A's
!     >61.28,   122.6,   94.32,    42.7,   24.01,   13.84,   6.461, !R.L.s in g/cm^2.
!     > 0.00,    0.14,    0.16,    0.25,    0.25,    0.25,    0.25, !RESOL's
!     > 0.00, 0.00225,    0.00, 0.01727,  0.0099, 0.01060, 0.00693, !E_SEP's in GeV.
!     > 0.00,   0.010,   0.010,  0.0040,  0.0034,  0.0034,  0.0033, !f(0)'s
!     > 0.00, 0.00101, 0.00101, 0.00079, 0.00053, 0.00053, 0.00083, !BigB's
!     > 0.00, 0.00751, 0.00751, 0.00410,  0.0040,  0.0040, 0.00412, !a's
!     > 0.00, 0.00963, 0.00963,  0.0100,  0.0100,  0.0100,  0.0105, !b's
!     > 0.00,     45.,     45.,    140.,    140.,    140.,    140./ !alpha's

!!-------------------------------------------------------------------------------
!!       H     2H      He      C       Al      Fe      Au
!!-------------------------------------------------------------------------------
!     >   1.,     1.,     2.,     6.,    13.,    26.,    79.,  !Z's
!     >   1.,     2.,     4.,    12.,    26.,    56.,   197.,  !A's
!     >61.28,  122.6,  94.32,   42.7,  24.01,  13.84,  6.461,  !R.L.s in g/cm^2.
!     >  m_p,    .88,    .88,  .9383,  .9383,  .9383,  .9383,  !Recoil masses (GeV/c^2)
!     > 0.00, 0.0022,   0.00,  0.030,  0.040,  0.041,  0.049,  !E_SEP's in GeV.
!     > 0.00,  .0322, .07727, 0.1288, 0.1541, 0.1586, 0.1587,  !Y0's
!     > 0.00,11.2302,11.2302, 12.905, 12.652, 12.652, 12.019,  !B's
!     > 0.00,  .1058,  .1058,  .1131,  .1131,  .1131,  .1131,  !C's
!     > 0.00,  6.381,  6.381,  8.247,  8.085,  8.085,  7.681,  !D's
!     > 0.00,   0.14,   0.16,    0.2,    0.2,    0.2,    0.2/  !RESOL's

C ============================ Executable Code =================================

C Get target specific stuff from lookup table.

	rho = tgt_thick/tgt_len
	found = .false.
	do i = 1,10				!loop over known targets.
	  if (lookup(i,1).eq.z.and.float(int(a)).eq.lookup(i,2)) then
	    found = .true.
	    itar=i
	    tgt_rl = tgt_thick/lookup(i,3)	!Divide by radiation length.
	    do j = 1,7
	      aux(j) = lookup(i,j+3)
	    enddo
	  endif
	enddo
	if (.not.found) then
	  write(6,*) 'cant find target in lookup table!'
          return			!Quit if couldn't find info.
	endif

C Compute TA and TB for target.  Target Chamber windows added at end.

	thrad = theta*d_r		!Angle in radians.
	thrad_tar = theta_tar*d_r

	if (z.ge.4) then		!Solid target.
	  tb = tgt_rl/2./abs(cos(thrad_tar))
	  ta = tgt_rl/2./abs(cos(thrad+spect_sign*thrad_tar))

	elseif (z.eq.1..or.z.eq.2) then		!Liq. 1H or 2H target.
	  tb = tgt_rl/2.
     >       + t_can_front/X0_cm_Al

	  ta = tgt_rl/2. + t_can_back/X0_cm_Al  ! this is for tuna can targets
c	  ta = tgt_rl/2.16 + t_can_back/X0_cm_Al  ! this is for tuna can targets
cdg	  if (tgt_len.le.7) then		! 4 cm target
cdg	    if (theta.lt.32.) then		! Out the end.
cdg	      ta = tgt_rl/2.
cdg     >           + t_can_back/X0_cm_Al/cos(thrad) !NOT CORRECTED FOR CURVED ENDCAP
cdg	    else				! Out the side.
cdg	      ta = can_radius*rho/lookup(itar,3)/sin(thrad)
cdg     >           + t_can_side/X0_cm_Al/sin(thrad)
cdg	    endif
cdg	  else                                  !15 cm target, always side.
cdg	      ta =  can_radius*rho/lookup(itar,3)/sin(thrad)
cdg     >           + t_can_side/X0_cm_Al/sin(thrad)
cdg	  endif
	endif

! The parameterization for the pair production cross section is cross
! section per nucleon, per effective tgt len.  teff is (5+t1)*(5+t2)/10
! where t1,t2 are twice the thicknesses in radiation lengths (%) before and
! after the center of the target. (i.e. t1=tgt_rl/cos(theta_tar))

	t1 = 100.*(2.*tb)	!double and convert to %
	t2 = 100.*(2.*ta)	!	"	"
	teff = (5.+t1)*(5+t2)/10.

! Contribution due to vacuum chamber walls, air, and entrance to spectrometer.

	ta = ta + t_exitwin/X0_cm_Al		!Al exit window.
     >          + t_air/X0_cm_air 		!air before magnets 
     >          + t_mylar/X0_cm_mylar		!mylar entrance windows.
     >          + t_kevlar/X0_cm_kevlar		!kevlar entrance windows.


	return	

1001	format(a)
	end
