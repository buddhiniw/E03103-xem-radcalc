	program radcalc

!hacked up version of radcor to just calculate radiative corrections
! w/o comparing to data.

	implicit none
	include 'math_physics.inc'
	include 'include.inc'
	include 'spect.inc'

C Input data values.

	real*8 ep

C Model data values.
	real*8 mod_rad		!radiated model cross section.
	real*8 mod_raw		!raw model cross section.

C General Kinematics
	real*8 e1,deg
	real*8 thrad,snsq
	real*8 e2,nu,x,q2		!temporary kinematical variables.

C Target Info
	real*8 aux(7)
	real*8 zp,zn,a_nuc,m_nuc,m_atom,teff

C Radiative correction parameters.
	real*8 dee
	real*8 tb,ta
	real*8 tgt_len,tgt_thick,theta_tar

C Function declarations.
	real*8 radiate_calc

C ============================ Executable Code =================================

!low eps
!	e1=0.843
!	ep=0.393
!	deg=66.875

!high eps
!	e1=1.645
!	ep=1.206
!	deg=26.03

!Roy
	e1=5.766
	ep=1.14
	deg=40.0

	theta_tar=0

!For Aluminum.
c	zp=13
c	zn=14
c	tgt_len=0.1		!cm
c	tgt_thick=0.263*2	!g/cm^2

!For Deuterium
	zp=1
	zn=1
	tgt_len=4.0
	tgt_thick=0.668		!g/cm^2



	a_nuc = zp + zn
	thrad = deg*d_r
	snsq = sin(thrad/2.)**2.

       call target_info(theta_tar,deg,tb,ta,zp,a_nuc,tgt_len,tgt_thick,aux,teff)

	if (a_nuc.eq.2) then
	  m_atom=2.0150
	else if (a_nuc.eq.12) then
          m_atom=12.0110
        else if (a_nuc.eq.27) then
          m_atom=26.98
        else if (a_nuc.eq.56) then
          m_atom=55.8470
        else if (a_nuc.eq.197) then
          m_atom=196.9665
	else
	  write(6,*) 'CANT FIGURE OUT M_ATOM.  GET OFF YOUR LAZY ASS AND GET RID OF THIS HARDWIRED CRAP!!'
	endif
        m_nuc = m_atom*m_amu-zp*m_e

C initialize some stuff.
	dee = 0.005	   			!For calls to RADIAT.

C Initial output to log file (for this target).
	write (6,*) e1,' GeV = Incident energy'
	write (6,*) deg,' Degrees = scattering angle'
	write (6,*) a_nuc,' = A (target)'
	write (6,*) zp,' = Z (target)'
	write (6,*) tb,' = Radiator before in radiation lengths'
	write (6,*) ta,' = Radiator after  in radiation lengths'
	write (6,*) aux(1),' = RESOL for DEEPSIG smearing (not used)'
	write (6,*) aux(2),' = E_SEP in GeV'
	write (6,*) aux(3),' = F(0)'
	write (6,*) aux(4),' = BigB'
	write (6,*) aux(5),' = a'
	write (6,*) aux(6),' = b'
	write (6,*) aux(7),' = alpha'
	write (6,*) zp,' = Number of protons in target nucleus'
	write (6,*) zn,' = Average number of neutrons in target nucleus'
	write (6,*) dee,' = Dee in GeV'
	write (6,*) ' '

	e2 = ep
	write(6,*) 'calling sigmodel_calc'
	call sigmodel_calc(e1,e2,deg,a_nuc,zp,m_nuc,aux,mod_raw,.false.)
	write(6,*) 'back from sigmodel_calc'

	write(6,*) 'calling radiate_init_calc'
	call radiate_init_calc(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	write(6,*) 'back from radiate_init_calc'

! If there is a data point at given E' (sigfact>0) then compute unradiated
! model cross section [nb/(MeV-sr)].  Then compute radiated model using
! peaking approximation. The radiative correction at this iteration is sig/rad.
! Cross section = counts * sigfact

	    e2 = ep
	write(6,*) 'calling radiate_calc'
	    mod_rad = radiate_calc(e2,.false.)
	write(6,*) 'back from radiate_calc'


! Banner to log file.
	write (6,*)' '
	write (6,*)'  e1     e2     q2     x    sigraw     sigcor  rad_corr'
	write (6,*) ' '

! Cross section results to log file.
	e2 = ep
	nu = e1 - e2
	q2 = 4.*e1*e2*snsq
	x = q2/(2.*m_p*nu)
	write (6,1003) e1,e2,q2,x,mod_raw,mod_rad,mod_raw/mod_rad

C End of file jump destinations.

992	continue
	stop


C =========================== Format Statements ================================

1001	format(a)
1002	format(1x,f6.3,x,f9.3,x,f9.3,x,e12.4,x,e12.4,x,e12.4,x,e12.4)
1003	format(f6.3,3f7.3,2e10.3,f10.4)
1004	format(' iteration = ',i3,x,' chisqr = ',e10.4,x,' oldchisqr = ',e10.4)
1005	format(f9.6,2x,f10.4,2x,1f10.4)
1006	format(f6.3,3e12.4)
	end
