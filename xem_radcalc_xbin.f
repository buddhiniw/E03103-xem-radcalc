	program xem_radcalc

! Hacked up version of John's hacked up version.
! Generate lookup tables in delta and theta for sticking radiative corrections
! in 1st pass XEM analysis

!hacked up version of radcor to just calculate radiative corrections
! w/o comparing to data.

	implicit none
	include 'math_physics.inc'
	include 'include.inc'
	include 'spect.inc'

C Input data values.

	real*8 epmin,epmax

C Model data values.
	real*8 mod_rad		!radiated model cross section.
	real*8 mod_raw		!raw model cross section.
	real*8 dis_raw,qe_raw   ! dis and qe pieces
	real*8 mod_cc		! coulomb boosted raw model cross section.
	real*8 dis_cc,qe_cc     ! coulomb boosted dis and qe pieces
	real*8 ccor             ! coulomb correction
	real*8 deltae, f1cc     ! boost and focus factor for coulomb correction.
C General Kinematics
	real*8 e1,deg
	real*8 e1cc, pitcc      ! boosted energies
	real*8 thrad,snsq
	real*8 e2,nu,x,q2,x_min		!temporary kinematical variables.
C Loop variables
	real*8 dth,thmin,thit
	real*8 dp,pmin,pit

	integer ith,ip,npbin, ix, nxbin
C Target Info
	real*8 aux(7)
	real*8 zp,zn,a_nuc,m_nuc,m_atom,teff

C Radiative correction parameters.
	real*8 dee
	real*8 tb,ta
	real*8 tgt_len,tgt_thick,theta_tar

C Function declarations.
	real*8 radiate_calc, sigma_peter
	real*8 x_max, dx, pix
	integer nn, na, eof

	logical just_this
	
C       ============================ Executable Code ===================
C       ==============
	
	just_this=.false.
!	if(do_exactly_this) then
				!read lines from file
!	  open(unit=11,file='manual_input.dat',status='old')
! 30	  READ(11,*,IOSTAT=EOF) e1, e2, deg, a_num
!	  if(EOF.ge.0) then
!	    manual_e1(index)=e1
!	    manual_e2(index)=e2
!	    manual_deg(index)=deg
!	    manual_anum(index)=a_num
!	  endif
!	  goto 30
!	endif
	
	
	write(6,*) '**Beam energy (GeV)'
	
	read(5,*) e1

	write(6,*) '**HMS central scattering angle (degrees)'
	read(5,*) deg

	write(6,*) '**Minimum HMS central momentum (GeV)'
	read(5,*) epmin
C Shift epmin to -12%

!       to get epmin -> 12% and th_min to get x_min, then solve for E' ! theta and x_min
C       at central
	
	x_min=4*e1*epmin*0.88*(sin((deg*d_r-0.034)/2))**2/(2*m_p*(e1
	1 -epmin*0.88))
	
	epmin=x_min*e1*m_p/(2*e1*(sin(deg*d_r/2))**2+m_p*x_min)
!	write (*,*) 'got a min ep of ', epmin


!	epmin = 0.88*epmin

	write(6,*) '**Maximum HMS central momentum (GeV)'
	read(5,*) epmax
	epmax = 1.12*epmax

	write(6,*) '**Target A (e.g. 12 for C)'
	read(5,*) a_nuc
	
	x_max=4*e1*epmax*(sin((deg*d_r+0.034)/2))**2/(2*m_p*(e1
	1 -epmax))
	
	if(x_max.gt.a_nuc) then
	  if (a_nuc.le.4.0) then
	    x_max=a_nuc
	  else
	    x_max=5.0
	  endif
	endif
	if(x_max.gt.5.0) x_max=5.0
	write (6,*) '**xmin and xmax are', x_min, x_max, a_nuc
C note that these deltae's are calculated as follows:
C deltae = (1.5)*(Z/R)*(hbar*c)*alpha*0.775
C the 0.775 comes from Aste EPJ A 26 167 (2005)
C also, all deltae's are computed for Z-1, not Z!
	if(a_nuc.eq.1) then	!Hydrogen
	   zp = 1
	   zn = 0
	   tgt_len=4.0
	   tgt_thick=0.289
	   deltae = 0.0
	elseif(a_nuc.eq.2) then !Deuterium
	   zp=1
	   zn=1
	   tgt_len=4.0
	   tgt_thick=0.655	!g/cm^2
	   deltae = 0.0
	elseif(a_nuc.eq.3) then !Helium-3
	   zp=2
	   zn=1
	   tgt_len=4.0
	   tgt_thick=0.280	!g/cm^2
	   deltae = 0.00085
	elseif(a_nuc.eq.4) then !Helium-4
	   zp=2
	   zn=2
	   tgt_len=4.0
	   tgt_thick=0.543	!g/cm^2
	   deltae = 0.0010
	elseif(a_nuc.eq.9) then !Beryllium
	   zp=4
	   zn=5
	   tgt_len=1.01206
	   tgt_thick=1.8703	!g/cm^2
	   deltae = 0.001875
	elseif(a_nuc.eq.12) then !Carbon
	   zp=6
	   zn=6
	   tgt_len=0.29434
	   tgt_thick=0.6667	!g/cm^2
	   deltae = 0.00292
	elseif(a_nuc.eq.27) then !Aluminum
	   zp=13
	   zn=14
	   tgt_len=0.1948
	   tgt_thick=0.5259	!g/cm^2
c	   tgt_thick=0.0371	!g/cm^2
	   deltae = 0.0061
	elseif(a_nuc.eq.63) then !Copper
	   zp=29
	   zn=34
	   tgt_len=0.08913
	   tgt_thick=0.7986	!g/cm^2
	   deltae = 0.0102
	elseif(a_nuc.eq.197) then !Gold
	   zp=79
	   zn=118
	   tgt_len=0.01964
	   tgt_thick=0.3795	!g/cm^2
	   deltae = 0.0199
	endif

!
!	do na=1,50
!	  slope(na)=0.
!	  offset(na)=0.
!	enddo

	
        do nn=1,197
	  fermip(nn)=.2
	  EPSN(nn)=0.01
	enddo

	fermip(197)=.264
	fermip(64)=.260
	fermip(27)=.240
	fermip(12)=.221
	fermip(4)=.160
	fermip(3)=.160D0
	fermip(2)=.160D0

	epsn(197)=.006
	epsn(64)=0.01
	epsn(27)=.0083
	epsn(12)=.016
	epsn(4)=.02
	epsn(3)=.0055
	epsn(2)=.0022
	

!	print *,'duuuude ', fermip(3), epsn(3)

	do na=1,40
!	  slope(na)=0.
!	  offset(na)=0.
	  param1(na)=0.
	  param2(na)=0.
	  param3(na)=0.
	  param4(na)=0.
	enddo


! As a function ofx, same param for all Q^2

	aa=0.093235
	bb=3.24433
	cc=-0.598669
	dd=-0.83881
c	a_nuc = zp + zn
	thrad = deg*d_r
	snsq = sin(thrad/2.)**2.
!	current_slope=slope(int(deg))
!	current_offset=offset(int(deg))
!	write(*,*) ' the theta I pass in is ', int(deg), current_slope, current_offset

C DJG Solid targets are not rotated.
	theta_tar = 0.0
!	write(*,*) 'duuuude before targetinfo ', fermip(2)
	call target_info(theta_tar,deg,tb,ta,zp,a_nuc,tgt_len,tgt_thick,aux,teff)

!	write(*,*) 'duuuude ', fermip(2)
	if (a_nuc.eq.1) then
	   m_atom=1.00794
	elseif(a_nuc.eq.2) then
	  m_atom=2.0141
	  aa=0.093235
	  bb=3.24433
	  cc=-0.598669
	  dd=-0.83881
	elseif(a_nuc.eq.3) then
	   m_atom=3.0160
	   aa  = 0.984549     
	   bb  = 0.0253108    
	   cc  = 0.000220649    
	   dd  = 0.0277503     
	elseif(a_nuc.eq.4) then
	   m_atom=4.002602
	elseif(a_nuc.eq.9) then
	   m_atom=9.012182
	else if (a_nuc.eq.12) then
          m_atom=12.0110
        else if (a_nuc.eq.27) then
          m_atom=26.98
        else if (a_nuc.eq.56) then
          m_atom=55.8470
	elseif(a_nuc.eq.63) then
	   m_atom=63.546
        else if (a_nuc.eq.197) then
          m_atom=196.9665
	else
	  write(6,*) 'CANT FIGURE OUT M_ATOM.  GET OFF YOUR LAZY ASS AND GET RID OF THIS HARDWIRED CRAP!!'
	endif
        m_nuc = m_atom*m_amu-zp*m_e

	if(m_atom.lt.2.0) then
	   m_nuc = m_p
	endif

C initialize some stuff.
	dee = 0.005	   			!For calls to RADIAT.

C Initial output to log file (for this target).
	write (6,*) '**',e1,' GeV = Incident energy'
	write (6,*) '**',deg,' Degrees = scattering angle'
	write (6,*) '**',a_nuc,' = A (target)'
	write (6,*) '**',zp,' = Z (target)'
	write (6,*) '**',tb,' = Radiator before in radiation lengths'
	write (6,*) '**',ta,' = Radiator after  in radiation lengths'
	write (6,*) '**',aux(1),' = RESOL for DEEPSIG smearing (not used)'
	write (6,*) '**',aux(2),' = E_SEP in GeV'
	write (6,*) '**',aux(3),' = F(0)'
	write (6,*) '**',aux(4),' = BigB'
	write (6,*) '**',aux(5),' = a'
	write (6,*) '**',aux(6),' = b'
	write (6,*) '**',aux(7),' = alpha'
	write (6,*) '**',zp,' = Number of protons in target nucleus'
	write (6,*) '**',zn,' = Average number of neutrons in target nucleus'
	write (6,*) '**',dee,' = Dee in GeV'
	write (6,*) '**',' '


! Banner to log file.
	write (6,*) '** '
	write (6,*) '**e1     e2   theta   q2     x    sig_dis     sig_qe   sigraw        sigcor     rad_corr   ccor'
	write (6,*) '**'


!	if (do_exactly_this) then
	  !call stuff for exact values
!	else
C       OK, now start looping do 10 MeV bins in p,+/-34 mrad in theta
	  dth = 0.001*180.0/pi
	  thmin=deg-0.034*180/pi
c       dg	thmin=deg
	  
c	dp=0.5*ep/100.0
	  dp = 0.010
	  pmin=epmin
	  npbin = int((epmax-epmin)/dp)+1 !+1 one for luck
	  
	  dx = 0.0050
	  nxbin = int((x_max-x_min)/dx)+1 !+1 one for luck
	  write(6,*) '** Number of momentum bins in this file'
!	write(6,*) npbin, nxbin, x_max, x_min, dx
	  write(6,*) nxbin
	do ith=1,69
C       do this just to get things at one angle only!
!	  do ith=1,1
!	    thit=deg
	  thit=thmin+(ith-1)*dth
c       DJG   Only need to call this once per angle (I think???)
	    call radiate_init_calc(e1,thit,tb,ta,zp,zn,m_nuc,aux,dee)
c       write(6,*) 'back from radiate_init_calc'
!       do ip=1,npbin
!       pit=pmin+(ip-1)*dp
	    do ix=1,nxbin
	      pix=x_min+(ix-1)*dx
	      pit=m_p*e1*pix/(m_p*pix+2*e1*(sin(thit*d_r/2))**2)
!       write(*,*) 'pit check', pix, pit, m_p, thit
!       write(6,*) 'calling sigmodel_calc',e1,pit,thit
!       write(*,*) 'duuuude ', fermip(2)
	      call sigmodel_calc(e1,pit,thit,a_nuc,zp,m_nuc,aux,dis_raw
	1	,qe_raw,mod_raw,.false.)
!       write(*,*) 'calling peter routine'
!       call quasiy8(e1,pit,thit,a_nuc,zp,m_atom,m_nuc,sigma_peter)
!       write(6,*) 'back from sigmodel_calc',mod_raw
!       write(*,*) 'and done' 
c       now get coulomb correction.  recipe is given in:
C       Aste et al., Eur. Phys. J. A 26, 167 (2005)
	      
	      e1cc  = e1 + deltae
	      pitcc = pit + deltae
!       write(*,*) 'calling it for cc purposes'
	      call sigmodel_calc(e1cc,pitcc,thit,a_nuc,zp,m_nuc,aux
	1	,dis_cc,qe_cc,mod_cc,.false.)
	      
	      f1cc = e1cc/e1
!       If there is a data point at given E' (sigfact>0) then compute unradiated
!       model cross section [nb/(MeV-sr)].  Then compute radiated model using
!       peaking approximation. The radiative correction at this iteration is sig/rad.
!       Cross section = counts * sigfact
	      e2 = pit
c       write(6,*) 'calling radiate_calc'
	      mod_rad = radiate_calc(e2,.false.)
c       write(6,*) 'back from radiate_calc',mod_rad
	      
!       Cross section results to log file.
	      thrad = thit*d_r
	      snsq = sin(thrad/2.)**2.
	      e2 = pit
	      nu = e1 - e2
	      q2 = 4.*e1*e2*snsq
	      x = q2/(2.*m_p*nu)
	      if(mod_raw.gt.1.0E-13.and.mod_rad.gt.1.0E-13) then
		ccor = mod_raw/mod_cc/f1cc/f1cc
!       write (6,1007) e1,e2,thit,q2,x,dis_raw,qe_raw,mod_raw,mod_rad,sigma_peter,mod_raw/mod_rad,ccor,mod_cc,dis_cc,qe_cc
		write (6,1003) e1,e2,thit,q2,x,dis_raw,qe_raw,mod_raw
	1	  ,mod_rad,mod_raw/mod_rad,ccor,mod_cc,dis_cc,qe_cc
	      else
!       write (6,1007) e1,e2,thit,q2,x,dis_raw,qe_raw,mod_raw,mod_rad,sigma_peter,10.00,1.00,10.00,1.00,1.00
		write (6,1003) e1,e2,thit,q2,x,dis_raw,qe_raw,mod_raw
	1	  ,mod_rad,10.00,1.00,10.00,1.00,1.00
	      endif
	    enddo
	  enddo
C       End of file jump destinations.
	  
 992	  continue
	  stop
!	endif


C =========================== Format Statements ================================

1001	format(a)
1002	format(1x,f6.3,x,f9.3,x,f9.3,x,e12.4,x,e12.4,x,e12.4,x,e12.4)
1003	format(f6.3,4f7.3,9e12.5)
1004	format(' iteration = ',i3,x,' chisqr = ',e10.4,x,' oldchisqr = ',e10.4)
1005	format(f9.6,2x,f10.4,2x,1f10.4)
1006	format(f6.3,3e12.4)
1007	format(f6.3,4f7.3,10e12.5)
	end
