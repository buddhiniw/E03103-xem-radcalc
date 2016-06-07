 	subroutine dis_smear(e_in,ep_in,theta,sig,anuc,znuc,smearflag)

C This routine calculates the convolution of F2_N(x/z,qsq) with f_A(z,beta)
C following the procedure of Benhar, et.al. (from Jun 1997 preprint).
C
C F2_A(x,qsq)=Integral[from x to A] of f_A(z,beta),F2_N(x/z,qsq)
C
C F2_N is the nucleon structure function for a nucleon with momentum z
C and a quark momentum fraction x/z.  Beta is defined as abs(q3)/nu.
C Here, we choose beta as a function of scattering angle only, such
C that the fit matches the data, and expected ratio of DIS/QE from other
C calculations.
C
C The f(z,beta) is taken from tables (generated by Benhar), for beta values
C of 1.0,1.1,...,1.5.  They are calculated for Nuclear Matter (?).
C
C F2_N(x/z,qsq) from bodel fits to proton and neutron.  Use N=Z=A/2 to
C get symmetric nucleus with A nucleons.  Correct for neutron excess and
C EMC effect with fit to EMC data.
C
C If smearflag=.false., return Z*f2_p + (A-Z)*f2_n - no smearing, no EMC.
C (Probably should have EMC, but I don't plan on using unsmeared).
C

	implicit none

	include 'math_physics.inc'

	real*8 e_in,ep_in,s_in,znuc,anuc,f2a,sig
	real*8 f2p,f2n,w2,w1,qsq_in,x_in,nu_in,q3_in
	real*8 xsectn,xsectp,x_test,qsq_test,xq,xtmp
	real*8 beta,amax,z,dz,fz
	real*8 e,ep,qsq,nu,sigmott
	real*8 nump,numn	!# of protons,neutrons to use in adding F2's
	real*8 theta
	real*8 r, dr
	real*8 sumef,rntd2,rautd2,emc_corr,zz
	real*8 x1,x2,frac,f2atmp,sigm

	logical smearflag
	logical initialized /.false./
	logical goodfit /.true./

	integer*4 ind,j
	integer*4 intpts
	parameter (intpts=150)

c These are the coefficients needed to correct for neutron excess
c and the emc effect. they were obtained by fitting the ratios of
c cross sections [cs(a>2)/cs(d2)] obtained in e139.the correcting
c function is a 6-th order polynomial in x.
c right now (3.6.6) ef() is not used to correct for neutron excess!

	real*8 ef(7) /-0.00136693,-0.00510425,-0.0375986,-0.0946004,
     >		    -0.122435,-0.0112751,0.406435/
	real*8 ecor(9) /163.2827,-2686.153,19152.23,-76860.48,189906.6,
     >              -295822.3,283716.9,-153174.6,35640.7/

	real*8 fz_table
	external fz_table

	real*8 slac_emc_func
	external slac_emc_func

*--------------------------------------------------------------------------

	if (.not.initialized) then
	  call fz_table_read
	  initialized=.true.
	endif

	s_in=sind(theta/2)**2
	nu_in=e_in-ep_in
	qsq_in=4*e_in*ep_in*s_in
	x_in=qsq_in/2/m_p/nu_in
	q3_in=sqrt(qsq_in+nu_in**2)

	beta=q3_in/(nu_in+0.5)

	qsq=qsq_in		!Want f2 at x=x/z, qsq=qsq

	if (.not.smearflag) then
	  call xsechd(e_in,ep_in,s_in,1,+1,f2p,w2,w1,xsectp,sigmott,x_test)
	  call xsechd(e_in,ep_in,s_in,1,-1,f2n,w2,w1,xsectn,sigmott,x_test)
	  if (abs(x_in-x_test)/x_in.gt.1e-4) then
	    write(6,*) 'xbodek,x_q=',x_test,xq
	  endif
	  f2a=znuc*f2p+(anuc-znuc)*f2n

	else 					!smearing

!(intpts)-point integration of f(z)*F2(x/z).  f(z) normalized to one, so no
! explicit check of normalization.  Truncate integral at z=5 (since f(z)=1e-10)
	  amax=min(5.,anuc)

	  if (x_in.ge.5) then		!truncated integral will give zero
	    sig=0
	    return
	  endif

	  f2a=0
	  dz=(amax-x_in)/dble(intpts)

	  do ind=1,intpts

!calculate modified kinematics to get F2(x/z,qsq).  Fix x_q=x/z, qsq, and theta
	    z=x_in+(ind-0.5)*(amax-x_in)/dble(intpts)	!nucleon has x=z
	    xq=x_in/z					!quark has x_q=x/z
	    nu=qsq/2/m_p/xq				!from x=qsq/2/m/nu
	    e=(nu+sqrt(nu**2+qsq/s_in))/2		!from qsq=2*e*(e-ep)*s
	    ep=e-nu

	    qsq_test=4*e*ep*s_in
	    x_test=qsq_test/2/m_p/nu

	   if (abs(qsq-qsq_test)/qsq.gt.1e-4.or.abs(xq-x_test)/xq.gt.1e-4) then
	      write(6,*) ' x , x_q =',x_test,xq
	      write(6,*) 'qsq,qsq_q=',qsq_test,qsq
	    endif

	    call xsechd(e,ep,s_in,1,+1,f2p,w2,w1,xsectp,sigmott,x_test)
	    call xsechd(e,ep,s_in,1,-1,f2n,w2,w1,xsectn,sigmott,x_test)
	    if (abs(xq-x_test)/xq.gt.1e-4) then
	      write(6,*) 'xbodek,x_q=',x_test,xq
	    endif

	    nump = ZNUC
	    numn = ANUC-ZNUC
	    fz=fz_table(z,beta)
	    if (anuc.lt.3.and.anuc.gt.1.) then 
	       fz=fz_table(1.3*(z-1)+1,beta) !20% narrower for deut.
	    elseif(anuc.lt.5.0.and.anuc.gt.2.) then
	       fz=fz_table(1.1*(z-1)+1,beta) !10% narrower for he3/he4.
	    endif

	    f2a=f2a + fz * (nump*f2p+numn*f2n) * dz

	    if (z.gt.1 .and. fz.lt.1.e-10) goto 200  !break loop, no more integrand

	  enddo
	endif		!if (smearflag)

200	continue


C DG The above convolution calculation does not reproduce the shape of the EMC Effect correctly
C DG Below is a correction to the shape - valid for 0.3<x<0.83. Beyond that I assume the correction
C DG is constant.

	if(x_in.lt.0.303) then
	   xtmp=0.303
	elseif(x_in.gt.0.832) then
	   xtmp=0.832
	else
	   xtmp=x_in
	endif
	
	emc_corr = 0.0
	do j=1,9
	   emc_corr = emc_corr + ecor(j)*xtmp**(j-1)
	enddo

	emc_corr = 1.0/emc_corr

C DJG Now correct for the fact that the above correction was fit for Au
C DJG For lighter nuclei, want a larger cross section, so take the ratio
C DJG of EMC(A)/EMC(Au)

	emc_corr = emc_corr * slac_emc_func(xtmp,ANUC)/slac_emc_func(xtmp,197.0d0)


        x1=.8		!x<x1 --> use xsechd (for deut.)
        x2=.9		!x1<x<x2 --> smooth transition from f2d to conv calc.
	if(anuc.gt.1.0.and.anuc.lt.3.0) then ! deuterium only
	   if (x_in.lt.x1) then
	      call xsechd(e_in,ep_in,s_in,1,+2,f2a,w2,w1,xsectp,sigm,x_test)
C---------------------
C	      call f2glob(x_in,qsq,'D',12,f2a)
C	      f2a = 2*f2a
C---------------------
	   else if (x_in.lt.x2) then
	      f2atmp = f2a
	      call xsechd(e_in,ep_in,s_in,1,+2,f2a,w2,w1,xsectp,sigm,x_test)
C---------------------
C	      call f2glob(x_in,qsq,'D',12,f2a)
C	      f2a = 2*f2a
C---------------------
	      frac = (x_in-x1)/(x2-x1)
	      f2a = f2atmp*frac + f2a*(1.-frac)
	   endif
	endif

	if(znuc.eq.1) emc_corr=1.0 !no correction for deuterium (or hydrogen)
	f2a = f2a*emc_corr

	sigmott=(19732./(2.0*137.0388*e_in*s_in))**2*(1-s_in)/1.d6

C  here's where the R = sigma_L/sigma_T goes!
C  0 = R_xem = 0.32/Q^2
C  1 = R_1990
C  2 = R_1998

C	call RLT(0,e_in,ep_in,theta,r)
C	call RLT(1,e_in,ep_in,theta,r)
C	call RLT(2,e_in,ep_in,theta,r)
C	r=0.32/qsq
	call R1990(X_IN,QSQ,R,DR,GOODFIT)
	w2=f2a/nu_in
	w1=(1.0+nu_in**2/qsq)/(nu_in*(1.0+r))*f2a
	sig=(w2+2.0*tand(theta/2)**2*w1)*sigmott

	return
	end

	function fz_table(z,beta)

	implicit none

C Table values
	integer*4 nzmax,nbmax
	parameter(nbmax=15)
	parameter(nzmax=500)
	integer*4 numb,numz(nbmax)
	real*8 bval(nbmax)
	real*8 zval(nzmax,nbmax),fzval(nzmax,nbmax)
	real*8 fznorm(nbmax)
	common/fztable/bval,fznorm,zval,fzval,numb,numz

	integer*4 ind,nb,nz
	real*8 z,beta
	real*8 mindiff
	real*8 z1,z2,fz1,fz2,logfz,w1,w2
	real*8 fz_table

!find nearest beta value
	mindiff=999.
	if (beta.ge.bval(numb)) then
	  nb=numb-1
	  w1=0
	  w2=1
	else if (beta.le.bval(1)) then
	  nb=1
	  w1=1
	  w2=0
	else
	  do ind=1,numb
	    if (beta.gt.bval(ind)) then
	      nb=ind
	      w2=(beta-bval(nb))/(bval(nb+1)-bval(nb))
	      w1=(bval(nb+1)-beta)/(bval(nb+1)-bval(nb))
	      if (abs(w1*bval(nb)+w2*bval(nb+1)-beta).gt.0.0001) then
	        write(6,*) 'w1,w2,beta,bval(nb)=',w1,w2,beta,bval(nb)
	        stop
	      endif
	    endif
	  enddo
	endif
	if (abs(w1+w2-1).gt.0.0001) then
	  write(6,*) 'nb,beta,w1+w2=',nb,beta,w1+w2
	  stop
	endif

	if (z.le.zval(1,nb)) then	!linear extrapolation of LOG(f(z)).
	  z1=zval(1,nb)
	  z2=zval(2,nb)
	  fz1=w1*fzval(1,nb)+w2*fzval(1,nb+1)
	  fz2=w1*fzval(2,nb)+w2*fzval(2,nb+1)

	else if (z.gt.zval(numz(nb),nb)) then	!linear extrapolation of LOG(f(z))
	  z1=zval(numz(nb)-1,nb)
	  z2=zval(numz(nb),nb)
	  fz1=w1*fzval(numz(nb)-1,nb)+w2*fzval(numz(nb)-1,nb+1)
	  fz2=w1*fzval(numz(nb),nb)+w2*fzval(numz(nb),nb+1)

	else
	  do nz=1,numz(nb)-1
	    if (z.ge.zval(nz,nb) .and. z.lt.zval(nz+1,nb)) then
	      z1=zval(nz,nb)
	      z2=zval(nz+1,nb)
	      fz1=w1*fzval(nz,nb)+w2*fzval(nz,nb+1)
	      fz2=w1*fzval(nz+1,nb)+w2*fzval(nz+1,nb+1)
	    endif
	  enddo
	endif

	logfz=(fz1 + (z-z1)*(fz2-fz1)/(z2-z1))
	fz_table=logfz/fznorm(nb)

	if (fz_table.lt.1.e-10) fz_table=0

	return
	end


	subroutine fz_table_read

	implicit none

C Table values
	integer*4 nzmax,nbmax
	parameter(nbmax=15)
	parameter(nzmax=500)
	integer*4 numb,numz(nbmax)
	real*8 bval(nbmax)
	real*8 zval(nzmax,nbmax),fzval(nzmax,nbmax)
	real*8 fznorm(nbmax)
	common/fztable/bval,fznorm,zval,fzval,numb,numz

	real*8 tmpb,tmpz,tmpfz
	real*8 oldbeta
	real*8 tmpsum
	real*8 z,dz
	integer*4 nb,nz

	real*8 fz_table
	external fz_table

	open(unit=55,file='benhar.data',status='old')

	nb=0
	oldbeta=0
100	read (55,*,end=999) tmpb,tmpz,tmpfz
	if (tmpb.ne.oldbeta) then	!reading next beta
	  oldbeta=tmpb
	  if (nb.ne.0) numz(nb)=nz-1	!keep # of z values for this beta
	  nb=nb+1
	  bval(nb)=tmpb
	  nz=1
	endif
	zval(nz,nb)=tmpz
	fzval(nz,nb)=tmpfz
	if (fzval(nz,nb).gt.1.e-8) then	!increment nz, otherwise overwrite 
	  nz=nz+1
	endif
	goto 100		!read next value

999	continue
	numb=nb
	numz(nb)=nz-1

	do nb=1,numb
	  fznorm(nb)=1		!must be =1 so that call to fz_table works.
	  tmpsum=0
	  do nz=1,1000	!integral from 0 to 5.
	    z=nz/200.
	    dz=1/200.
	    tmpsum=tmpsum+fz_table(z,bval(nb))*dz
	  enddo
	  fznorm(nb)=tmpsum
!	  write(6,*) 'beta=',bval(nb),' had normalization=',fznorm(nb)
	enddo

	return
	end

	real*8 function slac_emc_func(x,A)

	real*8 x,A
	real*8 alpha,C


	alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
	1  -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
	2  +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	slac_emc_func = C*A**alpha

	return 

	end
