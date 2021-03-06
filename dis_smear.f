 	subroutine dis_smear(e_in,ep_in,theta,sig,anuc,znuc,smearflag)

C       This routine calculates the convolution of F2_N(x/z,qsq) with f_A(z,beta)
C       following the procedure of Benhar, et.al. (from Jun 1997 preprint).
C       
C       F2_A(x,qsq)=Integral[from x to A] of f_A(z,beta),F2_N(x/z,qsq)
C       
C       F2_N is the nucleon structure function for a nucleon with momentum z
C       and a quark momentum fraction x/z.  Beta is defined as abs(q3)/nu.
C       Here, we choose beta as a function of scattering angle only, such
C       that the fit matches the data, and expected ratio of DIS/QE from other
C       calculations.
C       
C       The f(z,beta) is taken from tables (generated by Benhar), for beta values
C       of 1.0,1.1,...,1.5.  They are calculated for Nuclear Matter (?).
C       
C       F2_N(x/z,qsq) from bodel fits to proton and neutron.  Use N=Z=A/2 to
C       get symmetric nucleus with A nucleons.  Correct for neutron excess and
C       EMC effect with fit to EMC data.
C       
C       If smearflag=.false., return Z*f2_p + (A-Z)*f2_n - no smearing, no EMC.
C       (Probably should have EMC, but I don't plan on using unsmeared).
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
	real*8 r, dr, iso_cor
	real*8 sumef,rntd2,rautd2,emc_corr,zz
	real*8 x1,x2,frac,f2dconv,sigm,f2nf2p,f2deut,f2dconv_tmp
	real*8 corfac,betafac, emc_corr2, f2dslac
        real*8 pval
	logical smearflag
	logical initialized /.false./
	logical goodfit /.true./

	integer*4 ind,j
	integer*4 intpts
	parameter (intpts=150)

c       These are the coefficients needed to correct for neutron excess
c       and the emc effect. they were obtained by fitting the ratios of
c       cross sections [cs(a>2)/cs(d2)] obtained in e139.the correcting
c       function is a 6-th order polynomial in x.
c       right now (3.6.6) ef() is not used to correct for neutron excess!

	real*8 ef(7) /-0.00136693,-0.00510425,-0.0375986,-0.0946004,
	1    -0.122435,-0.0112751,0.406435/
	real*8 ecor(9) /163.2827,-2686.153,19152.23,-76860.48,189906.6,
	1    -295822.3,283716.9,-153174.6,35640.7/

	real*8 fz_table
	external fz_table

	real*8 slac_emc_func
	external slac_emc_func
	real*8 xem_emc_func
	external xem_emc_func
        real*8 xem_emc_func2
	external xem_emc_func2



c----------------------this is for debugging only--------------
c             open(unit=7,file='ajiborn.out',status='unknown')
c             open(unit=8,file='ajirad1.out',status='unknown')
c              open(unit=9,file='ajirad2.out',status='unknown')
c               open(unit=10,file='ajirad3.out',status='unknown')
c-----------------------------------------------------------------
c=======================================================================================
c       first get the deuterium f2 and later build nuc targets from that

	   if (.not.initialized) then
	      call fz_table_read
	      initialized=.true.
	   endif

	   s_in=sind(theta/2)**2
	   nu_in=e_in-ep_in
	   qsq_in=4*e_in*ep_in*s_in
	   x_in=qsq_in/2/m_p/nu_in
	   q3_in=sqrt(qsq_in+nu_in**2)

      	    qsq=qsq_in	

	   x1=0.8		!x<x1 --> use xsechd (for deut.)
	   x2=0.9		!x1<x<x2 --> smooth transition from f2d to conv calc.


	   if (x_in.lt.x2) then	!we need slac fit upto x=0.9 (in transition region also)
	      call xsechd(e_in,ep_in,s_in,1,+2,f2dslac,w2,w1,xsectp,sigm,x_test)
        
	   endif

	   
	   if (x_in.lt.x1) then
	      f2deut = f2dslac

            else if (x_in.ge.x1) then ! now start doing convolution calculation

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c       introduce a q2 dependent beta in the convolution integral
c       betafac= -3.31573 + (2.42234*qsq_in)
c       >     + (-0.518855*qsq_in**2)+ (0.0381449*qsq_in**3)
c       betafac= -0.320771 + (0.189242*qsq_in)

	      betafac=  -0.288897 + (0.177758*qsq_in)

	      if (betafac.lt.0.0) then
		 betafac=0.0
	      endif
 
c              betafac=0.5           
	      beta=q3_in/(nu_in+betafac)







c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				!(intpts)-point integration of f(z)*F2(x/z).  f(z) normalized to one, so no
				! explicit check of normalization.  Truncate integral at z=5 (since f(z)=1e-10)
C       amax=min(5.,anuc)
		 amax = 2.0	! aji since everything now depends only on ld2

		 if (x_in.ge.2.) then !truncated integral will give zero
		    sig=0.            !DIS contribution zero beyond x=anuc
		    return
		 endif

		 f2dconv_tmp=0.
		 dz=(amax-x_in)/dble(intpts)

		 do ind=1,intpts
				!calculate modified kinematics to get F2(x/z,qsq).  Fix x_q=x/z, qsq, and theta
		    z=x_in+(ind-0.5)*(amax-x_in)/dble(intpts) !nucleon has x=z
		    xq=x_in/z	!quark has x_q=x/z
		    nu=qsq/2/m_p/xq !from x=qsq/2/m/nu
		    e=(nu+sqrt(nu**2+qsq/s_in))/2 !from qsq=2*e*(e-ep)*s
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
		    
c       nump = ZNUC
c       numn = ANUC-ZNUC
		    
		    nump = 1.0
		    numn = 1.0
		    

		    fz=fz_table(1.10*(z-1)+1,beta) ! now we depend on f2d for all targets, so 10% reduction
		    f2dconv_tmp=f2dconv_tmp + fz * (nump*f2p+numn*f2n) * dz
		    if (z.gt.1 .and. fz.lt.1.e-10) goto 200 !break loop, no more integrand
		 enddo
	      
 200	      continue
	      f2dconv =f2dconv_tmp

c  ------------debug stuff------------------------------------ 
c              if (pval.lt.0) then
c                  write(7,1010) x_in,qsq_in,betafac,beta,f2dconv,sig ! this is born
c              
c              elseif (pval.eq.10) then
c                  write(8,1010) x_in,qsq_in,betafac,beta,f2dconv,sig ! this is radiated1
c              elseif (pval.eq.20) then
c                  write(9,1010) x_in,qsq_in,betafac,beta,f2dconv,sig ! this is radiated2
c               elseif (pval.eq.30) then
c                  write(10,1010) x_in,qsq_in,betafac,beta,f2dconv,sig ! this is radiated3
c              
c               endif
c 1010	format(f6.3,x,f6.3,x,f12.5,x,f12.5,x,e12.7,x,e12.7)
c--------------------------------------------------------------




	      if ((x_in.ge.x1).and.(x_in.lt.x2)) then ! now check whether in transition region 0.8<x<0.9
		 frac = (x_in-x1)/(x2-x1)
		 f2deut= f2dconv*frac + f2dslac*(1.-frac) 
	      elseif(x_in.ge.x2) then
		 f2deut= f2dconv 
	      endif
	   endif

c==========================================================================================
c       now check nuc or not

	if(ANUC.gt.2) then
	   emc_corr2 = xem_emc_func2(x_in,ANUC)
c           emc_corr2=slac_emc_func(x_in,ANUC)
c           emc_corr2=1.
	   f2a=f2deut*emc_corr2*(ANUC/2.) !for A>2 apply the emc correction, and this inelastic emc fit is from xem data
	elseif(ANUC.eq.2) then
	   f2a=f2deut
        else
	   stop			! we dont do hydrogen (yet!!)
 	endif


	sigmott=(19732./(2.0*137.0388*e_in*s_in))**2*(1-s_in)/1.d6
   
C       here's where the R = sigma_L/sigma_T goes!
C       0 = R_xem = 0.32/Q^2
C       1 = R_1990
C       2 = R_1998

C	call RLT(0,e_in,ep_in,theta,r)
C	call RLT(1,e_in,ep_in,theta,r)
C	call RLT(2,e_in,ep_in,theta,r)
C	r=0.32/qsq
	call R1990(X_IN,QSQ,R,DR,GOODFIT)
	w2=f2a/nu_in
	w1=(1.0+nu_in**2/qsq)/(nu_in*(1.0+r))*f2a
	sig=(w2+2.0*tand(theta/2)**2*w1)*sigmott

c       now apply my adhoc fit for different regions 
        call  polynomcor(x_in,corfac)
	sig = sig*corfac
	return
	end

c-----------------------------------------------------------------------------
	function fz_table(z,beta)
	implicit none

C       Table values
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

	if (z.le.zval(1,nb)) then !linear extrapolation of LOG(f(z)).
	   z1=zval(1,nb)
	   z2=zval(2,nb)
	   fz1=w1*fzval(1,nb)+w2*fzval(1,nb+1)
	   fz2=w1*fzval(2,nb)+w2*fzval(2,nb+1)

	else if (z.gt.zval(numz(nb),nb)) then !linear extrapolation of LOG(f(z))
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

c-----------------------------------------------------------------------------
	subroutine fz_table_read
	implicit none

C       Table values
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

!	open(unit=55,file='benhar.data',status='old')
        open(unit=55,file='/group/hallc_ana/xem/pass3/replay/radcor/working/benhar.data',status='old')

	nb=0
	oldbeta=0
 100	read (55,*,end=999) tmpb,tmpz,tmpfz
	if (tmpb.ne.oldbeta) then !reading next beta
	   oldbeta=tmpb
	   if (nb.ne.0) numz(nb)=nz-1 !keep # of z values for this beta
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
	   do nz=1,1000		!integral from 0 to 5.
	      z=nz/200.
	      dz=1/200.
	      tmpsum=tmpsum+fz_table(z,bval(nb))*dz
	   enddo
	   fznorm(nb)=tmpsum
				!	  write(6,*) 'beta=',bval(nb),' had normalization=',fznorm(nb)
	enddo

	return
	end
c-----------------------------------------------------------------------------
	real*8 function slac_emc_func(x,A)
	real*8 x,A,atemp
	real*8 alpha,C

	atemp = A
				!	if(A.eq.4.or.A.eq.3) then  ! emc effect is more like C for these 2...
				!	   atemp = 12
				!	endif

	alpha = -0.070+2.189*x - 24.667*x**2 + 145.291*x**3
	1    -497.237*x**4 + 1013.129*x**5 - 1208.393*x**6
	2    +775.767*x**7 - 205.872*x**8

	C = exp( 0.017 + 0.018*log(x) + 0.005*log(x)**2)
	
	slac_emc_func = C*atemp**alpha
	return 
	end
c-----------------------------------------------------------------------------

	real*8 function xem_emc_func(x,A) ! now compute the emc effect from our own fits.
	real*8 emc

	if(A.eq.3) then

	   emc = -0.726773 
	1	+16.5100*x
	2	-59.7528e*x**2
	3	+96.5573e*x**3
	4	-55.6012e*x**4
	5	-18.9869e*x**5
	6	+23.8024e*x**6

	else if(A.eq.4) then
	   emc =  -3.59542e+00 
	1	+5.31080e+01*x   
	2	-2.43262e+02*x**2
	3	+5.64535e+02*x**3
	4	-6.99660e+02*x**4
	5	+4.35500e+02*x**5
	6	-1.04688e+02*x**6

	else if(A.eq.12) then
	   emc = 2.80660e+00
	1	-2.43874e+01*x   
	2	+1.39945e+02*x**2
	3	-4.21566e+02*x**3
	4	+6.91573e+02*x**4
	5	-5.85967e+02*x**5
	6	+2.00842e+02*x**6


	else
	   emc = slac_emc_func(x,A)
	endif
	
	xem_emc_func = emc

	return
	end

C---------------------------------------------------------------
c       aji
c       changing beta -> beta(Q2) in convolution calculation for nuclear targets
c       After having a reasonable beta(Q2)  LD2 model,
c       we develop nuclear model like this 
c       take the inelastic non isoscaler correted emc ratios of A>2, do a polynomial it and multiply this to f2A
C       (TURNOFF all the existing emc corr and isoscalar corr)
c       
c       so simply f2a = f2d*xem_emc_func2 
c       


	real*8 function xem_emc_func2(x,A) ! now compute the emc effect from our own fits.
	real*8 emc2

	if(A.eq.3) then

           if (x.lt.1.02) then
	      
	      emc2 = 3.60574
	1	   -17.9930*x
	2	   +38.5132*x**2
	3	   -7.30356*x**3
	4	   -71.4791*x**4
	5	   +83.6995*x**5
	6	   -27.5481*x**6

	   else
	      emc2 = 1.55
	   endif

	else if(A.eq.4) then
	   if (x.lt.1.02) then

	      emc2 = 4.0905
	1	   -20.5864*x
	2	   +41.1811*x**2
	3	   -4.62533*x**3
	4	   -74.9766*x**4
	5	   +79.3210*x**5
	6	   -22.8134*x**6


	   else
	      emc2=1.5
	   endif


	else if(A.eq.9) then
	   if (x.lt.1.02) then
	      
	      emc2 = 4.39669
	1	   -21.8919*x
	2	   +41.8461*x**2
	3	   -2.69928*x**3
	4	   -75.3117*x**4
	5	   +75.1341*x**5
	6	   -19.8517*x**6


	   else
	      emc2=1.2
	   endif


	else if(A.eq.12) then
	   if (x.lt.1.01) then
	      
	      emc2 = 4.18521
	1	   -20.8696*x
	2	   +40.8226*x**2
	3	   -3.07714*x**3
	4	   -74.6169*x**4
	5	   +74.8572*x**5
	6	   -19.5719*x**6
	      
	   else
	      emc2=1.4
	   endif



	else if(A.eq.63) then 
	   if (x.lt.1.01) then
	      
	      emc2 = 4.01755
	1	   -19.8687*x
	2	   +39.0684*x**2
	3	   -4.74309*x**3
	4	   -71.9124*x**4
	5	   +78.8172*x**5
	6	   -23.8254*x**6

	   else
	      emc2=1.7
	   endif


	else if(A.eq.197) then

	   if (x.lt.1.01) then
	      
	      emc2 = 4.20323
	1	   -21.0955*x
	2	   +40.3644*x**2
	3	   -2.65750*x**3
	4	   -73.9456*x**4
	5	   +75.4186*x**5
	6	   -20.7732*x**6


	   else
	      emc2=1.4
	   endif
	else
	   emc2 = 1.0
	endif
	xem_emc_func2 = emc2
	return
	end

c-----------------------------------------------------------------------------

	subroutine polynomcor(x,corfac)

	real*8 x

	if (x.lt.0.7) then
	   corfac=1.
      	else if((x.ge.0.7).and.(x.lt.1.1)) then
	   corfac = 19.5848-(47.7807*x)-(0.207693*x**2)+(49.3179*x**3)+ 
	1	(31.2776*x**4)-(31.7735*x**5)-(71.185*x**6)+(51.6815*x**7) 
	else if((x.ge.1.1).and.(x.le.1.2)) then
	   corfac = 1.5
	else if (x.gt.1.2) then
	   corfac = 1.3
	endif

	return
	end
