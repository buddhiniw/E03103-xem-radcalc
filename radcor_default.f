	program radcor
C+______________________________________________________________________________
!
! DESCRIPTION:
!
! This program calculates radiative corrections to the inelastic continuum
! for SLAC experiment NE3.
!
! AUTHOR:
!
! D. H. Potterveld, Kellogg Lab, Caltech.
!
! MODIFICATION HISTORY:
!
! 1/5/88 - Modified to use new IMSL library routines for spline handling. (DHP)
! 
! 6/23/93 - modified to include ne18 analysis. (DLW)
!
! 1995-1996 - Modified beyond recognition for CEBAF x>1 analysis.
!
C-______________________________________________________________________________

	implicit none
	include 'math_physics.inc'
	include 'include.inc'
	include 'spect.inc'

C HARDWIRED STUFF!  BAD-EVIL-NAUGHTY!!!!!
	real*8 sys_err			!Syst. error applied before radcor,
	parameter (sys_err = .001)	!but removed at end. (accep corr)

C input/output files.
	integer*4 util				!main input and output files.
	integer*4 logfile,logfile2		!I/O channels to files.

C Status variables.
	logical converged
	logical modelhack	!.true. for radcor - tweak 15,23 degree model
				!.false. for tbincor - no tweak

C Input data values.
	real*8 tmpep,tmpcounts,tmpstaterr,tmpsyserr
	real*8 tmpsigfact,tmptb,tmpta,tmpctstosig

	real*8 ep(max_pnts)
	real*8 counts(max_pnts,max_sets)	!#/counts
	real*8 counts_sys(max_pnts,max_sets)	!sys. error
	real*8 counts_stat(max_pnts,max_sets)	!stat. error
	real*8 sigfact(max_pnts,max_sets)	!cts --> sigma conversion for each point.
	real*8 ctstosig(max_sets)		!overall (global) cts-->sig factor for each set.

C Summed values (and temporaries for the summing)
	real*8 tmperr2,sumstat,sumsys
	real*8 tmpxsec,sumxsec
	real*8 sumcounts(max_pnts)		!sum of counts for ep bin.
	real*8 dsumcounts(max_pnts)		!error (for fit) on sumcounts.
	real*8 sumcounts_sys(max_pnts)		!sys error.
	real*8 sumcounts_stat(max_pnts)		!stat error.


C Model data values.
	real*8 mod_rad(max_pnts,max_sets)	!radiated model cross section.
	real*8 mod_counts(max_pnts,max_sets)	!radiated model counts.
	real*8 mod_sumcounts(max_pnts)		!rad counts (sum of data sets).
	real*8 mod_raw(max_pnts)		!raw model cross section.
	real*8 mod_rawcounts(max_pnts,max_sets)	!raw model counts.
	real*8 mod_sumrawcounts(max_pnts)	!raw counts (sum of data sets).

C Output data values.
	real*8 xsecn(max_pnts)			!Radiatively unfolded
	real*8 xsecn_stat(max_pnts)		! cross section
	real*8 xsecn_sys(max_pnts)		! and errors.
	real*8 xsecn_raw(max_pnts)		!Measured xsec (no rad. corr)

C General Kinematics
	real*8 e1,deg
	real*8 thrad,snsq
	real*8 e2,nu,x,q2		!temporary kinematical variables.

C Target Info
	real*8 aux(7)
	real*8 zp,zn,a_nuc,m_nuc,m_atom,teff

C Radiative correction parameters.
	real*8 dee
	real*8 tb(max_sets),ta(max_sets)

C Iterative model corrections.
	integer*4 ier
	integer*4 n_good,ngood		!Data point counters.
	integer*4 iteration,iter_max	!Iteration stuff.
	real*8 corr(max_pnts)
	real*8 dcorr(max_pnts)
	real*8 new_corr(max_pnts)
	real*8 wsq_neg(max_pnts)
	real*8 max_wsq_neg, min_wsq_neg
	real*8 xfit(max_pnts),yfit(max_pnts),cfit(max_pnts,3)
	real*8 wk(7*(max_pnts+2)),dum(max_pnts)
	real*8 rfit(max_pnts),drfit(max_pnts)
	real*8 chisq
	real*8 rtst,chisq_old
	real*8 exp_rad(max_pnts),dexp_rad(max_pnts)
	real*8 rad_corr(max_pnts)

C Extra junk
	real*8 d1,d2			!dummy variables
	integer*4 filelen
 	integer*4 i,j,ind
	integer*4 ip,iset
	integer*4 numset		!# of different data sets.
	character*20 file		!File name.
	character*80 tmpfile
	character*80 line

	common/newmodfix/ xfit,yfit,cfit,ngood
	common/wsq_limits/ min_wsq_neg, max_wsq_neg

C Function declarations.

	real*8 radiate,evalspline

C ============================ Executable Code =================================

C Open a log file.

	util=81
	logfile=82
	logfile2=83

!	type 1001,'$Log file name: '
	read (5,1001,end=991) file
	do ind=20,1,-1
	  if (file(ind:ind).eq.' ') filelen=ind-1
	enddo

C open outputfile.
	tmpfile = file(1:filelen)//'.radout'
	open (unit=logfile,status='unknown',name=tmpfile)
	tmpfile = file(1:filelen)//'.radout2'
	open (unit=logfile2,status='unknown',name=tmpfile)

C Open data file, read in kinematics and target information (A/Z)
	tmpfile = file(1:filelen)//'.preradcor'
	open(unit=util,name=tmpfile,status='old',readonly,err=990)
	read(util,*) e1
	read(util,*) deg
	read(util,*) zp
	read(util,*) zn
	read(util,'(a80)') line

	a_nuc = zp + zn
	thrad = deg*d_r
	snsq = sin(thrad/2.)**2.

C Get target parameters (aux).  Don't need ta,tb,tgt_len,etc... so put in
C dummy lengths.  To get aux, just need to input A,Z.  Need theta<>0 to avoid
C error on cryotargets (for l/sin(theta))

!       call target_info(theta_tar,theta,tb,ta, z,    a,tgt_len,tgt_thick,aux,teff)
	call target_info(       0.,   0.,d1,d2,zp,a_nuc,     1.,       1.,aux,teff)

	if (a_nuc.eq.2) then
	  m_atom=2.0150
	else if (a_nuc.eq.12) then
          m_atom=12.0110
        else if (a_nuc.eq.56) then
          m_atom=55.8470
        else if (a_nuc.eq.197) then
          m_atom=196.9665
	else
	  write(6,*) 'CANT FIGURE OUT M_ATOM.  GET OFF YOUR LAZY ASS AND GET RID OF THIS HARDWIRED CRAP!!'
	endif
        m_nuc = m_atom*m_amu-zp*m_e

C initialize some stuff.
	min_wsq_neg=9999.
	max_wsq_neg=-9999.
	dee = 0.005	   			!For calls to RADIAT.

	ITER_MAX = 5

	do ip = 1,max_pnts
	  do iset = 1,max_sets
	    sigfact(ip,iset) = 0.0	!flag that no data point exists here
	    mod_rawcounts(ip,iset) = 0.0
	    mod_counts(ip,iset) = 0.0
	  enddo
	enddo

C Read in binned data.  Increment ip whenever E' (ep) changes.  Increment
C iset whenever a new data set is detected (by noting a new value of sigfact).

	ip=1
	numset=1
	read(util,*) ep(1),counts(1,1),counts_stat(1,1),counts_sys(1,1),
     &		sigfact(1,1),tb(1),ta(1),ctstosig(1)

	counts_stat(1,1) = counts_stat(1,1)*counts(1,1) !convert from frac err.
	nu = e1-ep(1)
	q2 = 4.*e1*ep(1)*snsq
	wsq_neg(1) =  q2 - m_p*m_p - 2.*m_p*nu	!-1*(Missing mass)^2
	xfit(1) = wsq_neg(1)
	min_wsq_neg = min( min_wsq_neg , wsq_neg(1) )
	max_wsq_neg = max( max_wsq_neg , wsq_neg(1) )

	do while (ip.le.max_pnts)
	  read(util,*,end=4) tmpep,tmpcounts,tmpstaterr,tmpsyserr,
     &		tmpsigfact,tmptb,tmpta,tmpctstosig

* Incremenet ip if new E' and determine which data set the point comes from.

	  if (abs(tmpep-ep(ip)).gt.0.001) ip=ip+1

	  iset = 0				!no set assigned yet.
	  do i = 1,numset			!check existing data sets.
	    if(abs(tmpctstosig-ctstosig(i))/((tmpctstosig+ctstosig(i)/2.))
     &        .lt. 0.0001) iset=i
	  enddo
	  if (iset.eq.0) then		!no set found --> new data set.
	    numset = numset+1		!increment number of sets.
	    iset = numset
	    tb(iset) = tmptb		!read in tb,ta,sigfact for this set.
	    ta(iset) = tmpta
	    ctstosig(iset) = tmpctstosig
	  endif

* Put input values into arrays once ip,iset known.
	  ep(ip) = tmpep
	  counts(ip,iset) = tmpcounts
	  counts_stat(ip,iset) = tmpstaterr*tmpcounts	!CONVERT FROM FRACTIONAL ERROR.
	  counts_sys(ip,iset) = tmpsyserr
	  sigfact(ip,iset) = tmpsigfact	!sigfact(ip,iset)<>0 --> real data.

* Calculated quantities.  
!!!!!!!!!!!!!!!!!!NEED STAT ERROR FROM COMBINE.PCM (FOR e+/dummy sub)!!!!!!!!!
!	  counts_stat(ip,iset) = sqrt(counts(ip,iset)+1)
	  nu = e1-ep(ip)
	  q2 = 4.*e1*ep(ip)*snsq
	  wsq_neg(ip) =  q2 - m_p*m_p - 2.*m_p*nu	!-1*(Missing mass)^2
	  xfit(ip) = wsq_neg(ip)
	  min_wsq_neg = min( min_wsq_neg , wsq_neg(ip) )
	  max_wsq_neg = max( max_wsq_neg , wsq_neg(ip) )
	enddo

4	continue		!reached end of input file
	n_good = ip
	ngood = n_good		!version for new spline routine

! Combine counts and errors for data.!!!!!!!ERRORS NOT DONE CORRECTLY YET!!!!!!
	do ip = 1,n_good
	  sumcounts(ip)=0.
	  dsumcounts(ip)=0.
	  sumstat=0.
	  sumsys=0.
	  sumxsec=0.
	  do iset=1,numset
	    if (sigfact(ip,iset).gt.0) then
	      sumcounts(ip) = sumcounts(ip) + counts(ip,iset)
	      tmperr2 = (counts_stat(ip,iset)/counts(ip,iset))**2
	      sumstat = sumstat + 1./tmperr2
	      sumsys  = sumsys  + counts_sys(ip,iset)/tmperr2

	      tmpxsec = counts(ip,iset)*sigfact(ip,iset)
	      sumxsec = sumxsec + tmpxsec/tmperr2
	    endif
	  enddo
	  sumcounts_stat(ip) = 1/sqrt(sumstat)		!fractional uncertainties.
	  sumcounts_sys(ip) = sumsys/sumstat
	  sumcounts_stat(ip) = sumcounts_stat(ip)*sumcounts(ip)
	  sumcounts_sys(ip) = sumcounts_sys(ip)*sumcounts(ip)
	  dsumcounts(ip) = sqrt( sumcounts_stat(ip)**2 + (sys_err*sumcounts(ip))**2 )
	  xsecn_raw(ip) = sumxsec/sumstat	!ave xsecn. (weighted by counts)
	enddo

C Setup cross section output file.
	close (unit=util)
	tmpfile = file(1:filelen)//'.radcor'
	open (unit=util,name=tmpfile,status='unknown')
	write (util,'(3f10.3,a)') e1,deg,a_nuc,' = E0, Theta, A'

C Initial output to log file (for this target).
	write (logfile,*) e1,' GeV = Incident energy'
	write (logfile,*) deg,' Degrees = scattering angle'
	write (logfile,*) a_nuc,' = A (target)'
	write (logfile,*) zp,' = Z (target)'
!	write (logfile,*) tb,' = Radiator before in radiation lengths'
!	write (logfile,*) ta,' = Radiator after  in radiation lengths'
	write (logfile,*) aux(1),' = RESOL for DEEPSIG smearing (not used)'
	write (logfile,*) aux(2),' = E_SEP in GeV'
	write (logfile,*) aux(3),' = F(0)'
	write (logfile,*) aux(4),' = BigB'
	write (logfile,*) aux(5),' = a'
	write (logfile,*) aux(6),' = b'
	write (logfile,*) aux(7),' = alpha'
	write (logfile,*) zp,' = Number of protons in target nucleus'
	write (logfile,*) zn,' = Average number of neutrons in target nucleus'
	write (logfile,*) dee,' = Dee in GeV'
	write (logfile,*) ' '
	write (logfile,*) ' nu      exp/rad     uncert cts-rawmodel cts-radmodel   cts-data     sig_iter'

C Initialize some things.
	do ind = 1,n_good
	  xsecn(ind) = 0.D00			!Zero unfolded cross section.
	  xsecn_stat(ind) = 0.D00
	  xsecn_sys(ind) = 0.D00
	  corr(ind) = 1.D00			!Assume xsec model is perfect.
	  new_corr(ind) = 1.D00
	  rfit(ind) = 1.D00
	  drfit(ind) = 1.D00
	enddo

	call cubgcv(wsq_neg,rfit,drfit,ngood,yfit,cfit,max_pnts,1.d0,0,dum,wk,ier)
	if (ier.ne.0) then
	  write(6,*) 'ier = ',ier,' for first call to cubgcv'
	  write(6,*) 'ngood=',ngood
	  write(6,*) 'wsq_neg(1:ngood)=',(wsq_neg(j),j=1,ngood)
	  write(6,*) 'rfit(1:ngood)=',(rfit(j),j=1,ngood)
	  write(6,*) 'drfit(1:ngood)=',(drfit(j),j=1,ngood)
	endif

	iteration = 0
	converged = .false.
	modelhack = .true.

C-----------------------------
! Main iterative loop.
C-----------------------------

!	do iteration = 1,iter_max
	do while (.not.converged .and. iteration.le.iter_max)
	  iteration=iteration+1

	  do ip = 1,n_good			!clear counts, calculate unradiated model.
	    mod_sumcounts(ip) = 0.D0
	    mod_sumrawcounts(ip) = 0.D0
	    e2 = ep(ip)
	    call sigmodel(e1,e2,deg,a_nuc,zp,m_nuc,aux,mod_raw(ip),modelhack)
	    call model_fix(wsq_neg(ip),mod_raw(ip))
	  enddo

	  do iset = 1,numset		!for each set, compute counts at each E'
	    call radiate_init(e1,deg,tb(iset),ta(iset),zp,zn,m_nuc,aux,dee)

! If there is a data point at given E' (sigfact>0) then compute unradiated
! model cross section [nb/(MeV-sr)].  Then compute radiated model using
! peaking approximation. The radiative correction at this iteration is sig/rad.
! Cross section = counts * sigfact

	    do ip = 1,n_good		!Loop over all data points.
	      e2 = ep(ip)
	      if (sigfact(ip,iset).gt.0) then
		mod_rawcounts(ip,iset) = mod_raw(ip) / sigfact(ip,iset)
		mod_sumrawcounts(ip) = mod_sumrawcounts(ip) + mod_rawcounts(ip,iset)

		mod_rad(ip,iset) = radiate(e2,modelhack)
		mod_counts(ip,iset) = mod_rad(ip,iset) / sigfact(ip,iset)
		mod_sumcounts(ip)=mod_sumcounts(ip) + mod_counts(ip,iset)
	      endif
	    enddo		  !loop over E' points
	  enddo			!loop over data sets (iset)

! Compare data to model, generate ratio, dratio,  and chisq of difference.
! Adjust model correcting parameters if EXP_RAD ratio significantly far from 1

	  chisq = 0.D00
	  do ip = 1,n_good
	    if (mod_sumcounts(ip).gt.0.) then
	      exp_rad(ip) = sumcounts(ip) / mod_sumcounts(ip)
	      dexp_rad(ip) = dsumcounts(ip) / mod_sumcounts(ip)
	      rad_corr(ip) = mod_sumcounts(ip) / mod_sumrawcounts(ip)

	      rtst = abs(exp_rad(ip)-1.D00)/dexp_rad(ip)
	      if (ip.le.n_good) chisq = chisq + rtst**2
	      corr(ip) = corr(ip)*exp_rad(ip)
	      dcorr(ip) = corr(ip)*dsumcounts(ip)/sumcounts(ip)

! Radiatively unfolded cross section (at this iteration).

	      xsecn(ip) = mod_raw(ip)*exp_rad(ip)
	      xsecn_sys(ip) = xsecn(ip)*sumcounts_sys(ip)/sumcounts(ip)
	      xsecn_stat(ip) = xsecn(ip)*sumcounts_stat(ip)/sumcounts(ip)
	      write(logfile,1002) e1-ep(ip),exp_rad(ip),dexp_rad(ip),mod_sumrawcounts(ip),mod_sumcounts(ip),sumcounts(ip),xsecn(ip)

C Radiated model is ZERO. If this was a 'GOOD' point, remove it from fit
C by adjusting N_GOOD. This also removes all points above it. The assumption
C here is that RAD is zero because we have reached too small of an energy loss.

	    else
	      write (logfile,'(1x,a,1x,f8.5)') 'RAD = zero! for Ep =',ep(ip)
	      n_good = min(ip-1,n_good)
	      ngood = n_good
	    endif

	  enddo		!loop over E' points.


C Check convergence criteria.

	  chisq = chisq/n_good			!Chisq per deg. of freedom.
	  write(6,*) ' chisq=',real(chisq),' chisq_old=',real(chisq_old)

	  if (iteration.eq.1) then
	    converged = .false.
	  elseif (chisq.le.1.0) then
	    converged = .true.
	  elseif ((chisq_old - chisq).lt..4D00 .and. (chisq_old-chisq).gt.0.d00) then
	    converged = .true.
	  elseif (iteration.ge.iter_max) then	!Maximum iterations reached
	    converged = .true.			!without true convergence.
	    write (logfile,*) ' WARNING: MAXIMUM ITERATION REACHED!'
	    type *,' WARNING: MAXIMUM ITERATION REACHED!'
CCC	    if (chisq.gt.chisq_old) converged = .false.
	  endif

C Output results of this iteration.

	  write (logfile,*) ' '
	  write (logfile,1004) iteration,chisq,chisq_old
	  write (logfile,*) ' '

	  if (chisq.gt.chisq_old .and. iteration.ne.1) then		!If model got worse:
	    write(6,*) 'CHISQUARED INCREASED AFTER PREVIOUS ITERATION'
	    stop
	  else
	    chisq_old = chisq			!Else, Save chisq per D.O.F. and restore correction factor
	  endif

C Smooth the new model correction parameters if we are going to iterate.

	  if (.not.converged) then
	    do ip = 1,n_good			!Form R*4 vectors.
	      rfit(ip) = corr(ip)
	      drfit(ip) = 2.*dcorr(ip)
	    enddo

	    call cubgcv(wsq_neg,rfit,drfit,ngood,yfit,cfit,max_pnts,1.d0,0,dum,wk,ier)
	    if (ier.ne.0) then
	      write(6,*) 'ier = ',ier,' for cubgcv call from radcor_cebaf.f'
	      write(6,*) 'ngood=',ngood
	      write(6,*) 'wsq_neg(1:ngood)=',(wsq_neg(j),j=1,ngood)
	      write(6,*) 'rfit(1:ngood)=',(rfit(j),j=1,ngood)
	      write(6,*) 'drfit(1:ngood)=',(drfit(j),j=1,ngood)
	    endif

	    do ip = 1,n_good
	      new_corr(ip) = evalspline(wsq_neg(ip),n_good,xfit,yfit,cfit)
	    enddo
	  endif

C Set correction factors to smoothed values (if iterating) or to the
C unadjusted values at the beginning of this iteration (if not iterating).

	  do ip = 1,n_good
	    corr(ip) = new_corr(ip)
	  enddo

	enddo		!Main iteration loop.

! Banner to log file.
	write (logfile,*)' '
	write (logfile,*)'  e1     e2     q2     x    sigexp     sigcor  dsig_stat dsig_sys   rad_corr'
	write (logfile,*) ' '

! Cross section results to log file.
	do ip = n_good,1,-1
	  e2 = ep(ip)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq
	  x = q2/(2.*m_p*nu)
	  write (logfile,1003) e1,e2,q2,x,xsecn_raw(ip),xsecn(ip),xsecn_stat(ip),xsecn_sys(ip),1./rad_corr(ip)
	enddo

! Second banner to log file.

	write (logfile,*) ' '
	write (logfile,*) '  Eloss        W**2     model_corr'
	write (logfile2,*) '  Eloss        W**2     model_corr'
	write (logfile,*) ' '

! Final model correction parameters.

	do ip = 1,n_good
	  nu = e1 - ep(ip)
	  write (logfile,1005) nu,-wsq_neg(ip),corr(ip)
	  write (logfile2,1005) nu,-wsq_neg(ip),corr(ip)
	enddo

! Results to cross section file. Cross section units are
! nB/(Mev-ster). Only 'good' data points are output.

	do ip = 1,n_good
	  e2 = ep(ip)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq
	  x = q2/(2.*m_p*nu)
	  if (x.lt.a_nuc) then
	    write (util,1006) e2,xsecn(ip),xsecn_stat(ip),xsecn_sys(ip)
	  endif
	enddo
	close (unit=util)


C End of file jump destinations.

992	write (logfile,*) ' '
	write (logfile,*) ' =========== End Of Analysis =========='
	close (unit=logfile)
	stop ' '

990	type *, ' RADCOR_NE18: Error reading file = ',tmpfile
991	stop 'Analysis aborted!'

C =========================== Format Statements ================================

1001	format(a)
1002	format(1x,f6.3,x,f9.3,x,f9.3,x,e12.4,x,e12.4,x,e12.4,x,e12.4)
1003	format(f6.3,3f7.3,4e10.3,f10.4)
1004	format(' iteration = ',i3,x,' chisqr = ',e10.4,x,' oldchisqr = ',e10.4)
1005	format(f9.6,2x,f10.4,2x,1f10.4)
1006	format(f6.3,3e12.4)
	end
