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
	real*8 pwid			!Expected momentum binsize (15 MeV)
	parameter (pwid = .015)
	real*8 e_err,ep_err,theta_err
	parameter (e_err=0.0003)	!0.03% uncertainty in ebeam.
	parameter (ep_err=0.0003)	!0.03% uncertainty in p.
	parameter (theta_err=0.0005*r_d) !0.5mr uncertainty in theta.

 	integer*4 i,j,ind		!Integer temporaries.
	integer*4 n_pts,n_good		!Data point counters.
	integer*4 ngood
	integer*4 iteration,iter_max	!Iteration stuff.
	integer*4 logfile,logfile2,logfile3,util,kine,accfile,logfile4 !I/O channels to files.
	integer*4 modelfile

	logical converged
	logical goodmodel	!is there a model xsecn. file for this run.
	logical foundsig	!is the point in the model file.
	logical enddata		!stop reading data (after first bin w/o counts)
	logical modelhack	!.true. for radcor - tweak 15,23 degree model
				!.false. for tbincor - no tweak

C Acceptance function.
	integer*4 mcnum		!mc events generated.
	integer*4 mcnumd,mcnumt	!number of delta/theta bins in accep fcn.

	real*8 ctsperbin	!#/events generated per delta-theta bin.
	real*8 mcdmin,mcdmax	!mc delta populated range.
	real*8 mctmin,mctmax	!mc theta populated range.
	real*8 mcpmin,mcpmax	!mc   phi populatedrange.
	real*8 mctlen		!mc target length
	real*8 mcdlo,mcdhi,mcdbin	!mc delta lowval,hival,binsize
	real*8 mctlo,mcthi,mctbin	!mc theta lowval,hival,binsize
	real*8 acceptance(150,100)	!acceptance function.
	real*8 mcpval(150)		!momentum values for bin left-edges.
	real*8 tmpd

 	common/accepblock/
     >     mcdlo, mcdhi, mcdbin,
     >     mctlo, mcthi, mctbin,
     >     mcpval, acceptance, mcnumd, mcnumt

	integer*4 mcnumtlo,mcnumthi,mcnumtgood

C Data variables.
	real*8 xibin,cnts,err,pbin,dbin	!xi,counts,error,p,delta
	real*8 xibinlo,xibinhi
	real*8 cts_to_sig,domega
	real*8 rxsecn(max_pnts)		!Raw (experimental)
	real*8 drxsecn(max_pnts)	!cross section & error.
	real*8 xsecn(max_pnts)		!Radiatively unfolded
	real*8 dxsecn(max_pnts)		!cross section & error.
	real*8 rad(max_pnts)		!radiated model xsection

	real*8 aux(7)

	real*8 xi_width

	integer*4 filelen
	character*20 file		!File name.
	character*80 tmpfile
	character*80 line,tmpline

	real*8 wsq_neg(max_pnts),wsqneg
	real*8 max_wsq_neg, min_wsq_neg
	real*8 rfit(max_pnts)
	real*8 drfit(max_pnts)

	real*8 xfit(max_pnts),yfit(max_pnts),cfit(max_pnts,3)
	real*8 wk(7*(max_pnts+2)),dum(max_pnts)
	integer ier

	real*8 corr(max_pnts)
	real*8 dcorr(max_pnts)
	real*8 new_corr(max_pnts)
	real*8 ep(max_pnts),xi(max_pnts)

	real*8 tmpaccep,accep_max
	real*8 rawx_cent   !for calculating binning corrections
	real*8 rawx2_cent,rawx2_sum
	real*8 waccep_sum
	real*8 wrawx_sum
	real*8 wxsec_sum	    !weighted by accep and domega(=sin(theta))
	real*8 domega_sum	    !sum of domega (=sin(theta))
	real*8 binaccep,bincorr
	real*8 tbincorr(max_pnts),pbincorr(max_pnts)
	real*8 tbincorr2(max_pnts)
	real*8 degbin

	real*8 t1,t2,t3		!Temporaries
	real*8 junk		!Very temporary

	real*8 e1,ep_cent,deg,theta_tar
	real*8 tgt_len,tgt_thick,corr_factor,zp,a_nuc,m_atom,q_tot,eta
	real*8 thrad,snsq,xi2,e2,fit_cut
	real*8 tb,ta,zn,m_nuc
	real*8 dee,ep_elas,nu,nu_elas,q2,x,q3sq
	real*8 chisq,difsum
	real*8 sig,rad_corr,exp_rad,dexp_rad
	real*8 rtst,chisq_old
	real*8 syserr(max_pnts)	!fractional systematic uncertanties from accep,binning,...
	real*8 xsecnorm         !sig=sig*xsecnorm from inefficiencies and such.
	real*8 sig0,dsig

	integer*4 maxnumpts
	parameter (maxnumpts=3000)
	integer*4 numpts

	real*8 ximod(maxnumpts)
	real*8 degmod(maxnumpts)
	real*8 sigmod(maxnumpts)
	real*8 modelsig

	integer*4 runno

	common/newmodfix/ xfit,yfit,cfit,ngood
	common/wsq_limits/ min_wsq_neg, max_wsq_neg

C Function declarations.

	real*8 radiate,evalspline
	real*8 evalp,evalxi,dxidp
	real*8 evalaccep

C ============================ Executable Code =================================

C Open a log file.

	accfile=84
	logfile2=83
	logfile3=85
	logfile4=87
	modelfile=86
	logfile=82
	util=81
	kine=80

	modelhack=.false.

	type 1001,'$Log file name: '
	read (5,1001,end=991) file
	do ind=20,1,-1
	  if (file(ind:ind).eq.' ') filelen=ind-1
	enddo

C Get setup info.
!	FIT_CUT = .05
!	FIT_CUT = .10		!cut 100 MeV from elastic peak of NUCLEUS.

	write(6,*) 'USING 10 MeV CUTOFF FOR RADCOR!!!!!!!!!!'
	write(6,*) 'USING 10 MeV CUTOFF FOR RADCOR!!!!!!!!!!'
	write(6,*) 'USING 10 MeV CUTOFF FOR RADCOR!!!!!!!!!!'
	write(6,*) 'USING 10 MeV CUTOFF FOR RADCOR!!!!!!!!!!'
	write(6,*) 'USING 10 MeV CUTOFF FOR RADCOR!!!!!!!!!!'

	FIT_CUT = .01		!cut 10 MeV from elastic peak of NUCLEUS.
	ITER_MAX = 12

	type 1001,'$Working...'

	read (file(1:4),'(i4)') runno
	if (file(5:5).eq.'a') then
	  write (6,*) 'run= ',file(1:5)
	else
	  write (6,*) 'run= ',file(1:4)
	endif

C Read in kinematics from .kine file.

	tmpfile = 'kine/'//file(1:filelen)//'.kine'

	call getkine(runno,e1,ep_cent,deg,theta_tar,tgt_len,tgt_thick,
     >		     corr_factor,zp,a_nuc,m_atom,q_tot,tmpfile)

* CONVERT NEGATIVE TGT THICKNESS (FLAG FOR DUMMY RUNS) BACK TO +.

	tgt_thick=abs(tgt_thick)

	thrad = deg*d_r
	snsq = sin(thrad/2.)**2.
	zn = a_nuc - zp
	m_nuc = m_atom*m_amu-zp*m_e
	eta = tgt_thick/m_atom*6.022d23/abs(cosd(theta_tar))

C Get tb,ta - effective radiation lengths (before/after interaction).
	call target_info(theta_tar,deg,tb,ta,zp,a_nuc,tgt_len,tgt_thick,aux)
	dee = 0.005	   			!For calls to RADIAT.

C Read in acceptance function.
	if (zp.gt.3) then
	  write(tmpfile,'(''~/cebaf/cebaf2/johna/NTUPLE_HMS/acceptance.solid.'',i2)')
     >							nint(deg)
	else if (tgt_len.le.6.) then
	  write(tmpfile,'(''~/cebaf/cebaf2/johna/NTUPLE_HMS/acceptance.short.'',i2)')
     >							nint(deg)
	else
	  write(tmpfile,'(''~/cebaf/cebaf2/johna/NTUPLE_HMS/acceptance.long.'',i2)')
     >							nint(deg)
	endif

	open(unit=accfile,file=tmpfile,status='old')

C read in acceptance function, set good theta limits (mcnumtlo,hi,good)
	read(accfile,*) mcnum		!mc events generated.
	read(accfile,*) mcdmin,mcdmax	!mc delta populated range.
	read(accfile,*) mctmin,mctmax	!mc theta populated range.
	read(accfile,*) mcpmin,mcpmax	!mc   phi populatedrange.
	read(accfile,*) mctlen		!mc target length
	read(accfile,*) mcdlo,mcnumd,mcdbin	!mc delta lowval,numbins,binsize
	read(accfile,*) mctlo,mcnumt,mctbin	!mc theta lowval,numbins,binsize

!Recon.f subtracts 0.4mr from yptar, which adds 0.4mr to theta.  Add same
! here, since it's a physics offset, not a data-vs-MC offset.
	MCTLO=MCTLO+0.0004

	mcdhi=mcdlo+mcnumd*mcdbin
	mcthi=mctlo+mcnumt*mctbin
	ctsperbin=mcnum * mcdbin/(mcdmax-mcdmin) * mctbin/(mctmax-mctmin)

	mcnumtlo=999			!lowest bin with Accep<>0.
	mcnumthi=0			!highst bin with Accep<>0.
	do j=1,mcnumt		!read in acceptance function
	  degbin = r_d * (mctlo + mctbin*(j-0.5))
	  do i=1,mcnumd
	    read(accfile,*) acceptance(i,j)
	    acceptance(i,j)=acceptance(i,j)/ctsperbin
	    if (acceptance(i,j).gt.0.001) then
	      mcnumtlo=min(mcnumtlo,j)
	      mcnumthi=max(mcnumthi,j)
	    endif
	  enddo
	enddo
	mcnumtgood=mcnumthi-mcnumtlo+1

C get acceptance at central delta bin.
	waccep_sum=0
	domega_sum=0
	i=int(mcnumd/2)
	do j=mcnumtlo,mcnumthi
	  degbin = r_d * (mctlo + mctbin*(j-0.5))
	  waccep_sum=waccep_sum+acceptance(i,j)*sind(degbin)
	  domega_sum=domega_sum+sind(degbin)
	enddo
	accep_max=waccep_sum/domega_sum	!accep. for central delta bin.
	write(6,*) 'mcnumtlo,mcnumthi,accep_max=',mcnumtlo,mcnumthi,accep_max

	do i=1,mcnumd+1		!generate vector of momentum lower edges.
	  tmpd = mcdlo + (i-1)*mcdbin
	  mcpval(i) = ep_cent*(1+tmpd/100.)
	enddo

C initialize some stuff.
	i = 1
	min_wsq_neg=9999.
	max_wsq_neg=-9999.
	ep_elas = e1/(1.+2.*e1*snsq/m_nuc)		!E-prime of elastic pk.
	enddata=.false.

C open outputfile.
	tmpfile = file(1:filelen)//'.out'
	open (unit=logfile,status='unknown',name=tmpfile)
	tmpfile = file(1:filelen)//'.out2'
	open (unit=logfile2,status='unknown',name=tmpfile)
	tmpfile = file(1:filelen)//'.modelout'
	open (unit=logfile3,status='unknown',name=tmpfile)
	tmpfile = file(1:filelen)//'.log'
	open (unit=logfile4,status='unknown',name=tmpfile)

C read in model xsecn. at fixed xi,theta values if they exist in a file.
	goodmodel=.false.
	tmpfile = file(1:filelen)//'.model'
	open (unit=modelfile,status='old',name=tmpfile,err=980)

	goodmodel=.true.
	do ind=1,maxnumpts
	  read(modelfile,*,end=980) ximod(ind),degmod(ind),sigmod(ind)
          numpts=ind
        enddo

980	continue		!no model
	write(6,*) 'GOODMODEL=',goodmodel
	write(logfile4,*) 'GOODMODEL=',goodmodel

C Open data file, read in solid angle=(phi region in MC)*
C sum of (sin(theta)*dtheta) for the theta bins in recon.f
	tmpfile = 'xi/'//file(1:filelen)//'.xi_raw'
	open(unit=util,name=tmpfile,status='old',readonly,err=990)
	read(util,*) domega

	domega=0
	do j=mcnumtlo,mcnumthi
	  degbin = r_d * (mctlo + mctbin*(j-0.5))
	  domega=domega+sind(degbin)*mctbin*(mcpmax-mcpmin)
	enddo

C Read in binned data.  File contains xi,counts,error.  Convert xi to
C E' (GeV).  Discard points with nu less than FIT_CUT from elastic peak.
C Note that this will crash for zero counts, so some cut on counts or 
C dcounts/counts is required.  Require at least one count, and as soon
C as we hit a bin with no counts, stop reading in file (mask all later
C points)

	do while (i.le.max_pnts)
44	  read(util,*,end=4) xibin,cnts,err
	  if (enddata) cnts=0
	  pbin=evalp(xibin,deg,e1)
	  dbin=100.*(pbin/ep_cent-1)
	  if (dbin.ge.-5 .and. cnts.eq.0) enddata=.true.
	  xibinlo=evalxi((pbin-pwid/2.),deg,e1)
	  xibinhi=evalxi((pbin+pwid/2.),deg,e1)
	  if (dbin.lt.mcdlo .or. dbin.gt.mcdhi) goto 44
!	  IF (DBIN.LT.-6.5) GOTO 44 	!CUTTING BELOW -6.5% FOR NOW.
!	  IF (DBIN.GT.11.5) GOTO 44 	!CUTTING ABOVE 11.5% FOR NOW.
	  IF (DBIN.LT.-9.0) GOTO 44 	!CUTTING BELOW -9.0% FOR NOW.
	  IF (DBIN.GT.11.5) GOTO 44 	!CUTTING ABOVE 11.5% FOR NOW.

!!	  IF (DBIN.GT.3 .AND. DBIN.LT.13) THEN
!!	    CNTS=CNTS*(1-(4-4/5*ABS(DBIN-8))/100)	!4% AT DEL=8, DOWN TO 0% AT 3%,13%
!!	    ERR=ERR*(1-(4-4/5*ABS(DBIN-8))/100)
!!	  ENDIF

	  nu=e1-pbin

	  if (abs( nu - (pwid*nint(nu/pwid)) ).gt.0.0001) then
	    write(6,*) 'NU is NOT a multiple of',pwid*1000.,'MeV'
	    write(6,*) 'xi,p,nu=',xibin,pbin,nu
	    stop
	  endif

	  t1 = pbin		!E'
	  t2 = cnts		!sig
	  t3 = err		!dsig
	  t3 = sqrt( t3**2 + (t2*sys_err)**2 )

C convert counts per xi bin to cross section (dsigma/dxi/domega)
	  xi_width=xibinhi-xibinlo
	  cts_to_sig=corr_factor/(q_tot/1.602d-19)/eta/xi_width/domega
	  cts_to_sig=cts_to_sig*1.d+33		! cm^2/sr --> nb/sr
C convert from dsigma/dxi/dOmega to dsigma/dE'/dOmega
	  cts_to_sig=cts_to_sig*dxidp(t1,deg,e1)

	  if ( abs(xi_width/dxidp(t1,deg,e1)/1000.-pwid)/pwid .gt. 0.01) then
	    write(6,*) 'xi bin=',xi_width
	    write(6,*) 'corresponds to p bin of',xi_width/dxidp(t1,deg,e1)
	    stop
	  endif


C get acceptance for xi bin.

	  waccep_sum = 0.
	  domega_sum = 0.

	  do j = mcnumtlo,mcnumthi
	    degbin = r_d * (mctlo + mctbin*(j-0.5))
	    waccep_sum = waccep_sum + 
     >         evalaccep(xibin,xibinlo,xibinhi,degbin,e1,ep_cent)*sind(degbin)
	    domega_sum = domega_sum + sind(degbin)
	  enddo

C Keep points with (Ep < Ep_elastic) and dsig/sig<50%.
C Also require that acceptance of xi bin be > about 50% by requireing
C that the unweighted AND weighted acceptances are >50% of the peak accep.

	  if (t1.ge.(ep_elas-fit_cut)) goto 4		!break input loop
	  if (t2.ge.1 .and. (waccep_sum/domega_sum)/accep_max.gt.0.5) then
	    ep(i) = t1
	    xi(i) = xibin
	    nu_elas = ep_elas - t1
	    rxsecn(i) = t2*cts_to_sig
	    drxsecn(i) = t3*cts_to_sig
	    q2 = 4.*e1*t1*snsq
	    wsq_neg(i) =  q2 - m_p*m_p - 2.*m_p*nu	!-1*(Missing mass)^2
	    xfit(i) = wsq_neg(i)
	    min_wsq_neg=min(min_wsq_neg,wsq_neg(i))
	    max_wsq_neg=max(max_wsq_neg,wsq_neg(i))
	    if (nu_elas.ge.fit_cut) n_good = i
	    i = i + 1
	  endif
	enddo

C Check that data is binned appropriatly
	if (i.gt.1) then
	  if (abs( (t1-ep(i-1)) - pwid ).gt.0.0001) then !not 15 MeV step.
	    write(6,*) 'Input data not binned in ',pwid*1000.,'MeV bins'
	    stop
	  endif
	endif


4	continue		!reached end of input file
	n_pts = i - 1
	ngood = n_good		!version for new spline routine
	write(6,*) 'ngood,n_good,n_pts=',ngood,n_good,n_pts
	write(logfile4,*) 'ngood,n_good,n_pts=',ngood,n_good,n_pts


C Setup cross section output file.
	close (unit=util)
	tmpfile = file(1:filelen)//'.radcor'
	open (unit=util,name=tmpfile,status='unknown')
	write (util,'(3f10.3,a)') e1,deg,a_nuc,' = E0, Theta, A'

C Initial output to log file (for this target).
	write (logfile,*) ' '
	write (logfile,*) 'CEBAF RADIATIVE CORRECTIONS'
	write (logfile,*) ' '
	write (logfile,*) e1,' GeV = Incident energy'
	write (logfile,*) deg,' Degrees = scattering angle'
	write (logfile,*) ' '
	write (logfile,*) 'A,Z (target) = ',a_nuc,zp
	write (logfile,*) ' '
	write (logfile,*) tb,' = Radiator before in radiation lengths'
	write (logfile,*) ta,' = Radiator after  in radiation lengths'
	write (logfile,*) aux(1),' = RESOL for DEEPSIG smearing (not used)'
	write (logfile,*) aux(2),' = E_SEP in GeV'
	write (logfile,*) aux(3),' = F(0)'
	write (logfile,*) aux(4),' = BigB'
	write (logfile,*) aux(5),' = a'
	write (logfile,*) aux(6),' = b'
	write (logfile,*) aux(7),' = alpha'
	write (logfile,*) zp,' = Number of protons in target nucleus'
	write (logfile,*) zn,' = Average number of neutrons in target nucleus'
	write (logfile,*) m_nuc,' = Mass of target nucleus in GeV/c^2'
	write (logfile,*) dee,' = Dee in GeV'
	write (logfile,*) fit_cut,' = Fit cut-off (GeV from elastic peak)'
	write (logfile,*) ' '
	write (logfile,*) ' '
	write (logfile,*) '   Eloss    exp/rad      uncert   modelxsecn  radiated  data'
	write (logfile,*) ' '

C Initialize some things.
	do i = 1,n_pts
	  xsecn(i) = 0.D00			!Zero unfolded cross section.
	  dxsecn(i) = 0.D00
	  corr(i) = 1.D00			!Assume xsec model is perfect.
	  new_corr(i) = 1.D00
	  rfit(i) = 1.D00
	  drfit(i) = 1.D00
	enddo

	call cubgcv(wsq_neg,rfit,drfit,ngood,yfit,cfit,max_pnts,1.d0,0,dum,wk,ier)
	if (ier.ne.0) then
	  write(6,*) 'ier = ',ier,' for first call to cubgcv'
	  write(6,*) 'ngood=',ngood
	  write(6,*) 'wsq_neg(1:ngood)=',(wsq_neg(j),j=1,ngood)
	  write(6,*) 'rfit(1:ngood)=',(rfit(j),j=1,ngood)
	  write(6,*) 'drfit(1:ngood)=',(drfit(j),j=1,ngood)
	endif




C Theta binning correction (apply BEFORE radiative corrections).
C For each xi bin, take mcnumt theta bins.
C Generate dsigma/dE'/dOmega (from model) for each xi,theta point.
C Convert to dsigma/dxi/dOmega, and take ratio of sigma(xi,thetacent) to
C integral of sigma(xi,thetabin)*accep(xi,thetabin) over the theta range
C (as is measured in the experiment).

	do i=1,n_pts

	  wxsec_sum = 0.
	  wrawx_sum = 0.
	  domega_sum = 0.
	  waccep_sum = 0.
	  syserr(i) = 0.	!FRACTIONAL uncertainty.

	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Accep = 1.0% & 4%*(1-A)
	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Bincent = 1% & .1*corr
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !PID Effic.
	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Q measurement
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !Track eff.
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !Tracking cuts.
	  syserr(i) = sqrt(syserr(i)**2 + 0.025**2) !Radcorr = 2.5% 

	  if (abs(tgt_thick-0.2129).le.0.001) then	!Fe_2%
	    syserr(i) = sqrt(syserr(i)**2 + 0.010**2)
	  else if (abs(tgt_thick-0.8034).le.0.001) then	!Fe_6%
	    syserr(i) = sqrt(syserr(i)**2 + 0.020**2)
	  else if (abs(tgt_thick-0.3768).le.0.001) then	!Au_6%
	    syserr(i) = sqrt(syserr(i)**2 + 0.010**2)
	  else if (abs(tgt_thick-0.8915).le.0.001) then	!C_2%
	    syserr(i) = sqrt(syserr(i)**2 + 0.005**2)
	  else if (abs(tgt_thick-2.5096).le.0.001) then	!C_6%
	    syserr(i) = sqrt(syserr(i)**2 + 0.005**2)
	  else
	    syserr(i) = sqrt(syserr(i)**2 + 0.010**2)
	    write(6,*) 'DID NOT FIND TGT_THICK IN LIST OF SOLID TARGETS, using 1% syserr'
!	    write(6,*) 'THIS HAD BETTER BE A CRYOTARGET RUN'
	  endif

! consitant within systematics ==> don't increase error for it.
!	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Hydrogen Norm.

          xsecnorm = 1.000              !elastic normalization
          xsecnorm = xsecnorm / .998    !Ave Cerenkov effic.  ??????
          xsecnorm = xsecnorm / 1.00	!Ave Calorim effic.

	  e2 = ep(i)
	  rxsecn(i) = rxsecn(i)*xsecnorm
	  drxsecn(i) = drxsecn(i)*xsecnorm
	  xibin = xi(i)
	  pbin = evalp(xibin,deg,e1)
	  xibinlo = evalxi((pbin-pwid/2.),deg,e1)
	  xibinhi = evalxi((pbin+pwid/2.),deg,e1)

C get sigma at central angle
	  degbin=deg
	  call radiate_init(e1,degbin,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig = radiate(pbin,modelhack)				!Radiated model
*	  call sigmodel(e1,pbin,degbin,a_nuc,zp,m_nuc,aux,sig,modelhack) !Raw model

	  sig = sig / dxidp(pbin,degbin,e1)	!dsig/dE' -> dsig/dxi

	  rawx_cent = sig
	  do j = mcnumtlo,mcnumthi
	    degbin = r_d * (mctlo + mctbin*(j-0.5))
	    pbin = evalp(xibin,degbin,e1)

	    foundsig=.false.
	    if (goodmodel) then
	      do ind=1,numpts
		if ( abs(ximod(ind)-xibin).lt.0.0001 .and. 
     >		     abs(degmod(ind)-degbin).lt.0.0001 ) then
	          modelsig = sigmod(ind)
	          foundsig=.true.
		endif
	      enddo
	    endif

	    if (.not.foundsig .or. abs(degbin-deg).le.(r_d*mctbin)) then
	      call radiate_init(e1,degbin,tb,ta,zp,zn,m_nuc,aux,dee)
	      sig = radiate(pbin,modelhack)				!Radiated model
*	      call sigmodel(e1,pbin,degbin,a_nuc,zp,m_nuc,aux,sig,modelhack) !Raw model

	      sig = sig / dxidp(evalp(xibin,degbin,e1),degbin,e1)	!dsig/dE' -> dsig/dxi
	    endif


	    if (foundsig) then
	      if (abs(degbin-deg).le.(r_d*mctbin)) then
	        if (abs(sig-modelsig)/((sig+modelsig)/2).gt.0.001) then ! >.1% off
		  write(6,*) 'BAD MODEL VALUE!!!!!!!!!!'
		  write(6,*) 'BAD MODEL VALUE!!!!!!!!!!'
		  write(6,*) 'BAD MODEL VALUE!!!!!!!!!!'
		  write(6,*) 'BAD MODEL VALUE!!!!!!!!!!'
		  write(6,*) 'xi,deg,sig,modelsig=',xibin,degbin,sig,modelsig
		  stop
	        endif
	      else
		sig=modelsig
	      endif
	    endif

	    tmpaccep = evalaccep(xibin,xibinlo,xibinhi,degbin,e1,ep_cent)
	    wxsec_sum = wxsec_sum + sig*tmpaccep*sind(degbin)
	    wrawx_sum = wrawx_sum + sig*sind(degbin)
	    domega_sum = domega_sum + sind(degbin)
	    waccep_sum = waccep_sum + tmpaccep*sind(degbin)

	    write(logfile3,'(2f10.5,e15.6)') xibin,degbin,sig

	  enddo

! accepcorr is APPROXIMATLY the acceptance correction.  xsecn_cent/xsecn_ave
! combines BOTH the CROSS SECTION WEIGHTED acceptance and the theta binning
! correction.  accepcorr divides out the UNWEIGHTED acceptance, as a rough
! approximation.  The product of the two is still xsecn_cent/xsecn_ave.

	  if (wxsec_sum.gt.0.0) then
	    bincorr = min(1000. , (rawx_cent/(wxsec_sum/domega_sum)) )
	    binaccep = wxsec_sum/wrawx_sum
	    tbincorr(i)=min(1000. , bincorr*binaccep)	!totalcorr/accepcorr
	    tbincorr2(i)=rawx_cent/(wxsec_sum/(waccep_sum+.0000001)+.0000001)
	  else
	    bincorr=1000.
	    binaccep=0.
	    tbincorr(i)=0.
	    tbincorr2(i)=0.
	  endif

C xi binning correction (apply BEFORE radiative corrections).
C For each xi bin, take 40 xi bins.  We already took the acceptance
C into account in the theta binning correction, so we just take the
C cross section variation into account.

C WRONG!!!!  Should also take acceptance VARIATIONS over xi into account.
C THIS IS NOT DONE YET, but should be a small effect.  Applying the
C acceptance variation over such small steps just introduces the single-bin
C statistics of the acceptance MC.

	  rawx2_sum=0
	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  do j = 0,20
	    xi2 = xibinlo + j/20.*(xibinhi-xibinlo)
	    e2 = evalp(xi2,deg,e1)

	    foundsig=.false.
	    if (goodmodel) then
	      do ind=1,numpts
		if ( abs(ximod(ind)-xi2).lt.0.0001 .and. 
     >		     abs(degmod(ind)-deg).lt.0.0001 ) then
	          modelsig = sigmod(ind)
	          foundsig=.true.
		endif
	      enddo
	    endif

!!!	    if (.not.foundsig .or. abs(e2-ep(i)).le.(0.0001)) then
	    if (.not.foundsig) then
	      call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	      sig = radiate(e2,modelhack)				!Radiated model
*	      call sigmodel(e1,e2,deg,a_nuc,zp,m_nuc,aux,sig,modelhack)	!Raw Model

	      sig = sig / dxidp(e2,deg,e1)	!dsig/dE' -> dsig/dxi
	    endif


	    if (foundsig) then
!!	      if (abs(e2-ep(i)).le.(0.0001)) then
!!	        if (abs(sig-modelsig)/((sig+modelsig)/2).gt.0.001) then ! >.1% off
!!		  write(6,*) 'xi,deg,sig,modelsig=',xibin,degbin,sig,modelsig
!!	        endif
!!	      else
		sig=modelsig
!!	      endif
	    endif

	    rawx2_sum=rawx2_sum+sig
	    if (j.eq.10) rawx2_cent=sig

	    write(logfile3,'(2f10.5,e15.6)') xi2,deg,sig

	  enddo

	  if (rawx2_cent.ne.0) then
	    pbincorr(i)=rawx2_cent/(rawx2_sum/21.)
	  else
	    pbincorr(i)=999.
	  endif

C assign uncertainty to acceptance correction.  Base 2% offset (already
C added in) + .04*(1-A). Gives 2.5% uncertainty at A=.5 (sqrt(1.5**2+2.0**2).
	  junk = (waccep_sum/domega_sum)/accep_max !approx. accep/accep_max
	  syserr(i)=sqrt( syserr(i)**2 + (0.04*(1-junk))**2 )

C assign 10%*(1-tbincorr2) uncertainty to tbincorr.  Max. correction is 20%
C except at 15 degress, extremely lo nu (nu <= 105 MeV).
	  syserr(i)=sqrt( syserr(i)**2 + (0.10*(1-tbincorr2(i)))**2 )

	  rxsecn(i)=rxsecn(i)*bincorr*pbincorr(i)
	  drxsecn(i)=drxsecn(i)*bincorr*pbincorr(i)

C Find error in cross section due to uncertainties in e,ep,theta.

!Raw model:
!	  call sigmodel(e1,e2,deg,a_nuc,zp,m_nuc,aux,sig0,modelhack)
!
!	  call sigmodel(e1,e2*(1+ep_err),deg,a_nuc,zp,m_nuc,aux,sig,modelhack)
!	  dsig = abs( sig0 - sig )
!	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from ep offset
!	write(6,'(a,f6.2,a)') 'dsig/sig (ep)=',nint(10000*dsig/sig)/100.,'%'
!
!	  call sigmodel(e1*(1+e_err),e2,deg,a_nuc,zp,m_nuc,aux,sig,modelhack)
!	  dsig = abs( sig0 - sig )
!	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from e offset
!	write(6,'(a,f6.2,a)') 'dsig/sig (e)=',nint(10000*dsig/sig)/100.,'%'
!
!	  call sigmodel(e1,e2,deg+theta_err,a_nuc,zp,m_nuc,aux,sig,modelhack)
!	  dsig = abs( sig0 - sig )
!	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from theta offset
!	write(6,'(a,f6.2,a)') 'dsig/sig (th)=',nint(10000*dsig/sig)/100.,'%'

!	Radiated model.

	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig = radiate(e2,modelhack)			!Radiated model

	  dsig = abs( sig - radiate(e2*(1+ep_err),modelhack) )
	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from ep offset
	write(6,'(a,f6.2,a)') 'dsig/sig (ep)=',nint(10000*dsig/sig)/100.,'%'

	  call radiate_init(e1*(1+e_err),deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  dsig = abs( sig - radiate(e2,modelhack) )
	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from e offset
	write(6,'(a,f6.2,a)') 'dsig/sig (e)=',nint(10000*dsig/sig)/100.,'%'

	  call radiate_init(e1,deg+theta_err,tb,ta,zp,zn,m_nuc,aux,dee)
	  dsig = abs( sig - radiate(e2,modelhack) )
	  syserr(i)=sqrt(syserr(i)**2 + (dsig/sig)**2)	!error from theta offset
	write(6,'(a,f6.2,a)') 'dsig/sig (th)=',nint(10000*dsig/sig)/100.,'%'

	  write(6,'(a,8f8.4)') 'xi,rawA,weiA,tc,tc2,pc,totc,sys=',xibin,
     >		(waccep_sum/domega_sum)/accep_max,	! accepcorr(nowei)
     >		binaccep/accep_max,			! accepcorr(wei)
     >		tbincorr(i),				! tbincorr
     >		tbincorr2(i),				! tbincorr
     >		pbincorr(i),				! pbincorr
     >		1/(bincorr*pbincorr(i)*accep_max),	! totcorr
     >		syserr(i)
	  write(logfile4,'(a,8f8.4)') 'xi,rawA,weiA,tc,tc2,pc,totc,sys=',xibin,
     >		(waccep_sum/domega_sum)/accep_max,	! accepcorr(nowei)
     >		binaccep/accep_max,			! accepcorr(wei)
     >		tbincorr(i),				! tbincorr
     >		tbincorr2(i),				! tbincorr
     >		pbincorr(i),				! pbincorr
     >		1/(bincorr*pbincorr(i)*accep_max),	! totcorr
     >		syserr(i)

	enddo

	iteration = 0
	converged = .false.
	modelhack = .true.
	write(6,*) '----------beginning-correction-procedure--------------------'

C-----------------------------
! Main iterative loop.
C-----------------------------

	do while (.not.converged)
	iteration = iteration + 1
	chisq = 0.D00
	difsum = 0.D00
	call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)

	do i = 1,n_pts				!Loop over all data points.

! Compute unradiated model cross section. Units are nb/(MeV-ster).
	  e2 = ep(i)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq

	  call sigmodel(e1,e2,deg,a_nuc,zp,m_nuc,aux,sig,modelhack)

	  call model_fix(wsq_neg(i),sig)

! Compute radiated model using peaking approximation. The radiative correction
! for this data point at this iteration is sig/rad.

	  rad(i) = radiate(e2,modelhack)		!Radiated model.
	  if (rad(i).gt.0.) then
	    rad_corr = sig/rad(i)			!Radiative correction.
	    exp_rad = rxsecn(i)/rad(i)			!Ratio of exp/rad
	    dexp_rad = drxsecn(i)/rad(i)		!and its uncert.

! Compute sum (over data points) of square of difference between radiative
! correction at this iteration and previous iteration.

	    difsum = difsum + (rad_corr - xsecn(i)/rxsecn(i))**2

! Adjust model correcting parameters if EXP_RAD ratio significantly
! differs from unity.

	    rtst = abs(exp_rad-1.D00)/dexp_rad
	    corr(i) = corr(i)*exp_rad
	    dcorr(i) = corr(i)*drxsecn(i)/rxsecn(i)

! Radiatively unfolded cross section (at this iteration).

	    xsecn(i) = sig*exp_rad
	    dxsecn(i) = xsecn(i)*sqrt((drxsecn(i)/rxsecn(i))**2-sys_err**2)

	    write(logfile,1002) nu,exp_rad,dexp_rad,sig,rad(i),rxsecn(i)

	  else

C Radiated model is ZERO. If this was a 'GOOD' point, remove it from fit
C by adjusting N_GOOD. This also removes all points above it. The assumption
C here is that RAD is zero because we have reached too small of an energy
C loss.

	    write (logfile,'(1x,f8.5,a,3x)') nu, 'RAD = zero!'
	    corr(i) = corr(i-1)
	    dcorr(i) = dcorr(i-1)
	    xsecn(i) = 0.D00
	    dxsecn(i) = 0.D00
	    n_good = min(i-1,n_good)
	    ngood = n_good
	  endif

! Compute a 'chi-squared' for deviation of "exp_rad's" from unity. Only
! Include data points which were not flagged as too close to incident energy.

	  if (i.le.n_good) chisq = chisq + ((exp_rad-1.D00)/dexp_rad)**2
	enddo

C Check convergence criteria.

	chisq = chisq/n_good			!Chisq per deg. of freedom.
	write(6,*) ' chisq=',real(chisq),' chisq_old=',real(chisq_old)
	write(logfile4,*) ' chisq=',real(chisq),' chisq_old=',real(chisq_old)


	if (iteration.eq.1) then
	  converged = .false.
	  difsum = 0.D00
	elseif (chisq.le.1.0) then
	  converged = .true.
	elseif ((chisq_old - chisq).lt..4D00 .and. (chisq_old-chisq).gt.0.d00) then
	  converged = .true.
	elseif (difsum.lt.n_pts/1000000.D00) then
	  converged = .true.
	elseif (iteration.ge.iter_max) then	!Maximum iterations reached
	  converged = .true.			!without true convergence.
	  write (logfile,*) ' WARNING: MAXIMUM ITERATION REACHED!'
	  type *,' WARNING: MAXIMUM ITERATION REACHED!'
CCC	  if (chisq.gt.chisq_old) converged = .false.
	endif

C Output results of this iteration.

	write (logfile,*) ' '
	write (logfile,1004) iteration,chisq,chisq_old,difsum/n_pts
	write (logfile,*) ' '

	if (chisq.gt.chisq_old .and. iteration.ne.1) then		!If model got worse:
	  write(6,*) 'CHISQUARED INCREASED AFTER PREVIOUS ITERATION'
	else
	  chisq_old = chisq			!Else, Save chisq per D.O.F. and restore correction factor
	endif

C Smooth the new model correction parameters if we are going to iterate.

	if (.not.converged) then
	  do i = 1,n_good			!Form R*4 vectors.
	    rfit(i) = corr(i)
	    drfit(i) = dcorr(i)
	  enddo

	  call cubgcv(wsq_neg,rfit,drfit,ngood,yfit,cfit,max_pnts,1.d0,0,dum,wk,ier)
	  if (ier.ne.0) then
	    write(6,*) 'ier = ',ier,' for cubgcv call from radcor_cebaf.f'
	    write(6,*) 'ngood=',ngood
	    write(6,*) 'wsq_neg(1:ngood)=',(wsq_neg(j),j=1,ngood)
	    write(6,*) 'rfit(1:ngood)=',(rfit(j),j=1,ngood)
	    write(6,*) 'drfit(1:ngood)=',(drfit(j),j=1,ngood)
	  endif

	  do i = 1,n_good
	    new_corr(i) = evalspline(wsq_neg(i),n_good,xfit,yfit,cfit)
	  enddo
	endif

C Set correction factors to smoothed values (if iterating) or to the
C unadjusted values at the beginning of this iteration (if not iterating).

	do i = 1,n_good
	  corr(i) = new_corr(i)
	enddo

	enddo

C------------------------------
!End of main iteration loop.
C------------------------------

C Final output.

! Banner to log file.

	write (logfile,*)' '
	write (logfile,*)'   e1     e2     nu     q2     x     '//
     >	'sigexp     sigcor     dsigcor ratio(expcor/exp)'
	write (logfile,*) ' '

! Cross section results to log file.

	do i = n_good,1,-1
	  e2 = ep(i)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq
	  x = q2/(2.*m_p*nu)
	  q3sq = nu**2. + q2
	  exp_rad = xsecn(i)/rxsecn(i)
	  write (logfile,1003) e1,e2,nu,q2,x,rxsecn(i),xsecn(i),
     >	             dxsecn(i),exp_rad
	enddo

! Second banner to log file.

	write (logfile,*) ' '
	write (logfile,*) '  Eloss        W**2     model_corr tbincorr pbincorr'
	write (logfile2,*) '  Eloss        W**2     model_corr tbincorr pbincorr'
	write (logfile,*) ' '

! Final model correction parameters.

	do i = 1,n_good
	  nu = e1 - ep(i)
	  write (logfile,1005) nu,-wsq_neg(i),corr(i),tbincorr2(i),pbincorr(i)
	  write (logfile2,1005) nu,-wsq_neg(i),corr(i),tbincorr2(i),pbincorr(i)
	enddo

! Results to cross section file. Cross section units are
! nB/(Mev-ster). Only 'good' data points are output.

	do i = 1,n_good
	  e2 = ep(i)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq
	  x = q2/(2.*m_p*nu)
	  if (x.lt.a_nuc) then
	    write (util,1006) e2,xsecn(i),dxsecn(i),syserr(i)*xsecn(i),rad(i),rxsecn(i),tbincorr(i)
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
1002	format(1x,f8.5,2x,f9.3,2x,f9.3,2x,e12.4,2x,e12.4,2x,e12.4)
1003	format(5f7.3,3e11.3,f10.4)
1004	format(' iteration = ',i3,x,' chisqr = ',e10.4,x,' oldchisqr = ',e10.4,x,' difsum = ',e10.4)
1005	format(f9.6,2x,f10.4,2x,3f10.4)
1006	format(f6.3,5e12.4,f7.3)
1007	format(1x,f8.5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4)
	end


C evalaccep generates acceptance(xi,theta) by interpolating between
C acceptance(delta,theta) bins.  Delta is %, theta is degrees.

	function evalaccep(xi,xilo,xihi,theta,e1,e2)

	implicit none

	include 'math_physics.inc'

	real*8 xi,theta,e1,e2	!xi,theta,beam energy,central momentum.

	integer*4 mcnumd,mcnumt	!number of delta/theta bins in accep fcn.
	real*8 mcdlo,mcdhi,mcdbin	!mc delta lowval,hival,binsize
	real*8 mctlo,mcthi,mctbin	!mc theta lowval,hival,binsize
	real*8 acceptance(150,100)	!acceptance function.
	real*8 mcpval(150)		!momentum values for bin left-edges.

	common/accepblock/
     >     mcdlo, mcdhi, mcdbin,
     >     mctlo, mcthi, mctbin,
     >     mcpval, acceptance, mcnumd, mcnumt

	integer*4 id,it

	real*8 thr
	real*8 xilo,xihi	!central,lo,hi xi values of XI BIN.
	real*8 pcent,plo,phi	!central,lo,hi p values of XI BIN.
	real*8 accepsum		!for taking weighted ave. of bins.
	real*8 lowp,highp	!tmp variables used to find frac.
	real*8 frac		!fraction of delta bin containing xi bin.
	real*8 evalaccep
	real*8 totfrac
	real*8 dcent,del

	real*8 evalp,evalxi

	thr=theta*d_r

	pcent=evalp(xi,theta,e1)
	plo=evalp(xilo,theta,e1)
	phi=evalp(xihi,theta,e1)
	totfrac=(phi-plo)/(e2*mcdbin/100.)	!size of xi bin (in p bins)

C apply delta cut to acceptance identical to values in recon.f
	plo=max(plo,mcpval(1))			!delta > -14%
	phi=max(phi,mcpval(1))

	phi=min(phi,mcpval(mcnumd+1))		!delta < +14%
	plo=min(plo,mcpval(mcnumd+1))

C For delta bins that contain part of the xi bin, take frac*accep for
C numerator, and totfrac for denominator.  frac=fraction of delta bin in xi bin.

	accepsum=0
	it = nint( (thr-mctlo)/mctbin - 0.5 ) + 1	!theta bin index
	do id=1,mcnumd
	  if (plo.lt.mcpval(id+1) .and. phi.gt.mcpval(id)) then !xibin in pbin
	    lowp=max(plo,mcpval(id))
	    highp=min(phi,mcpval(id+1))
	    frac=(highp-lowp)/(mcpval(id+1)-mcpval(id))
	    accepsum=accepsum+frac*acceptance(id,it)
	  endif
	enddo

	evalaccep=0
	if (totfrac.gt.0) evalaccep=accepsum/totfrac

	return
	end
