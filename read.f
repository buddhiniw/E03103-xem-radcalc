	program read

! TEST version with cut of theta<5mr (roughly, based on nearest bin of
! acceptance function).

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
! 1999 - Converted into two programs.  read.f reads in the raw data, 
!    determines the systematics, acceptance and bin centering corrections,
!    etc...  Everything necessary to convert the counts to cross section,
!    except for the radiative corrections.  This way, seperate runs can
!    be combined BEFORE the radiative corrections.
C-______________________________________________________________________________

	implicit none
	include 'math_physics.inc'
	include 'include.inc'
	include 'spect.inc'

C HARDWIRED STUFF!  BAD-EVIL-NAUGHTY!!!!!
	real*8 pwid			!Expected momentum binsize (15 MeV)
	parameter (pwid = .015)
	real*8 e_err,ep_err,theta_err
	parameter (e_err=0.0005)	!0.05% uncertainty in ebeam.
	parameter (ep_err=0.0003)	!0.03% uncertainty in p.
	parameter (theta_err=0.0005*r_d) !0.5mr uncertainty in theta.
	real*8 dsig_e,dsig_ep,dsig_theta

 	integer*4 i,j,ind			!Integer temporaries.
	integer*4 n_good		!Data point counters.
	integer*4 ngood			!For (destructive) fitting routine

	integer*4 util				!main input and output files.
	integer*4 kine,accfile			!kinematics, acceptance files.
	integer*4 logfile			!I/O channels to files.
	integer*4 modelinfile,modeloutfile      !input/output files for sigmodel
        
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
	real*8 domega,acccut
	real*8 cts_to_sig		!convert counts to dsigma/dxi/domega
	real*8 sigfact(max_pnts)	!counver counts to dsigma/dEp/domega
	real*8 numcts(max_pnts)

	real*8 aux(7)

	real*8 xi_width

	integer*4 filelen
	character*20 file		!File name.
	character*80 tmpfile

	real*8 wsq_neg(max_pnts)
	real*8 rfit(max_pnts),drfit(max_pnts)
	real*8 xfit(max_pnts),yfit(max_pnts),cfit(max_pnts,3)
	real*8 wk(7*(max_pnts+2)),dum(max_pnts)
	integer ier

	real*8 ep(max_pnts),xi(max_pnts)

	real*8 tmpaccep,accep_max
	real*8 rawx_cent   !for calculating binning corrections
	real*8 rawx2_cent,rawx2_sum
	real*8 waccep_sum
	real*8 wrawx_sum
	real*8 wxsec_sum	    !weighted by accep and domega(=sin(theta))
	real*8 domega_sum	    !sum of domega (=sin(theta))
	real*8 binaccep(max_pnts),bincorr(max_pnts)
	real*8 tbincorr(max_pnts),pbincorr(max_pnts)
	real*8 tbincorr2(max_pnts)
	real*8 degbin

!stuff for coulomb corrections
	real*8 coulomb_de,coulomb_pt
	real*8 coulomb_theta_e,coulomb_theta_ep,coulomb_theta
	real*8 coulcor_e(max_pnts),coulcor_th(max_pnts)

!stuff for charge-symmetric background subtraction (positron subtraction)
	real*8 teff,posxsec
	real*8 fracerr
	real*8 num_pos,dnum_pos
	real*8 lumin,tmpdomegade
	real*8 dnumstat,dnumsys

	real*8 junk		!Very temporary

	real*8 e1,ep_cent,deg,theta_tar
	real*8 tgt_len,tgt_thick,corr_factor,zp,a_nuc,m_atom,q_tot,eta
	real*8 thrad,snsq,xi2,e2,fit_cut
	real*8 tb,ta,zn,m_nuc
	real*8 dee,ep_elas,nu,nu_elas,q2,x
	real*8 sig0,sig1
	real*8 syserr(max_pnts)	!fractional systematic uncertanties from accep,binning,...
	real*8 staterr(max_pnts) !fractional statistical uncertainty

	integer*4 maxmodpts
	parameter (maxmodpts=3000)
	integer*4 modpts
	real*8 emod
	real*8 pmod(maxmodpts)
	real*8 degmod(maxmodpts)
	real*8 sigmod(maxmodpts)
	common/modelblock/pmod,degmod,sigmod,emod,modpts

	integer*4 runno

	common/newmodfix/ xfit,yfit,cfit,ngood

C Function declarations.

	real*8 radiate
	real*8 evalp,evalxi,dxidp
	real*8 evalaccep

C ============================ Executable Code =================================

C Open a log file.
	kine=80
	util=81
	logfile=82
	accfile=84
	modeloutfile=85
	modelinfile=86

	modelhack=.false.

	type *,'$Log file name: '
	read (5,'(a)',end=991) file
	do ind=20,1,-1
	  if (file(ind:ind).eq.' ') filelen=ind-1
	enddo

C Get setup info.
!	FIT_CUT = .05
!	FIT_CUT = .10		!cut 100 MeV from elastic peak of NUCLEUS.
	FIT_CUT = .02		!cut 20 MeV from elastic peak of NUCLEUS.

	write(6,*) 'USING 20 MeV CUTOFF FROM NUCLEAR ELASTIC!!!!!!!!!!'
	write(6,*) 'USING 20 MeV CUTOFF FROM NUCLEAR ELASTIC!!!!!!!!!!'

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

! Coulomb energy shift, max. p_t kick (scaled from thesis values by (Z-1)/Z.
	if (abs(a_nuc-2.).le.0.0001) then
	  coulomb_de = 0.
	  coulomb_pt = 0.
	else if (abs(a_nuc-12.).le.0.001) then
	  coulomb_de = 0.0032*(zp-1.)/zp
	  coulomb_pt = 0.0027*(zp-1.)/zp
	else if (abs(a_nuc-56.).le.0.001) then
	  coulomb_de = 0.0093*(zp-1.)/zp
	  coulomb_pt = 0.0080*(zp-1.)/zp
	else if (abs(a_nuc-197.).le.0.001) then
	  coulomb_de = 0.0198*(zp-1.)/zp
	  coulomb_pt = 0.0172*(zp-1.)/zp
	else
	  stop 'Bad value of a_nuc - Cant determine coulomb boost.'
	endif

C Get tb,ta - effective radiation lengths (before/after interaction).
	call target_info(theta_tar,deg,tb,ta,zp,a_nuc,tgt_len,tgt_thick,aux,teff)
	dee = 0.005	   			!For calls to RADIAT.

C Read in acceptance function.
	if (zp.gt.3) then
	  write(tmpfile,'(''mc_accept/acceptance.solid.'',i2)') nint(deg)
	else if (tgt_len.le.6.) then
	  write(tmpfile,'(''mc_accept/acceptance.short.'',i2)') nint(deg)
	else
	  write(tmpfile,'(''mc_accept/acceptance.long.'',i2)') nint(deg)
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
!	write(6,*) 'mctlo,mctbin=',mctlo,mctbin
	do j=1,mcnumt			!read in acceptance function
	  degbin = r_d * ((mctlo-0.0004) + mctbin*(j-0.5))
!	write(6,*) 'degbin=',degbin
	  do i=1,mcnumd
	    read(accfile,*) acceptance(i,j)
	    acceptance(i,j)=acceptance(i,j)/ctsperbin
	    if (acceptance(i,j).gt.0.001) then
!!!!! apply a cut of roughly 10mr (up to nearest bin inside of 10mr).
	      if (abs(degbin-deg).lt.(0.002*r_d)) then
	        mcnumtlo=min(mcnumtlo,j)
	        mcnumthi=max(mcnumthi,j)
	      endif
	    endif
	  enddo
	enddo
	mcnumtgood=mcnumthi-mcnumtlo+1
	write(6,*) 'lo/hi/dtheta=',mcnumtlo,mcnumthi,mcnumtgood*mctbin

C get acceptance at central delta bin.
	waccep_sum=0
	domega_sum=0
	i=int(mcnumd/2)
	do ind=mcnumtlo,mcnumthi
	  degbin = r_d * (mctlo + mctbin*(ind-0.5))
	  waccep_sum=waccep_sum+acceptance(i,ind)*sind(degbin)
	  domega_sum=domega_sum+sind(degbin)
	enddo
	accep_max=waccep_sum/domega_sum	!accep. for central delta bin.
	write(6,*) 'mcnumtlo,mcnumthi,accep_max=',mcnumtlo,mcnumthi,accep_max

	do ind=1,mcnumd+1	!generate vector of momentum lower edges.
	  tmpd = mcdlo + (ind-1)*mcdbin
	  mcpval(ind) = ep_cent*(1+tmpd/100.)
	enddo

C initialize some stuff.
	i = 1
	ep_elas = e1/(1.+2.*e1*snsq/m_nuc)		!E-prime of elastic pk.

C open outputfile.
	tmpfile = file(1:filelen)//'.out'
	open (unit=logfile,status='unknown',name=tmpfile)
	tmpfile = file(1:filelen)//'.modelout'
	open (unit=modeloutfile,status='unknown',name=tmpfile)

C read in model xsecn. at fixed xi,theta values if they exist in a file.
	tmpfile = file(1:filelen)//'.model'
	open (unit=modelinfile,status='old',name=tmpfile,err=980)

	modpts=0
	read(modelinfile,*,end=980) emod
	do ind=1,maxmodpts
	  read(modelinfile,*,end=980) pmod(ind),degmod(ind),sigmod(ind)
	  modpts=ind
	enddo
980	continue		!no model

	if (modpts.gt.0 .and. abs(emod-e1).gt.0.0001) then
	  write(6,*) 'e_model=',emod
	  write(6,*) 'e_beam =',e1
	  stop 'BAD!!!!!!!!'
	endif
	write(modeloutfile,*) e1

	write(6,*) 'GOODMODEL=',(modpts.gt.0)
	write(logfile,*) 'GOODMODEL=',(modpts.gt.0)

C Open data file, read in solid angle=(phi region in MC)*
C sum of (sin(theta)*dtheta) for the theta bins in recon.f
	tmpfile = 'xi/'//file(1:filelen)//'.xi_raw'
	open(unit=util,name=tmpfile,status='old',readonly,err=990)
	read(util,*) domega

	domega=0
	do ind=mcnumtlo,mcnumthi
	  degbin = r_d * (mctlo + mctbin*(ind-0.5))
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
	  pbin = evalp(xibin,deg,e1)
	  nu = e1 - pbin
	  dbin = 100.*(pbin/ep_cent-1)
	  xibinlo = evalxi((pbin-pwid/2.),deg,e1)
	  xibinhi = evalxi((pbin+pwid/2.),deg,e1)
	  xi_width = xibinhi - xibinlo

! throw out point if xi=100, dbin too high or low, and all points after N=0
	  if (xibin.gt.99)   goto 44	!set to 100 for bins with E'>E_beam
	  if (xibinhi.gt.99) goto 44	!evalxi=100 if E'>E_beam
!!!	  if (dbin.lt.-9.0)  goto 44 	!cutting below -9.0% for now.
!!!	  if (dbin.gt.11.5)  goto 44 	!cutting above 11.5% for now.
	  if (dbin.lt.-11.0)  goto 44 	!cutting below -9.0% for now.
	  if (dbin.gt.13.5)  goto 44 	!cutting above 11.5% for now.
	  if (dbin.lt.mcdlo .or. dbin.gt.mcdhi) goto 44
	  if (dbin.ge.-5 .and. cnts.eq.0) goto 4      !stop at first n=0 point.

! check for proper (15 MeV) binning
	  if (abs( nu - (pwid*nint(nu/pwid)) ).gt.0.0001) then
	    write(6,*) 'NU is NOT a multiple of',pwid*1000.,'MeV'
	    write(6,*) 'xi,p,nu=',xibin,pbin,nu
	    stop
	  endif
	  if ( abs(xi_width/dxidp(pbin,deg,e1)/1000.-pwid)/pwid .gt. 0.01) then
	    write(6,*) 'xi bin=',xi_width
	    write(6,*) 'corresponds to p bin of',xi_width/dxidp(pbin,deg,e1)
	    stop
	  endif

! convert counts per xi bin to cross section (dsigma/dxi/domega)
	  cts_to_sig = corr_factor/(q_tot/1.602d-19)/eta/xi_width/domega
	  cts_to_sig = cts_to_sig*1.d+33	! cm^2/sr --> nb/sr
          cts_to_sig = cts_to_sig / 1.000	!elastic normalization
          cts_to_sig = cts_to_sig / 0.998	!Ave Cerenkov effic.  ??????
          cts_to_sig = cts_to_sig / 1.000	!Ave Calorim effic.

C convert from dsigma/dxi/dOmega to dsigma/dE'/dOmega
	  sigfact(i)=cts_to_sig
	  sigfact(i)=sigfact(i)*dxidp(pbin,deg,e1)

C get acceptance for this xi bin.
	  waccep_sum = 0.
	  domega_sum = 0.
	  do ind = mcnumtlo,mcnumthi
	    degbin = r_d * (mctlo + mctbin*(ind-0.5))
	    waccep_sum = waccep_sum + 
     >         evalaccep(xibin,xibinlo,xibinhi,degbin,e1,ep_cent)*sind(degbin)
	    domega_sum = domega_sum + sind(degbin)
	  enddo

C Keep points with (Ep < Ep_elastic) and dsig/sig<50%.
C Also require that acceptance of xi bin be > about 50% by requireing
C that the unweighted AND weighted acceptances are >50% of the peak accep.

	  if (pbin.ge.(ep_elas-fit_cut)) goto 4		!break input loop
	  acccut=0.5		!require A/A_max > 50%
	  if (abs(a_nuc-2.).le.0.0001)	acccut=0.3	!require > 30%
	write(6,*) 'acc,acccut=',(waccep_sum/domega_sum)/accep_max,acccut
	  if (cnts.ge.1 .and. (waccep_sum/domega_sum)/accep_max.gt.acccut) then
	    ep(i) = pbin
	    xi(i) = xibin
	    nu_elas = ep_elas - pbin
	    numcts(i) = cnts
	    staterr(i) = sqrt(cnts+0.5)/cnts	!+0.5???????? LAME COMPRIMIZE, that I can't even spell correctly
	    if (abs(a_nuc-2.).le.0.0001) then
	      staterr(i)=err/cnts
	    endif
	    wsq_neg(i) =  q2 - m_p*m_p - 2.*m_p*nu      !-1*(Missing mass)^2
	    xfit(i) = wsq_neg(i)
	    if (nu_elas.ge.fit_cut) n_good = i
	    i = i + 1
	  endif
	enddo
	close (unit=util)

C Check that data is binned appropriatly
	if (i.gt.1) then
	  if (abs( (pbin-ep(i-1)) - pwid ).gt.0.0001) then !not 15 MeV step.
	    write(6,*) 'Input data not binned in ',pwid*1000.,'MeV bins'
	    stop
	  endif
	endif


4	continue		!reached end of input file
	if (n_good .ne. i-1) then
	  write(6,*) 'n_pts =',i-1,' but n_good=',n_good
	  write(6,*) 'THIS IS BAD, since I replaced n_pts with n_good in the code'
	endif
	ngood = n_good		!version for new spline routine

C Initial output to log file (for this target).
	write (logfile,*) e1,' = Incident energy (GeV)'
	write (logfile,*) deg,' = scattering angle (degrees)'
	write (logfile,*) a_nuc,' = A (target)'
	write (logfile,*) zp,' = Z (target)'
	write (logfile,*) tb,' = Radiator before in radiation lengths'
	write (logfile,*) ta,' = Radiator after  in radiation lengths'
	write (logfile,*) teff,' = t_eff (effective tgt len for pair production'
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
	write (logfile,*) e_err,' = E_beam uncertainty'
	write (logfile,*) ep_err,' = Ep uncertainty'
	write (logfile,*) theta_err*d_r*1000.,' = theta uncertainty (mr)'
	write (logfile,*) ' '
	write (logfile,*) '  xi   rawA   weiA     tc    tc2     pc   totc    sys'

C Initialize some things.
	do ind = 1,n_good
	  rfit(ind) = 1.D00
	  drfit(ind) = 1.D00
	enddo

        call cubgcv(wsq_neg,rfit,drfit,ngood,yfit,cfit,max_pnts,1.d0,0,dum,wk,ier)
        if (ier.ne.0) then
          write(6,*) 'ier = ',ier,' for first call to cubgcv'
          write(6,*) 'ngood=',ngood
          write(6,*) 'wsq_neg(1:ngood)=',(wsq_neg(ind),ind=1,ngood)
          write(6,*) 'rfit(1:ngood)=',(rfit(ind),ind=1,ngood)
          write(6,*) 'drfit(1:ngood)=',(drfit(ind),ind=1,ngood)
        endif

	do i=1,n_good	!Apply all corrections/uncertainties for each E' point.


C ***** THETA BINNING CORRECTION *************************
C For each xi bin, take mcnumt theta bins. Generate dsigma/dE'/dOmega (model)
C for each xi,theta point. Convert to dsigma/dxi/dOmega, and take ratio
C of sigma(xi,thetacent) to integral of sigma(xi,thetabin)*accep(xi,thetabin)
C over the theta range (as is measured in the experiment).

! initialize some stuff.
	  xibin = xi(i)
	  pbin = evalp(xibin,deg,e1)
	  xibinlo = evalxi((pbin-pwid/2.),deg,e1)
	  xibinhi = evalxi((pbin+pwid/2.),deg,e1)

! get sigma at central angle.
	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig0 = radiate(pbin,modelhack)				!Radiated model
	  sig0 = sig0 / dxidp(pbin,deg,e1)	!dsig/dE' -> dsig/dxi
	  rawx_cent = sig0

! get integral sigma*acceptance over theta bins.
	  wxsec_sum = 0.
	  wrawx_sum = 0.
	  domega_sum = 0.
	  waccep_sum = 0.
	  do ind = mcnumtlo,mcnumthi
	    degbin = r_d * (mctlo + mctbin*(ind-0.5))
	    pbin = evalp(xibin,degbin,e1)
	    call getmodelsig(e1,degbin,tb,ta,zp,zn,m_nuc,aux,dee,
     &		pbin,modelhack,sig0,modeloutfile)

	    tmpaccep = evalaccep(xibin,xibinlo,xibinhi,degbin,e1,ep_cent)
	    wxsec_sum = wxsec_sum + sig0*tmpaccep*sind(degbin)
	    wrawx_sum = wrawx_sum + sig0*sind(degbin)
	    domega_sum = domega_sum + sind(degbin)
	    waccep_sum = waccep_sum + tmpaccep*sind(degbin)
	  enddo

! accepcorr is APPROXIMATLY the acceptance correction.  xsecn_cent/xsecn_ave
! combines BOTH the CROSS SECTION WEIGHTED acceptance and the theta binning
! correction.  accepcorr divides out the UNWEIGHTED acceptance, as a rough
! approximation.  The product of the two is still xsecn_cent/xsecn_ave.

	  bincorr(i) = min(1000. , (rawx_cent/(wxsec_sum/domega_sum)) )
	  binaccep(i) = wxsec_sum/wrawx_sum
	  tbincorr(i)=min(1000.,bincorr(i)*binaccep(i)) !totalcorr/accepcorr
	  tbincorr2(i)=rawx_cent/(wxsec_sum/(waccep_sum+.0000001)+.0000001)

	  sigfact(i)=sigfact(i)*bincorr(i)


C ***** XI BINNING CORRECTION ****************************
C For each xi bin, compare sig-central to ave. over 21 xi bins.
C NOTE!!!!!!  Should also take acceptance VARIATIONS over xi into account.
C THIS IS NOT DONE YET, but should be a small effect.  Applying the
C acceptance variation over such small steps just introduces the single-bin
C statistics of the acceptance MC.  Treat it as if accep. flat over xi bin.

	  rawx2_sum=0
	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  do ind = 0,20
	    xi2 = xibinlo + ind/20.*(xibinhi-xibinlo)
	    e2 = evalp(xi2,deg,e1)
	    call getmodelsig(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2,modelhack,sig0,modeloutfile)
	    rawx2_sum=rawx2_sum+sig0
	    if (ind.eq.10) rawx2_cent=sig0
	  enddo

	  pbincorr(i)=rawx2_cent/(rawx2_sum/21.)
	  sigfact(i)=sigfact(i)*pbincorr(i)

C ***** COULOMB CORRECTION - ENERGY BOOST ****************
C for each bin, calculate cross section at measured kinematics and at
C the vertex kinematics (+dE to both e and e').  Measured value is at offset
C kinematics, so correct by ratio of central kinematics to offset kinematics.

	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig0 = radiate(ep(i),modelhack)
	  call radiate_init(e1+coulomb_de,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig1 = radiate(ep(i)+coulomb_de,modelhack)
	  coulcor_e(i)=sig0/sig1
	  sigfact(i)=sigfact(i)*coulcor_e(i)

C ***** COULOMB CORRECTION - DEFLECTION ******************
C For each xi bin, average cross section from 21 theta bins (due to range
C of deflection of electron).  Combining e,e' in quadriture (just a guess).

	  coulomb_theta_e = coulomb_pt / e1	!max. angular deflection of beam
	  coulomb_theta_ep = coulomb_pt / ep(i)	!max. angular deflection of e'
	  coulomb_theta = sqrt(coulomb_theta_e**2 + coulomb_theta_ep**2)
	  coulomb_theta = coulomb_theta * r_d

	  rawx2_sum=0.
	  e2=ep(i)
	  do ind = -10,10
	    degbin = deg + (float(ind)/10.)*coulomb_theta
	    call getmodelsig(e1,degbin,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2,modelhack,sig0,modeloutfile)
	    rawx2_sum=rawx2_sum+sig0
	    if (ind.eq.0) rawx2_cent=sig0
	  enddo

	  coulcor_th(i)=rawx2_cent/(rawx2_sum/21.)
	  sigfact(i)=sigfact(i)*coulcor_th(i)


C****** UNCERTAINTY IN CROSS SECTION DUE TO UNCERTAINTY IN E1,E2,THETA ******
C Take difference between central cross section and cross section with a
C one sigma offset in E,E',theta as the one sigma uncertainty in the xsec.

	  e2=ep(i)
! Central Cross Section.
	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig0 = radiate(e2,modelhack)
	  sig0 = sig0 / dxidp(e2,deg,e1)
	  call getmodelsig(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2,modelhack,sig1,modeloutfile)

	  if (abs(sig0-sig1)/((sig0+sig1)/2.).gt.0.00001) then
	    write(6,*) 'sig0,sig1=',sig0,sig1
	    write(6,*) 'e1,e2,deg=',e1,e2,deg
	    stop
	  endif
	  sig0 = sig1	!getmodelsig value

	  if (sig0.eq.0) sig0=1.e-22	!aviod divide by zero.

! Beam Energy.
	  call radiate_init(e1*(1+e_err),deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  sig1 = radiate(e2,modelhack)
	  sig1 = sig1 / dxidp(e2,deg,e1*(1+e_err))
!	  call getmodelsig(e1*(1+e_err),deg,tb,ta,zp,zn,m_nuc,aux,dee,
!     &		e2,modelhack,sig1,modeloutfile)
	  dsig_e = abs(sig1-sig0)/sig0

! Spectrometer Momentum.
!	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
!	  sig1 = radiate(e2*(1+ep_err),modelhack) 
!	  sig1 = sig1 / dxidp(e2*(1+ep_err),deg,e1*(1+e_err))
	  call getmodelsig(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2*(1+ep_err),modelhack,sig1,modeloutfile)
	  dsig_ep = abs(sig1-sig0)/sig0

! Spectrometer Angle.
!	  call radiate_init(e1,deg+theta_err,tb,ta,zp,zn,m_nuc,aux,dee)
!	  sig1 = radiate(e2,modelhack)
!	  sig1 = sig1 / dxidp(e2,deg+theta_err,e1)
	  call getmodelsig(e1,deg+theta_err,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2,modelhack,sig1,modeloutfile)
	  dsig_theta = abs(sig1-sig0)/sig0

C ***** APPLY SYSTEMATIC UNCERTAINTIES & EFFICIENCY CORRECTIONS ******
C syserr is fractional uncertainty, all pieces added in quadriture.

	  syserr(i) = 0.	!FRACTIONAL uncertainty.
	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Accep = 1.0% & 4%*(1-A)
	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Bincent = 1% & .1*corr
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !PID Effic.
	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Q measurement
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !Track eff.
	  syserr(i) = sqrt(syserr(i)**2 + 0.005**2) !Tracking cuts.
	  syserr(i) = sqrt(syserr(i)**2 + 0.025**2) !Radcorr = 2.5% 
	  syserr(i) = sqrt(syserr(i)**2 + dsig_e**2)	 !error from e offset
	  syserr(i) = sqrt(syserr(i)**2 + dsig_ep**2)	 !error from ep offset
	  syserr(i) = sqrt(syserr(i)**2 + dsig_theta**2) !error from theta offset

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
	  endif

! Hydrogen is consitant within systematics ==> don't increase error for it.
!	  syserr(i) = sqrt(syserr(i)**2 + 0.010**2) !Hydrogen Norm.

C assign uncertainty to acceptance correction.  Base 1% offset (already
C added in) + .04*(1-A). Gives 2.2% uncertainty at A=.5 (sqrt(1.0**2+2.0**2).
	  junk = (waccep_sum/domega_sum)/accep_max !approx. accep/accep_max
	  syserr(i)=sqrt( syserr(i)**2 + (0.04*(1-junk))**2 )

C assign 10%*(1-tbincorr2) uncertainty to tbincorr.  Max. correction is 20%
C except at 15 degress, extremely lo nu (nu <= 105 MeV).
	  syserr(i)=sqrt( syserr(i)**2 + (0.10*(1-tbincorr2(i)))**2 )


C ***** SUBTRACTION OF CHARGE-SYMMETRIC BACKGROUND *************************
C Subtract off counts based on fit to e+ production cross section, assuming
C the cross section is function of xi, and flat in theta.  Convolve xsec
C with acceptance over theta bins.  ONLY FOR 55 DEGREE DATA!!!

	  if (abs(deg-55).lt.0.01) then

! get integral sigma*acceptance over theta bins. sigma=exp(-17.85-11.22*xi)
! is dsigma/dOmega/dE'

	    xibin = xi(i)
	    pbin = evalp(xibin,deg,e1)
	    xibinlo = evalxi((pbin-pwid/2.),deg,e1)
	    xibinhi = evalxi((pbin+pwid/2.),deg,e1)

	    posxsec = exp(-17.85-11.22*(xibin-1))*a_nuc*teff
	    fracerr = 1.5*(2.0425+exp(6.0477*(xibin-.5534)))/100
	    fracerr = max(.06,fracerr)

	    num_pos = 0.
	    lumin = eta*(q_tot/1.602d-19)/corr_factor/1.d+33*0.998
	    do ind = mcnumtlo,mcnumthi
	      degbin = r_d * (mctlo + mctbin*(ind-0.5))
	      tmpaccep = evalaccep(xibin,xibinlo,xibinhi,degbin,e1,ep_cent)
	      tmpdomegade = (1000.*pwid)*sind(degbin)*mctbin*(mcpmax-mcpmin)
	      num_pos = num_pos + tmpaccep*posxsec*tmpdomegade*lumin
	    enddo
	    dnum_pos = num_pos*fracerr
	write(6,*) 'e-,e+,de+=',nint(numcts(i)),num_pos,dnum_pos

!subtract counts, and take care of stat/sys errors.
	    dnumstat = staterr(i)*numcts(i)
	    dnumsys = syserr(i)*numcts(i)
	    dnumsys = sqrt ( dnumsys**2 + dnum_pos**2)
	    numcts(i) = numcts(i) - num_pos
	    staterr(i) = dnumstat / numcts(i)
	    syserr(i) = dnumsys / numcts(i)

	  endif 		!if 55 degrees.


C ***** OUTPUT CORRECTIONS/UNCERTAINTIES *****
	  write(6,'(a,5(f6.2,a))') 'dsig/sig (ep,e,theta,coul_dE,coul_th)=',
     >		nint(10000*dsig_ep)/100.,'%',
     >		nint(10000*dsig_e)/100.,'%',
     >		nint(10000*dsig_theta)/100.,'%',
     >		nint(10000*(coulcor_e(i)-1))/100.,'%',
     >		nint(10000*(coulcor_th(i)-1))/100.,'%'
	  write(logfile,'(a,5(f6.2,a))') 'dsig/sig (ep,e,theta,coul_dE,coul_th)=',
     >		nint(10000*dsig_ep)/100.,'%',
     >		nint(10000*dsig_e)/100.,'%',
     >		nint(10000*dsig_theta)/100.,'%',
     >		nint(10000*(coulcor_e(i)-1))/100.,'%',
     >		nint(10000*(coulcor_th(i)-1))/100.,'%'

	  write(6,'(a,f5.3,7f7.4)') 'p,rawA,weiA,tc,tc2,pc,totc,sys=',ep(i),
     >		(waccep_sum/domega_sum)/accep_max,	! accepcorr(nowei)
     >		binaccep(i)/accep_max,			! accepcorr(wei)
     >		tbincorr(i),				! tbincorr
     >		tbincorr2(i),				! tbincorr
     >		pbincorr(i),				! pbincorr
     >		1/(bincorr(i)*pbincorr(i)*accep_max),	! totcorr
     >		syserr(i)
	  write(logfile,'(f5.3,7f7.4)') ep(i),
     >		(waccep_sum/domega_sum)/accep_max,	! accepcorr(nowei)
     >		binaccep(i)/accep_max,			! accepcorr(wei)
     >		tbincorr(i),				! tbincorr
     >		tbincorr2(i),				! tbincorr
     >		pbincorr(i),				! pbincorr
     >		1/(bincorr(i)*pbincorr(i)*accep_max),	! totcorr
     >		syserr(i)

	enddo

! Output counts and correction factors needed for conversion to cross section
! and for radiative corrections

	tmpfile = file(1:filelen)//'.preradcor'
	open(unit=util,name=tmpfile,status='unknown')

	write(util,*) e1,'	= Beam Energy (GeV)'
	write(util,*) deg,'	= Scattering Angle (degrees)'
	write(util,*) zp,'	= Z'
	write(util,*) zn,'	= N'
!	write(util,*) tb,'	= tb (% radiation length)'
!	write(util,*) ta,'	= ta (% radiation length)'

! Results to cross section file. Cross section units are
! nB/(Mev-ster). Only 'good' data points are output.

	write(logfile,*) ' '
	write(logfile,*) '  Ep           N dN/Nstat dN/Nsys      sigfact       tb       ta       ctstosig'
	write(util,*)    '  Ep           N dN/Nstat dN/Nsys      sigfact       tb       ta       ctstosig'
	do i = 1,n_good
	  e2 = ep(i)
	  nu = e1 - e2
	  q2 = 4.*e1*e2*snsq
	  x = q2/(2.*m_p*nu)
	  if (x.lt.a_nuc) then
	    write(util,1006) e2,numcts(i),staterr(i),syserr(i),sigfact(i),tb,ta,cts_to_sig
	    write(logfile,1006) e2,numcts(i),staterr(i),syserr(i),sigfact(i),tb,ta,cts_to_sig
	  endif
	enddo
1006	format (f6.3,f12.3,2f8.5,f15.12,2f9.6,f15.12)
	close (unit=util)

C End of file jump destinations.

992	close (unit=logfile)
	stop ' '

990	type *, ' RADCOR_NE18: Error reading file = ',tmpfile
991	stop 'Analysis aborted!'

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

	real*8 evalp

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


C getmodelsig gets the radiated sigma.  It first checks the lookup table.  If it
C doesn't find the cross section for the given kinematics, it calls
C radiate_init and then radiate.  This is somewhat slower due to calling
C radiate_init if the model file does not exist (how much slower - one call
C to ddilog - how slow is that?).  Also checks (checkfrac) of the model values
C by calculating by hand.
C NOTE: model is dsigma/dOmega/dxi, NOT dsigma/dOmega/dE'
C                              ^^^                    ^^^
C though it is given as a function of E',theta.

	subroutine getmodelsig(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee,
     &		e2,modelhack,sig,modeloutfile)

	real*8 e1,deg,tb,ta,zp,zn,m_nuc,aux(7),dee	!args. for radiate_init
	real*8 e2,modelhack				!args. for radiate
	real*8 sig					!output
	integer*4 modeloutfile				!model output file.

	integer*4 maxmodpts
	parameter (maxmodpts=3000)
	integer*4 modpts
	real*8 emod
	real*8 pmod(maxmodpts)
	real*8 degmod(maxmodpts)
	real*8 sigmod(maxmodpts)

	common/modelblock/pmod,degmod,sigmod,emod,modpts

	real*8 tablesig,calcsig

	logical foundsig,checksig

	integer*4 i,istart
	data i /1/

	integer*4 iran/19857141/
	real*8 ran,dxidp

C Look for cross section in lookup table.  Skip table if e1 <> emod.
	istart=i
	foundsig=.false.
	if (modpts.gt.0 .and. abs(e1-emod).le.0.00001) then	!lookup table exists.
	  do while (.not.foundsig)
	    if (abs(pmod(i)-e2).lt.0.00001) then
	      if (abs(degmod(i)-deg).lt.0.00001) then	!found sigma.
	        tablesig = sigmod(i)
	        foundsig = .true.
	      endif
	    endif
	    i=i+1
	    if (i.gt.modpts) i=1	!reset counter
	    if (i.eq.istart) goto 55
	  enddo
	endif
55	continue		!break loop if didn't find

C If value is not in lookup (or if checking this value) then calculate sigma.

	checksig = (ran(iran).le.0.02)	!check 2% of lookup table values
	if (.not.foundsig .or. checksig) then
	  call radiate_init(e1,deg,tb,ta,zp,zn,m_nuc,aux,dee)
	  calcsig = radiate(e2,modelhack)
	  calcsig = calcsig / dxidp(e2,deg,e1)	!dsig/dE' -> dsig/dxi
	endif

	if (foundsig .and. checksig) then	!check against real model.
	  if (calcsig.gt.0 .or. tablesig.gt.0) then   !avoid div. by zero.
	    if (abs(calcsig-tablesig)/((calcsig+tablesig)/2).gt.0.0001) then ! >.1% off
	      write(6,*) 'BAD MODEL VALUE!!!!!!!!!!'
	      write(6,*) 'xi,deg,calcsig,tablesig=',xibin,deg,calcsig,tablesig
	      stop
	    endif
	  endif
	endif

	if (foundsig) then
	  sig=tablesig
	else
	  sig=calcsig
	endif

	if (abs(e1-emod).lt.0.00001) write(modeloutfile,'(2f10.5,e15.6)') e2,deg,sig

	return
	end
