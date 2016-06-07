	subroutine model_fix(mmsq_neg,sigma)
C+______________________________________________________________________________
!
! MODEL_FIX - Multiply SIGMA by (interpolated) function of missing
!   mass squared. The function is a cubic spline, evaluated by the
!   IMSL routine ICSEVU, which must have ordinates arranged in
!   increasing order. Because of this, we actually use the NEGATIVE
!   of the missing mass, squared.
!
! ARGUMENTS:
!
!   MMSQ_NEG:	(R*8):	The NEGATIVE of the missing mass squared at
!	    		which SIGMA was evaluated. (Input).
!   SIGMA:	(R*8):	The cross section. (Input and output).
!
! Modification history:
!
!   1/5/88 - (DHP) Reworked code to use new IMSL library routines.
!
!   10/11/93 - (JRA) use spline value at edge of data instead of spline
!		fit when W squared is outside of the data range. (otherwise,
!		use values from a cubic fit, which go to +/- large numbers.
C-______________________________________________________________________________

	implicit none
	include 'include.inc'

	real*8 tmp_mmsq_neg,min_wsq_neg,max_wsq_neg
	real*8 mmsq_neg,sigma

C!!	real*8 tmpfactor,dcsval
C!!	real*8 brkpts(max_pnts),cs(4,max_pnts)

	real*8 factor,evalspline
	real*8 xfit(max_pnts),yfit(max_pnts),cfit(max_pnts,3)

	integer*4 ngood
C!!	integer*4 n_good
C!!	common /modfix/ brkpts,cs,n_good
	common /newmodfix/ xfit,yfit,cfit,ngood
	common /wsq_limits/ min_wsq_neg,max_wsq_neg

C ============================ Executable Code =================================

C Evaluated cubic spline function at mmsq.

	tmp_mmsq_neg=min(mmsq_neg,max_wsq_neg)		!the fit can be < 0 outside of the range with
	tmp_mmsq_neg=max(tmp_mmsq_neg,min_wsq_neg)	! data and this can mess up the integrals.

	factor = evalspline(tmp_mmsq_neg,ngood,xfit,yfit,cfit)
C!!	tmpfactor = dcsval(tmp_mmsq_neg,n_good-1,brkpts,cs)
C!!	write(67,*) tmp_mmsq_neg,tmpfactor,factor
C!!	if (abs(tmpfactor-factor).ge.0.02*(tmpfactor+factor)) then
C!!	  write(69,*) '>4% diff. in model_fix: wsq,dcsval,eval=',real(tmp_mmsq_neg),
C!!     &       real(tmpfactor),real(factor)
C!!	endif

C Modify sigma, and go home.

	sigma = sigma*factor

!	if (mmsq_neg.lt.(-1.3) .or. mmsq_neg.gt.(0.3)) then
!	  write(99,*) tmp_mmsq_neg,mmsq_neg,factor
!	else if (abs(100.*mmsq_neg-int(100.*mmsq_neg)) .lt. 0.1) then
!	  write(99,*) tmp_mmsq_neg,mmsq_neg,factor
!	endif
	return
	end
