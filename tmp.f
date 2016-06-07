	function fz_table(z,beta)

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
	        stop
	      endif
	    endif
	  enddo
	endif
	if (abs(w1+w2-1).gt.0.0001) then
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
==============================================
	subroutine fz_table_read

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
	enddo
	return
	end
