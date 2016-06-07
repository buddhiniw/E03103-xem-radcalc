	subroutine getkine(runno,e1,ep_cent,deg,theta_tar,tgt_len,tgt_thick,
     >			   corr_factor,zp,a_nuc,m_atom,q_tot,kinefile)

	implicit none

	integer*4 runno

	real*8 e1,ep_cent,deg,theta_tar,tgt_len,tgt_thick
	real*8 corr_factor,zp,a_nuc,m_atom,q_tot

	character*80 line,tmpline
	character*80 tmpfile,kinefile

C ============================ Executable Code =================================

	write(tmpfile,'(''kine/'',i4,''.kine'')') runno
!	open(unit=80,name=tmpfile,status='old',readonly)

	open(unit=80,name=kinefile,status='old',readonly)

55      read(80,'(a)',end=66) line

	if (line(1:6).eq.'E_beam') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) e1
!	  write(6,*) 'e1=',e1
* gebeam from engine has energy loss alread done, must set back to 4.045
!	  write(6,*) '*****SETTING E_BEAM TO 4.045 GeV********'
	  e1=4.045
	else if (line(1:1).eq.'P') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) ep_cent
!	  write(6,*) 'ep_cent=',ep_cent
	else if (line(1:5).eq.'Theta') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) deg
!	  write(6,*) 'deg=',deg
	  if (abs(deg-nint(deg)).gt.0.011) then
	    write(6,*) 'KINEMATICS FILE HAS THETA=',deg
	    deg=nint(deg)
	    write(6,*) '*****SPECTROMETER ANGLE BEING FORCED TO ',deg,'********'
	  endif
	else if (line(1:3).eq.'tgt') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) zp
	  tmpline=tmpline(index(tmpline,':')+1:)
	  read (tmpline,*) a_nuc
	  tmpline=tmpline(index(tmpline,':')+1:)
	  read (tmpline,*) m_atom
!	  write(6,*) 'zp,a_nuc,m_atom=',zp,a_nuc,m_atom
	else if (line(1:5).eq.'thick') then
	  tmpline=line(index(line,',')+1:)
	  read (tmpline,*) tgt_len
	  tmpline=tmpline(index(tmpline,',')+1:)
	  read (tmpline,*) tgt_thick
!	  write(6,*) 'tgt_len,tgt_thick=',tgt_len,tgt_thick
	else if (line(1:11).eq.'gtarg_theta') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) theta_tar
!	  write(6,*) 'targ_ang=',theta_tar
	else if (line(3:6).eq.'CORR') then
	  tmpline=line(index(line,'=')+1:)
	  read (tmpline,*) corr_factor
!	  write(6,*) 'corr_factor=',corr_factor
	else if (line(1:6).eq.'Q_tot(') then
	  tmpline=line(index(line,':')+1:)		!bcm1
	  tmpline=tmpline(index(tmpline,':')+1:)	!bcm2
	  read (tmpline,*) q_tot
	  q_tot=q_tot/1.e6		!microC --> Coulombs
!	  write(6,*) 'q_tot=',q_tot
	endif
	goto 55			!read next line

66	close(80)

* Fix up the target lengths that were incorrect in original replay.

	if (abs(tgt_thick-.890).le.0.001) then		!fix up C2
	  tgt_thick=.8915
	else if (abs(tgt_thick-2.50).le.0.001) then	!fix up C6
	  tgt_thick=2.5096
	else if (abs(tgt_thick-0.376).le.0.001) then	!fix up AU6
	  tgt_thick=0.3768
	else if (abs(tgt_len-3.9972).lt.0.001) then	!fix up 4cm hydrogen.
	  tgt_len=4.36
	  tgt_thick=0.315228    			!4.38cm*.0723g/cm^3
	endif

* Theta_tar is usually given as 0, +/- 20 degrees. Sometimes, it's 90. treat
* 90 as zero.  Also, set to zero for cryotargets.

	if (abs(theta_tar-90).lt.0.01) theta_tar=0
	if (a_nuc.le.4) theta_tar=0     !some cryo runs have theta<>0 in kine.

	return
	end
