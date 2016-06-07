C function evalxi.  gives xi(p,theta,ebeam), with ebeam&p in GeV, theta in degrees

	function evalxi(p,theta,ebeam)

	implicit none

	real*8 ebeam,p,theta
	real*8 nu,qsq,x,evalxi

	real*8 m_p
	parameter (m_p = .93827231)

	if (p.ge.ebeam) then
	  evalxi=100
	else
	  nu=ebeam-p
	  qsq=4.*ebeam*p*sind(theta/2)**2
	  x=qsq/2./m_p/nu
	  evalxi=2.*x/(1+sqrt(1+4.*m_p**2*x**2/qsq))
	endif

	return

	end

C evalp gives p(xi,theta,ebeam). ebeam&p in GeV, theta in degrees.
C p=m*xi/2 * (m*xi+2*e)/(m*xi+2*e*sin(theta/2)**2)
	function evalp(xi,theta,ebeam)

	implicit none

	real*8 xi,ebeam,theta
	real*8 evalp
	real*8 tmp1,tmp2

	real*8 m_p
	parameter (m_p = .93827231)

	tmp1 = m_p*xi + 2.*ebeam
	tmp2 = m_p*xi + 2.*ebeam*sind(theta/2.)**2
	evalp = m_p*xi/2. * tmp1/tmp2
	return

	end

C dxidp gives dxi/dp(p,theta,ebeam). ebeam&p in GeV, theta in degrees.
C   since this is always used to convert cross sections, convert units to
C   dxi(unitless)/dp(MeV).
	function dxidp(p,theta,ebeam)

	real*8 ebeam,p,theta
	real*8 nu,s,qsq,z,dxidp

	real*8 m_p
	parameter (m_p = .93827231)

	nu = ebeam - p
	s = sind(theta/2)**2
	qsq = 4.*ebeam*p*s
	z = nu + sqrt(nu**2 + qsq)

	dxidp = qsq/m_p/z**2 * ( z/p + (z - 2.*ebeam*s)/(z - nu) )

	dxidp = dxidp/1000.

	return
	end
