! Amir Shahmoradi, Tuesday July 10, 2012, 3:57 PM, IFS, The University of Texas at Austin.
! luminosity distance in units of Mpc.
FUNCTION ldisWickram(z)
  USE COSMOparameters
  USE constants, ONLY: CoverH
  DOUBLE PRECISION, INTENT(IN) :: z
  DOUBLE PRECISION :: ldisWickram,zz,term,alpha,alpha0,x,x0,psi,psi0
  alpha=1.d0+2.d0*omega_DE/(omega_DM*(1.d0+z)**3)
  alpha0=1.d0+2.d0*omega_DE/omega_DM
  x=dlog(alpha+dsqrt(alpha*alpha-1.d0))
  x0=dlog(alpha0+dsqrt(alpha0*alpha0-1.d0))
  psi=x**0.33333333*(1.5874010519682-6.2992105236833d-3*x*x+7.5375168659459d-5*x**4)
  psi0=x0**0.33333333*(1.5874010519682-6.2992105236833d-3*x0*x0+7.5375168659459d-5*x0**4)
  ldisWickram=CoverH*(1.d0+z)*(psi0-psi)/(omega_DE**0.16666666*omega_DM**0.33333333)
END FUNCTION ldisWickram