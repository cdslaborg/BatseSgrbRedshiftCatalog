SUBROUTINE qromo(func,a,b,ss,choose)
  INTEGER JMAX,JMAXP,K,KM
  DOUBLE PRECISION :: a,b,ss,EPS
  DOUBLE PRECISION, EXTERNAL :: func
  PARAMETER (EPS=1.0d-3, JMAX=50, JMAXP=JMAX+1, K=5, KM=K-1)
  INTEGER j
  DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
  h(1)=1.0d0
  do j=1,JMAX
    call choose(func,a,b,s(j),j)
    if (j.ge.K) then
      call polint(h(j-KM),s(j-KM),K,0.0d0,ss,dss)
      if (DABS(dss).le.EPS*DABS(ss)) return
    endif
    s(j+1)=s(j)
    h(j+1)=h(j)/9.0d0
	  end do
  WRITE(*,*) 'too many steps in qromo.f90'
	  STOP
END
