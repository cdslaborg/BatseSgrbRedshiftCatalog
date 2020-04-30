! Amir Shahmoradi, Sunday 6:12 AM, April 22, 2012, IFS, UTEXAS
SUBROUTINE qromb(func,a,b,ss)
  INTEGER JMAX,JMAXP,K,KM
  DOUBLE PRECISION a,b,func,ss,EPS
  EXTERNAL func
  PARAMETER (EPS=1.d-3, JMAX=20, JMAXP=JMAX+1, K=4, KM=K-1)
  INTEGER j
  DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
  h(1)=1.d0
  do j=1,JMAX
    call trapzd(func,a,b,s(j),j)
    if (j.ge.K) then
      call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
      if (dabs(dss).le.EPS*dabs(ss)) return
    endif
    s(j+1)=s(j)
    h(j+1)=0.25d0*h(j)
  end do
  write(*,*) 'too many steps in qromb, Press enter to stop...'
  read(*,*)
  STOP
END SUBROUTINE qromb
