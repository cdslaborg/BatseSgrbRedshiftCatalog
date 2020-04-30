SUBROUTINE polint(xa,ya,n,x,y,dy)
  INTEGER n
  DOUBLE PRECISION dy,x,y,xa(n),ya(n)
  INTEGER, PARAMETER :: NMAX=10
  INTEGER i,m,ns
  DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=dabs(x-xa(1))
  do i=1,n
    dift=dabs(x-xa(i))
    if (dift.lt.dif) then
      ns=i
      dif=dift
    endif
    c(i)=ya(i)
    d(i)=ya(i)
  end do
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    do i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
      if(den.eq.0.d0) then
			write(*,*) 'failure in polint'
			STOP
	  end if
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
		end do
    if (2*ns.lt.n-m)then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    endif
    y=y+dy
  end do
END SUBROUTINE polint
