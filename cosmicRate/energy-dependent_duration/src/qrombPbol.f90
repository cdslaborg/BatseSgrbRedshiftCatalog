! Last modified by: Amir Shahmoradi, Friday 7:47 PM, Dec 13, 2013, IFS/ICMB, UT Austin
SUBROUTINE qrombPbol(func,a,b,ss)

      INTEGER JMAX,JMAXP,K,KM
      DOUBLE PRECISION a,b,func,ss,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-5, JMAX=20, JMAXP=JMAX+1, K=6, KM=K-1)
      INTEGER j
      DOUBLE PRECISION dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do j=1,JMAX
        call trapzdPbol(func,a,b,s(j),j)
        if (j.ge.K) then
          call polintPbol(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (dabs(dss).le.EPS*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
	  end do
      write(*,*) 'too many steps in qrombPbol, Press enter to stop...'
	  read(*,*)
	  STOP
END SUBROUTINE qrombPbol

!	--------------------------

SUBROUTINE polintPbol(xa,ya,n,x,y,dy)
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
          if(den.eq.0.) then
			write(*,*) 'failure in polintPbol'
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
      END SUBROUTINE polintPbol
!	-----------------------
      SUBROUTINE trapzdPbol(func,a,b,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.d0
        do j=1,it
          sum=sum+func(x)
          x=x+del
		end do
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif

END SUBROUTINE trapzdPbol
