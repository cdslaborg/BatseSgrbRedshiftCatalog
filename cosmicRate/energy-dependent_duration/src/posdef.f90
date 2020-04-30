FUNCTION posdef(a,n)
  INTEGER n
  DOUBLE PRECISION a(n,n),p(n)
  INTEGER i,j,k
  DOUBLE PRECISION sum
  INTEGER :: posdef
  posdef=1
  do i=1,n
    do j=i,n
      sum=a(i,j)
      do k=i-1,1,-1
        sum=sum-a(i,k)*a(j,k)
      end do
      if(i.eq.j)then
        if(sum.le.0.) then 
          write (*,*) 'choldc failed'
          posdef=0
        end if
        p(i)=dsqrt(sum)
      else
        a(j,i)=sum/p(i)
      endif
    end do
  end do
END FUNCTION posdef
