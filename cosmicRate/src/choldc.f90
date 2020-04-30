SUBROUTINE choldc(a,n,np,p)

    INTEGER n,np
    DOUBLE PRECISION a(np,np),p(n)
    INTEGER i,j,k
    DOUBLE PRECISION sum
    do i=1,n
		do j=i,n
			sum=a(i,j)
			do k=i-1,1,-1
				sum=sum-a(i,k)*a(j,k)
			end do
			if(i.eq.j)then
				if(sum.le.0.) then 
					write (*,*) 'choldc failed'
					STOP
				end if
				p(i)=dsqrt(sum)
			else
				a(j,i)=sum/p(i)
			endif
		end do
	end do

END SUBROUTINE choldc
