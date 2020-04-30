!	Solves the set of n linear equations A.X = B. Here a is input, not as the matrix A but
!	rather as its LU decomposition, determined by the routine ludcmp. indx is input as the
!	permutation vector returned by ludcmp. b(1:n) is input as the right-hand side vector B,
!	and returns with the solution vector X. a, n, np, and indx are not modiffied by this routine
!	and can be left in place for successive calls with different right-hand sides b. This routine
!	takes into account the possibility that b will begin with many zero elements, so it is efficient
!	for use in matrix inversion.
!	Amir Shahmoradi, Oct 18, 2009, 1:54 PM, MTU
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER, INTENT(IN) :: n,np
	  INTEGER, DIMENSION(n), INTENT(IN) :: indx
      DOUBLE PRECISION, DIMENSION(n), INTENT(INOUT) :: b
      DOUBLE PRECISION, DIMENSION(np,np), INTENT(IN) :: a
      INTEGER :: i,ii,j,ll
      DOUBLE PRECISION :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.0d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      END SUBROUTINE lubksb