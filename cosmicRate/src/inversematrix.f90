!	This code replaces returns the inverse matrix Y of a n*n matrix S of physical dimension np*np.
!	Amir Shahmoradi, Oct 18, 2009, 1:54 AM, MTU
!	USES ludcmp.f90 & lubksb.f90
	SUBROUTINE inversematrix(n,np,S,Y)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,np
	DOUBLE PRECISION, DIMENSION(np,np), INTENT(IN) :: S
	DOUBLE PRECISION, DIMENSION(np,np), INTENT(OUT) :: Y
	DOUBLE PRECISION, DIMENSION(np,np) :: Sdummy
	INTEGER :: i,j,indx(np)
	DOUBLE PRECISION :: d
	do i = 1,n
		do j = 1,n
			Y(i,j) = 0.0d0
		end do
		Y(i,i) = 1.0d0
	end do
	Sdummy = S
	call ludcmp(Sdummy,n,np,indx,d)
	do j = 1,n
		call lubksb(Sdummy,n,np,indx,Y(1:n,j))
	end do
	END SUBROUTINE inversematrix