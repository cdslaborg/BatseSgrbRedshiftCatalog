!	This subroutine finds the determinant of a given n*n matrix S with physical dimension of np*np
!	Amir Shahmoradi, Oct 18, 2009, 4:10 PM, MTU
!	USES ludcmp.f90
	FUNCTION determinant(n,np,S)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,np
	DOUBLE PRECISION, DIMENSION(np,np), INTENT(IN) :: S
	DOUBLE PRECISION :: determinant
	INTEGER, DIMENSION(np) :: indx
	DOUBLE PRECISION, DIMENSION(np,np) :: Sdummy
	INTEGER :: j
	Sdummy = S
	call ludcmp(Sdummy,n,np,indx,determinant) !	This returns determinant as +-1.
	do j=1,n
		determinant = determinant*Sdummy(j,j)
	end do
	END FUNCTION determinant