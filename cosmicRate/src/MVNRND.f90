!	Given the mean vector MU and the covariance matrix COV, this subroutine generates a random vector X (of length n>=2)
!	from an n-dimensional multivariate normal distribution.
!	First a vector of n standard normal random deviates is generated,
!	then this vector is multiplied by the Cholesky decomposition of the covariance matrix,
!	then the vector MU is added to the resulting vector, and is stored in the output vector X.
!	ATTENTION: Only the upper half of the covariance matrix (plus the diagonal terms) need to be given in the input.
!	In the ouput, the upper half and diagonal part will still be the covariance matrix, while the lower half will be
!	the Cholesky decomposition elements (excluding its diagonal terms that are provided only in the vector DIAG).

!	USES choldc.f90, gasdevran2.f90
!	Amir Shahmoradi, March 22, 2012, 2:21 PM, IFS, UTEXAS
	SUBROUTINE MVNRND(idum,n,MU,COV,X)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(INOUT) :: idum
	DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: MU
	DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: X
	DOUBLE PRECISION, DIMENSION(n,n), INTENT(INOUT) :: COV
	DOUBLE PRECISION, DIMENSION(n) :: DIAG,VECTOR
	DOUBLE PRECISION :: gasdevran2
	INTEGER :: i,j
	call choldc(COV,n,n,DIAG)
	do i=1,n
		VECTOR(i)=gasdevran2(idum)
		X(i)=VECTOR(i)*DIAG(i)
	end do
	do i=2,n
		X(i)=X(i)+dot_product(COV(i,1:i-1),VECTOR(1:i-1))
	end do
	X=X+MU
	END SUBROUTINE MVNRND