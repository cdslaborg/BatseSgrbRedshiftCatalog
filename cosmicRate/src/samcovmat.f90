!	This subroutine finds the elements of the symmetric p*p sample covariance matrix of a given set of
!	N observations, each one with p parameters in the format of a "" N*p "" matrix.
!	For a review refer to Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
!	Also refer to Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
!	Amir Shahmoradi, Oct 16, 2009, 11:14 AM, MTU
	SUBROUTINE samcovmat(n,p,D,S)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n,p
!	n is the number of observations, p is the number of parameters for each observation
	DOUBLE PRECISION, INTENT(IN) :: D(n,p)
	DOUBLE PRECISION, INTENT(OUT) :: S(p,p)
!	D is the matrix of the data, S contains the elements of the sample covariance matrix
	INTEGER :: i,j,k
!	DOUBLE PRECISION :: sigma
	DOUBLE PRECISION, DIMENSION(n,p) :: avg
	DOUBLE PRECISION, DIMENSION(n) :: dummy
!	avg is the average of each of the p parameters
!	write(*,*) 'Hello!'
	do i = 1,n
		do j = 1,p
			avg(i,j) = sum(D(1:n,j))/dble(n)
		end do
	end do
	do i = 1,p
		do j = 1,p
			s(i,j) = dot_product(D(1:n,i)-avg(1:n,i),D(1:n,j)-avg(1:n,j))/dble(n-1)
		end do
	end do
	END SUBROUTINE samcovmat