!  Amir Shahmoradi, Oct 16, 2009, 11:14 AM
SUBROUTINE samcovmat(n,p,D,S)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,p
  DOUBLE PRECISION, INTENT(IN) :: D(n,p)
  DOUBLE PRECISION, INTENT(OUT) :: S(p,p)
  INTEGER :: i,j,k
  DOUBLE PRECISION, DIMENSION(n,p) :: avg
  DOUBLE PRECISION, DIMENSION(n) :: dummy
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