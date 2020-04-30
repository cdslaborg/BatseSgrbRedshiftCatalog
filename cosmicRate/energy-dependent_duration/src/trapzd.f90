SUBROUTINE trapzd(func,a,b,s,n)
  INTEGER n
  DOUBLE PRECISION a,b,s,func
  EXTERNAL func
  INTEGER it,j
  DOUBLE PRECISION del,sum,tnm,x
  if (n.eq.1) then
    s=0.5*(b-a)*(func(a)+func(b))
  else
    it=2**(n-2)
    tnm=it
    del=(b-a)/tnm
    x=a+0.5*del
    sum=0.
    do j=1,it
      sum=sum+func(x)
      x=x+del
		end do
    s=0.5*(s+(b-a)*sum/tnm)
  endif
END SUBROUTINE trapzd
