SUBROUTINE midexp(funk,aa,bb,s,n)
  INTEGER n
  DOUBLE PRECISION aa,bb,s,funk
  EXTERNAL funk
  INTEGER it,j
  DOUBLE PRECISION ddel,del,sum,tnm,x,func,a,b
  func(x)=funk(-dlog(x))/x
  b=dexp(-aa)
  a=0.d0
  if (n.eq.1) then
    s=(b-a)*func(0.5d0*(a+b))
  else
    it=3**(n-2)
    tnm=it
    del=(b-a)/(3.d0*tnm)
    ddel=del+del
    x=a+0.5d0*del
    sum=0.d0
    do j=1,it
      sum=sum+func(x)
      x=x+ddel
      sum=sum+func(x)
      x=x+del
		end do
    s=(s+(b-a)*sum/tnm)/3.d0
  endif
  return
END
