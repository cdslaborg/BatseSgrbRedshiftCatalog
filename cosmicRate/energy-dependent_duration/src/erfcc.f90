! Amir Shahmoradi, Tuesday Dec 11, 2013, 1:39 PM, IFS, The University of Texas at Austin.
FUNCTION erfcc(x)
  DOUBLE PRECISION erfcc,x
  DOUBLE PRECISION t,z
  z=DABS(x)
  t=1.0d0/(1.0d0+0.5d0*z)
  erfcc=t*DEXP(-z*z-1.26551223d0+t*(1.00002368d0+t*(0.37409196d0+t*&
        (0.09678418d0+t*(-0.18628806d0+t*(0.27886807d0+t*(-1.13520398d0+t*&
        (1.48851587d0+t*(-0.82215223d0+t*0.17087277d0)))))))))
  if (x.lt.0.0d0) erfcc=2.0d0-erfcc
END FUNCTION erfcc
