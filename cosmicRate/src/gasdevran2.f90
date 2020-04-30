!	April 10, 2009, 00:27', Amir Shahmoradi, Michigan Tech University
	function gasdevran2(idum)
	integer idum
	double precision gasdevran2
!	ATTENTION--ATTENTION--ATTENTION--ATTENTION
!	ATTENTION: This code USES "ran2.f90"
!	This code returns a normally distributed deviate with zero mean and unit variance, using ran2(idum)
!	as the source of uniform deviates.
	integer iset
	double precision fac,gset,rsq,v1,v2,ran2
	save iset,gset
	data iset/0/
	if (iset .eq. 0) then
1		v1 = 2.0d0*ran2(idum) - 1.0d0
		v2 = 2.0d0*ran2(idum) - 1.0d0
		rsq = v1**2 + v2**2
		if (rsq .ge. 1.0d0 .or. rsq .eq. 0) goto 1
		fac = dsqrt(-2.0d0*dlog(rsq)/rsq)
		gset = v1*fac
		gasdevran2 = v2*fac
		iset = 1
	else
		gasdevran2 = gset
		iset = 0
	endif
	return
	end
		