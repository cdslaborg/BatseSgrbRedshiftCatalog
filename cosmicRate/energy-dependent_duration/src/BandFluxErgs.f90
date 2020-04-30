function bandfluxergs(emin,emax,alpha,beta,epk)

	!	Amir Shahmoradi, Monday Dec 9, 2013, 9:04, IFS, The University of Texas at Austin.
	implicit none
	double precision, intent(in) :: emin,emax,alpha,beta,epk
	double precision :: bandfluxergs
	double precision :: ebrk	! ebreak
	double precision :: coef,de,e
	double precision, parameter :: erg2kev = 6.2415d8	! 1 (erg) = erg2kev (kev)

	!	---------e_break calculation-------------------------------------------------------
		ebrk = epk*(alpha-beta)/(2.d0+alpha)
	!	 ---------coef is the multiplicant of e^beta in phi/a when e>ebrk-------------------
		coef = ebrk**(alpha-beta)*dexp(beta-alpha)
	!	---------calculation of the flux in units of kev/cm2s -------------------
	
	de=0.5d0
	e=emin
	bandfluxergs=0.d0
	
	if (ebrk<emax) then
		if (emin<ebrk) then
			do while(e<=ebrk)
				e=e+de
				bandfluxergs=bandfluxergs+e**(alpha+1.d0)/dexp(e*(2.d0+alpha)/epk)
			end do
			bandfluxergs=bandfluxergs*de
			bandfluxergs=bandfluxergs+coef*(emax**(beta+2.d0)-ebrk**(beta+2.d0))/(beta+2.d0)
		else
			bandfluxergs=bandfluxergs+coef*(emax**(beta+2.d0)-emin**(beta+2.d0))/(beta+2.d0)
		end if
	else
		do while(e<=emax)
			e=e+de
			bandfluxergs=bandfluxergs+e**(alpha+1.d0)/dexp(e*(2.d0+alpha)/epk)
		end do
		bandfluxergs=bandfluxergs*de
	end if
	bandfluxergs=bandfluxergs/erg2kev

end function bandfluxergs