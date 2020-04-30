function bandfluxph(emin,emax,alpha,beta,epk)

	!	Amir Shahmoradi, Monday Dec 9, 2013, 9:04, IFS, The University of Texas at Austin.
	implicit none
	double precision, intent(in) :: emin,emax,alpha,beta,epk
	double precision :: bandfluxph
	double precision :: ebrk	! ebreak
	double precision :: coef,de,e

	!	---------e_break calculation-------------------------------------------------------
	ebrk = epk*(alpha-beta)/(2.d0+alpha)
	!	 ---------coef is the multiplicant of e^beta in phi/a when e>ebrk-------------------
	coef = ebrk**(alpha-beta)*dexp(beta-alpha)
	!	---------calculation of the flux in units of kev/cm2s -------------------
	
	de=0.5d0
	e=emin
	bandfluxph=0.d0
	
	if (ebrk<emax) then
		if (emin<ebrk) then
			do while(e<=ebrk)
				e=e+de
				bandfluxph=bandfluxph+e**alpha/dexp(e*(2.d0+alpha)/epk)
			end do
			bandfluxph=bandfluxph*de
			bandfluxph=bandfluxph+coef*(emax**(beta+1.d0)-ebrk**(beta+1.d0))/(beta+1.d0)
		else
			bandfluxph=bandfluxph+coef*(emax**(beta+1.d0)-emin**(beta+1.d0))/(beta+1.d0)
		end if
	else
		do while(e<=emax)
			e=e+de
			bandfluxph=bandfluxph+e**alpha/dexp(e*(2.d0+alpha)/epk)
		end do
		bandfluxph=bandfluxph*de
	end if

end function bandfluxph