! Given Log10(Epk [KeV]) of an LGRB and its bolometric (0.0001-20000 KeV) peak flux [in units of Ergs/s], Log(Pbol), 
! this function calculates the corresponding Log10(peak photon flux) in the BATSE detection energy range [50,300] KeV.

! Amir Shahmoradi, last updated on Friday 7:42 PM, Dec 13, 2013, IFS/ICMB, UT Austin

function pbolepk2p1024ph_default(logepk)

	implicit none
	double precision :: pbolepk2p1024ph_default
	double precision, intent(in) :: logepk
	if (logepk<-2.915056638230699d0) then
		PbolEpk2P1024ph_default= 4.92d0
	else if (logepk>=-2.915056638230699d0 .and. logepk<1.5d0) then
		PbolEpk2P1024ph_default= &
		5.73612d0 + 0.30936d0*logepk + 0.00456d0*logepk**2 + 0.00159d0*logepk**3 + 1.53336d-4*logepk**4 - 3.5748d-4*logepk**5
	else if (logepk>=1.5d0 .and. logepk<2.5d0) then
		PbolEpk2P1024ph_default= &
		1.91128d0 + 39.71039d0*logepk - 96.60628d0*logepk**2 + 109.24696d0*logepk**3 - &
		67.2718d0*logepk**4 + 23.40239d0*logepk**5 - 4.34544d0*logepk**6 + 0.33606d0*logepk**7
	else if (logepk>=2.5d0 .and. logepk<4.d0) then
		PbolEpk2P1024ph_default= &
		2.80206d0 + 4.56907d0*logepk - 1.92772d0*logepk**2 + 0.29381d0*logepk**3 - 0.01489d0*logepk**4
	else if (logepk>=4.d0 .and. logepk<5.4093868613659435d0) then
		PbolEpk2P1024ph_default= &
		-10.46533d0 + 26.70637d0*logepk - 14.47631d0*logepk**2 + 3.54041d0*logepk**3 - 0.40957d0*logepk**4 + 0.01831d0*logepk**5
	else if (logepk>=5.4093868613659435d0) then
		PbolEpk2P1024ph_default= 4.92d0
	end if

end function PbolEpk2P1024ph_default