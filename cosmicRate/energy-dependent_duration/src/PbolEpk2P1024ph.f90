! Amir Shahmoradi, last updated on Friday 7:42 PM, Dec 13, 2013, IFS/ICMB, UT Austin
function PbolEpk2P1024ph(logepk,logpbol)

  implicit none
  double precision :: pbolepk2p1024ph
  double precision, intent(in) :: logepk,logpbol
  if (logepk<-2.915056638230699d0) then
    PbolEpk2P1024ph=logpbol + 4.92d0
  else if (logepk>=-2.915056638230699d0 .and. logepk<1.5d0) then
    PbolEpk2P1024ph=logpbol + &
    5.73612d0 + 0.30936d0*logepk + 0.00456d0*logepk**2 + 0.00159d0*logepk**3 + 1.53336d-4*logepk**4 - 3.5748d-4*logepk**5
  else if (logepk>=1.5d0 .and. logepk<2.5d0) then
    PbolEpk2P1024ph=logpbol + &
    1.91128d0 + 39.71039d0*logepk - 96.60628d0*logepk**2 + 109.24696d0*logepk**3 - &
    67.2718d0*logepk**4 + 23.40239d0*logepk**5 - 4.34544d0*logepk**6 + 0.33606d0*logepk**7
  else if (logepk>=2.5d0 .and. logepk<4.d0) then
    PbolEpk2P1024ph=logpbol + &
    2.80206d0 + 4.56907d0*logepk - 1.92772d0*logepk**2 + 0.29381d0*logepk**3 - 0.01489d0*logepk**4
  else if (logepk>=4.d0 .and. logepk<5.4093868613659435d0) then
    PbolEpk2P1024ph=logpbol - &
    10.46533d0 + 26.70637d0*logepk - 14.47631d0*logepk**2 + 3.54041d0*logepk**3 - 0.40957d0*logepk**4 + 0.01831d0*logepk**5
  else if (logepk>=5.4093868613659435d0) then
    PbolEpk2P1024ph=logpbol + 4.92d0
  end if

end function PbolEpk2P1024ph