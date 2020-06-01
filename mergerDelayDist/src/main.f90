program DelayedMergerRate

use Constants_mod, only: IK, RK
use Cosmology_mod, only: getLogLumDisWicMpc
use StarFormation_mod, only: getbinaryMergerRate

#if defined H06
        use StarFormation_mod, only: getLogRate => getLogRateH06
#elif defined L08
        use StarFormation_mod, only: getLogRate => getLogRateL08
#elif defined B10
        use StarFormation_mod, only: getLogRate => getLogRateB10
#elif defined M14
        use StarFormation_mod, only: getLogRate => getLogRateM14
#elif defined M17
        use StarFormation_mod, only: getLogRate => getLogRateM17
#elif defined F18
        use StarFormation_mod, only: getLogRate => getLogRateF18
#else
#error "Unknown SFR model in main.f90"
#endif

implicit none

real(RK)    , parameter     :: DELTA_ZPLUS1 = 2.5e-2_RK
real(RK)    , parameter     :: ZPLUS1_MIN = 1.03_RK
real(RK)    , parameter     :: ZPLUS1_MAX = 21._RK    ! 9.9_RK
real(RK)    , parameter     :: NZPLUS1 = (ZPLUS1_MAX-ZPLUS1_MIN) / DELTA_ZPLUS1 

integer(IK)                 :: i
real(RK)                    :: zplus1(0:NZPLUS1), logzplus1, binaryMergerRate(NZPLUS1), starFormationRate(NZPLUS1)

!integer(IK)                 :: nargin                  ! number of input arguments
!integer(IK)                 :: outputFileLen           ! length of the output file string
integer(IK)                 :: outputFileUnit           ! unit number for the opened output file
character(:), allocatable   :: outputFile

! get the command line argument: outputFile

!nargin = command_argument_count()
!if (nargin/=1) then
!    write(*,*)
!    write(*,"(A)") "Fatal Error: The number of input command-line arguments to the executable should be 1."
!    write(*,"(A)") 
!    write(*,"(A)") "    Usage: "
!    write(*,"(A)") 
!    write(*,"(A)") "        main.exe ''"
!endif

!call get_command_argument(1,length=outputFileLen)
!allocate(character(outputFileLen) :: outputFile)
!call get_command_argument(1,value=outputFile)

! open outputFile

outputFile = "mergerDelayRate"

#if defined H06
        outputFile = outputFile // "H06"
#elif defined L08
        outputFile = outputFile // "L08"
#elif defined B10
        outputFile = outputFile // "B10"
#elif defined M14
        outputFile = outputFile // "M14"
#elif defined M17
        outputFile = outputFile // "M17"
#elif defined F18
        outputFile = outputFile // "F18"
#else
#error "Unknown SFR model in main.f90"
#endif

! compute rate

zplus1(0) = ZPLUS1_MIN - DELTA_ZPLUS1
do i = 1, nzplus1
    zplus1(i) = zplus1(i-1) + DELTA_ZPLUS1
    if (mod(i,100)==0) write(*,"(*(g0,:,' '))") "redshift =", zplus1(i) - 1._RK
    binaryMergerRate(i) = getbinaryMergerRate  ( zplus1 = zplus1(i) &
                                           !, zplus1Max = 100._RK &
                                           !, nRefinement &
                                           !, maxRelativeError &
                                            , getMergerDelayTimePDF = getMergerDelayTimeDistLognormal &
                                            , getStarFormationRateDensity = getStarFormationRateDensity &
                                            )
    starFormationRate(i) = exp( getLogRate(zplus1(i), log(zplus1(i)), 2*getLogLumDisWicMpc(zplus1(i))) )
    !if (zplus1(i) <= ZPLUS1_MAX) cycle
    !exit
end do

! write output

outputFile = outputFile // ".txt"

open( newunit = outputFileUnit &
    , file = outputFile &
    , status="replace"  &
    )

write(outputFileUnit,"(*(g0,:,','))") "redshift", "binaryMergerRate", "SFR"

binaryMergerRate = binaryMergerRate / sum(binaryMergerRate)     ! normalize rate
starFormationRate = starFormationRate / sum(starFormationRate)  ! normalize rate

do i = 1, nzplus1
    write(outputFileUnit,"(*(g0,:,','))") zplus1(i)-1._RK, binaryMergerRate(i), starFormationRate(i) 
end do

contains

    pure function getStarFormationRateDensity(zplus1) result(starFormationRateDensity)
#if defined H06
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityH06
#elif defined L08
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityL08
#elif defined B10
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityB10
#elif defined M14
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityM14
#elif defined M17
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityM17
#elif defined F18
        use StarFormation_mod, only: getLogRateDensity => getLogRateDensityF18
#else
#error "Unknown SFR model in main.f90"
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: starFormationRateDensity
        starFormationRateDensity = exp( getLogRateDensity(log(zplus1)) )
    end function getStarFormationRateDensity

    pure function getMergerDelayTimeDistLognormal(mergerDelayTime) result(mergerDelayTimeDistLognormal)
        use Statistics_mod, only: getLogProbLogNorm
        use Constants_mod, only: RK
        implicit none
        real(RK)    , parameter :: LOG_MEAN = -1._RK                ! mean of the lognormal merger delay time dist in GYrs.
        real(RK)    , parameter :: SIGMA = 1.11943638_RK            ! standard deviation of the lognormal merger delay time dist.
        real(RK)    , parameter :: INV_VARIANCE = 1._RK/SIGMA**2
        real(RK)    , parameter :: LOG_SQRT_INV_VARIANCE = log(sqrt(INV_VARIANCE))
        real(RK), intent(in)    :: mergerDelayTime
        real(RK)                :: mergerDelayTimeDistLognormal
        mergerDelayTimeDistLognormal = exp  ( getLogProbLogNorm ( logMean = LOG_MEAN &
                                                                , inverseVariance = INV_VARIANCE &
                                                                , logSqrtInverseVariance = LOG_SQRT_INV_VARIANCE &
                                                                , logPoint = log(mergerDelayTime) &
                                                                ) &
                                            )
    end function getMergerDelayTimeDistLognormal

end program DelayedMergerRate
