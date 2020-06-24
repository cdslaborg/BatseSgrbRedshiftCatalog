module BatseSgrbWorldModel_mod

#if defined H06
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormH06
#elif defined L08
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormL08
#elif defined B10
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormB10
#elif defined M14
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormM14
#elif defined M17
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormM17
#elif defined F18
    use StarFormation_mod, only: getLogBinaryMergerRate => getLogBinaryMergerRateLognormF18
#else
#error "Unknown SFR model in BatseSgrbWorldModel_mod.f90"
#endif
    use Constants_mod, only: IK, RK, SPR, PI, NEGINF_RK
    use Batse_mod, only: GRB

    implicit none

#if defined ERR_ESTIMATION_ENABLED
    real(RK)                :: zone_relerr, zgrb_relerr, durz_relerr, liso_relerr, epkz_relerr
    real(RK)                :: zone_neval, zgrb_neval, durz_neval, liso_neval, epkz_neval
    integer(IK)             :: zgrb_count, durz_count, liso_count, epkz_count
#endif

    character(*), parameter :: MODULE_NAME = "@BatseSgrbWorldModel_mod"

    ! *********************************************
    ! world model parameters
    ! *********************************************

    integer(IK) , parameter :: ERFK = SPR   ! the real kind of the input value to erf()
    integer(IK) , parameter :: NVAR = 4_IK  ! number of GRB attributes used in the world model
    integer(IK) , parameter :: NPAR = 16_IK ! number of world model's parameters

    ! the normalization factor of the multivariate log-normal distribution

    real(RK)    , parameter :: SQRT_TWOPI_POW_NVAR = sqrt((2._RK*PI)**NVAR)

    ! the exponent of zone in time-dilation translation of T90 to T90z

#if defined kfacOneThird
    real(RK)    , parameter :: TIME_DILATION_EXPO = 0.666666666666667_RK
#endif

    ! the half-width of efficiency curve from 0 to 1, in units of threshold standard deviation
    real(RK)    , parameter :: THRESH_SIGNIFICANCE = 4._RK

    ! The Epk below and boave which the detector threshold becomes independent of the Epk of the events.
    real(RK)    , parameter :: INTEGRATION_LIMIT_LOGEPK_MIN = -6.712165960423344_RK
    real(RK)    , parameter :: INTEGRATION_LIMIT_LOGEPK_MAX = 12.455573549219071_RK

    ! The T90 duration below and boave which the detector threshold becomes independent of the durations of the events.
    real(RK)    , parameter :: INTEGRATION_LIMIT_LOGDUR_MIN = -3.2_RK
    real(RK)    , parameter :: INTEGRATION_LIMIT_LOGDUR_MAX = +2.5_RK

    ! *********************************************
    ! variables to be read from the input file
    ! *********************************************

    !type :: IntegrationLimit_type
    !    real(RK) :: zonemin = 1.1_RK, zonemax = 1.01e2_RK
    !    real(RK) :: logepkmin = -6.712165960423344_RK, logepkmax = 12.455573549219071_RK
    !end type IntegrationLimit_type
    !type(IntegrationLimit_type) :: IntegrationLimit

    ! integration specifications

#if defined CAL_STAGE_0
    real(RK)    :: zoneMin = 1.99_RK
    real(RK)    :: zoneMax = 2.01_RK
#elif defined CAL_STAGE_1
    real(RK)    :: zoneMin = 1.9_RK
    real(RK)    :: zoneMax = 2.1_RK
#else
    real(RK)    :: zoneMin = 1.09_RK    ! 1.1e0_RK
    real(RK)    :: zoneMax = 21.0_RK    ! 2.1e1_RK
#endif
    real(RK)    :: zoneTol = 1.e-3_RK   ! 1.e-4_RK
    real(RK)    :: durzTol = 1.e-4_RK   ! 5.e-5_RK
    real(RK)    :: lisoTol = 5.e-5_RK   ! 1.e-5_RK
    real(RK)    :: epkzTol = 1.e-6_RK   ! 5.e-6_RK
    integer(IK) :: zoneRef = 4_IK
    integer(IK) :: durzRef = 4_RK
    integer(IK) :: lisoRef = 5_RK
    integer(IK) :: epkzRef = 5_RK

    ! *********************************************
    ! other shared variables used in this module
    ! *********************************************

    type :: Attribute_type
        real(RK) :: logLiso, logEpkz, logEiso, logDurz
    end type Attribute_type
    type(Attribute_type) :: mv_Avg, mv_Std

    type :: Threshold_type
        real(RK) :: avg, invStdSqrt2, logPbolMin, logPbolMax
    end type Threshold_type
    type(Threshold_type) :: mv_Thresh

    type :: ConditionalVariable_type
        real(RK):: tilt, bias, avg, std, invStdSqrt2, invStdSqrt2pi
    end type ConditionalVariable_type
    type(ConditionalVariable_type) :: mv_LogLisoGivenLogDurz
    type(ConditionalVariable_type) :: mv_LogEpkzGivenLogDurz
    type(ConditionalVariable_type) :: mv_LogEpkzGivenLogDurzLogLiso

    ! The covariance matrix of the world model

    real(RK)        :: mv_CholeskyLowerLogNormModel(NVAR,NVAR), mv_DiagonalLogNormModel(NVAR)
    real(RK)        :: mv_InvCovMatLogNormModel(NVAR,NVAR)

    real(RK)        :: mv_logLisoLogPbolDiff                ! this is log10(4*pi*dl^2) where dl is luminosity distance in units of mpc
    real(RK)        :: mv_effectivePeakPhotonFluxCorrection ! P64ms peak photon flux correction to get the effective peak flux of Shahmoradi and Nemiroff 2015
    real(RK)        :: mv_logZone, mv_logZoneKfacCorrected, mv_logEpkzMin, mv_logEpkzMax, mv_logPbol
    real(RK)        :: mv_logLisoInvStdSqrt2, mv_logLisoInvStdSqrt2pi
    real(RK)        :: mv_logDurzInvStdSqrt2, mv_logDurzInvStdSqrt2pi
    real(RK)        :: mv_logLisoAtFullEfficiency
    real(RK)        :: mv_thresh_erfc_avg_at_z
    integer(IK)     :: mv_igrb                      ! index for referencing GRBs in the computation of logPostProb
    integer(IK)     :: mv_counter = 0_IK            ! counter counting how many times the function is called
    integer(IK)     :: mv_ierr = 0_IK               ! flag indicating whether a lack-of-convergence error has occurred in integrations.

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! Amir Shahmoradi, Monday 4:43 PM, May 25, 2020, SEIR, UTA, Arlington, TX.
    function getLogPostProb(npar,Param) result(logPostProb)
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: npar
        real(RK)    , intent(in)    :: Param(npar)
        real(RK)                    :: logPostProb
        real(RK)                    :: logNormCoefLogPostProb(2)
        logical                     :: logPostProbNotNeeded
        logPostProbNotNeeded = .false.
        logNormCoefLogPostProb = getLogNormCoefLogPostProb(npar,Param,logPostProbNotNeeded)
        logPostProb = logNormCoefLogPostProb(2)
    end function getLogPostProb

    ! Amir Shahmoradi, Sunday 12:12 AM, Dec 17, 2018, SEIR, UTA, Arlington, TX.
    function getLogNormCoefLogPostProb(npar,Param,logPostProbNotNeeded) result(logNormCoefLogPostProb)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK, RK, SQRT2, SQRT2PI
#if defined quadpackDPR
        use QuadPackDPR_mod, only: dqag
#elif defined quadpackSPR
        use QuadPackSPR_mod, only: qag
#else
        use Integration_mod, only: doQuadRombOpen, midexp
        use Integration_mod, only: doQuadRombClosed
        use Integration_mod, only: ErrorMessage
#endif
        use Matrix_mod, only: getCholeskyFactor, getInvMatFromCholFac
        use Batse_mod, only: MIN_LOGPH53_4_LOGPBOLZERO, MAX_LOGPH53_4_LOGPBOLZERO, THRESH_LOGPBOL64_CORRECTION
        implicit none
        integer(IK)                 :: i, j, ierr
        integer(IK) , intent(in)    :: npar
        real(RK)    , intent(in)    :: Param(npar)
        logical     , intent(in)    :: logPostProbNotNeeded
        real(RK)                    :: logNormCoefLogPostProb(2)
        real(RK)                    :: normFac, probGRB
        real(RK)                    :: rhoLisoEpkz
        real(RK)                    :: rhoLisoEiso
        real(RK)                    :: rhoLisoDurz
        real(RK)                    :: rhoEpkzEiso
        real(RK)                    :: rhoEpkzDurz
        real(RK)                    :: rhoEisoDurz
        real(RK)                    :: rhoEpkzLisoGivenDurz
        real(RK)                    :: modelint ! integral of the model over the redshift range given by variable zone.

        ! integration variables

        real(RK)                :: relerr
        integer(IK)             :: neval

#if defined quadpackDPR
        integer(IK), parameter  :: limit = 1000_IK
        integer(IK), parameter  :: lenw = 4_IK * limit
        integer(IK)             :: last
        integer(IK)             :: iwork(limit)
        real(RK)                :: work(lenw)
#endif

        mv_ierr = 0_IK
        mv_counter = mv_counter + 1_IK

        ! mean and standard deviations

        mv_Avg%logLiso  = Param(1)
        mv_Avg%logEpkz  = Param(2)
        mv_Avg%logEiso  = Param(3)
        mv_Avg%logDurz  = Param(4)
        mv_Std%logLiso  = exp(Param(5))
        mv_Std%logEpkz  = exp(Param(6))
        mv_Std%logEiso  = exp(Param(7))
        mv_Std%logDurz  = exp(Param(8))

        ! do inverse Fisher-transform to get the correlation coefficients

        rhoLisoEpkz = tanh(Param(9))
        rhoLisoEiso = tanh(Param(10))
        rhoLisoDurz = tanh(Param(11))
        rhoEpkzEiso = tanh(Param(12))
        rhoEpkzDurz = tanh(Param(13))
        rhoEisoDurz = tanh(Param(14))

        ! covariance matrix of the LogNormal GRB world model

        mv_CholeskyLowerLogNormModel(1,1) = mv_Std%logLiso * mv_Std%logLiso
        mv_CholeskyLowerLogNormModel(2,2) = mv_Std%logEpkz * mv_Std%logEpkz
        mv_CholeskyLowerLogNormModel(3,3) = mv_Std%logEiso * mv_Std%logEiso
        mv_CholeskyLowerLogNormModel(4,4) = mv_Std%logDurz * mv_Std%logDurz
        mv_CholeskyLowerLogNormModel(1,2) = mv_Std%logLiso * mv_Std%logEpkz * rhoLisoEpkz
        mv_CholeskyLowerLogNormModel(1,3) = mv_Std%logLiso * mv_Std%logEiso * rhoLisoEiso
        mv_CholeskyLowerLogNormModel(1,4) = mv_Std%logLiso * mv_Std%logDurz * rhoLisoDurz
        mv_CholeskyLowerLogNormModel(2,3) = mv_Std%logEpkz * mv_Std%logEiso * rhoEpkzEiso
        mv_CholeskyLowerLogNormModel(2,4) = mv_Std%logEpkz * mv_Std%logDurz * rhoEpkzDurz
        mv_CholeskyLowerLogNormModel(3,4) = mv_Std%logEiso * mv_Std%logDurz * rhoEisoDurz

        ! BATSE detection threshold

        mv_Thresh%avg           = Param(15)
        mv_Thresh%invStdSqrt2   = exp(Param(16))    ! momentarily is the threshold's standard deviation for the needs below.

        ! logPbol below which no trigger happens, and above which trigger efficiency is 100%.

        mv_Thresh%logPbolMin    = mv_Thresh%avg - THRESH_SIGNIFICANCE*mv_Thresh%invStdSqrt2 - MAX_LOGPH53_4_LOGPBOLZERO   ! equivalent to eff_min_lpb
        mv_Thresh%logPbolMax    = mv_Thresh%avg + THRESH_SIGNIFICANCE*mv_Thresh%invStdSqrt2 + THRESH_LOGPBOL64_CORRECTION ! equivalent to glb_max_lpb
        mv_Thresh%invStdSqrt2   = 1._RK / (mv_Thresh%invStdSqrt2*SQRT2)

        ! The parameter rhoLE_given_Durz is the partial correlation of Liso & Epkz conditional on rest-frame duration:

        rhoEpkzLisoGivenDurz = ( rhoLisoEpkz - rhoLisoDurz * rhoEpkzDurz ) / sqrt( (1._RK-rhoLisoDurz**2) * (1._RK-rhoEpkzDurz**2) )

        ! conditional standard deviations of logEpkz and logLiso given logDurz.

        mv_LogLisoGivenLogDurz%std = sqrt(1._RK - rhoLisoDurz**2) * mv_Std%logLiso
        mv_LogEpkzGivenLogDurz%std = sqrt(1._RK - rhoEpkzDurz**2) * mv_Std%logEpkz
        mv_LogEpkzGivenLogDurzLogLiso%std = sqrt(1._RK - rhoEpkzLisoGivenDurz**2) * mv_LogEpkzGivenLogDurz%std

        mv_LogLisoGivenLogDurz%tilt = rhoLisoDurz * mv_Std%logLiso / mv_Std%logDurz
        mv_LogEpkzGivenLogDurz%tilt = rhoEpkzDurz * mv_Std%logEpkz / mv_Std%logDurz
        mv_LogLisoGivenLogDurz%bias = mv_Avg%logLiso - mv_Avg%logDurz * mv_LogLisoGivenLogDurz%tilt
        mv_LogEpkzGivenLogDurz%bias = mv_Avg%logEpkz - mv_Avg%logDurz * mv_LogEpkzGivenLogDurz%tilt
        mv_LogEpkzGivenLogDurzLogLiso%tilt = rhoEpkzLisoGivenDurz * mv_LogEpkzGivenLogDurz%std / mv_LogLisoGivenLogDurz%std
  
        ! terms used in the conditional distribution of logEpkz given logLiso / logDurz or logLiso given logDurz.

        mv_logLisoInvStdSqrt2   = 1._RK / (mv_Std%logLiso * SQRT2)   ! scale factor in the exponent of Gaussian distribution.
        mv_logLisoInvStdSqrt2pi = 1._RK / (mv_Std%logLiso * SQRT2PI) ! normalization constant of the univariate Gaussian function.

        mv_logDurzInvStdSqrt2   = 1._RK / (mv_Std%logDurz * SQRT2)   ! scale factor in the exponent of Gaussian distribution.
        mv_logDurzInvStdSqrt2pi = 1._RK / (mv_Std%logDurz * SQRT2PI) ! normalization constant of the univariate Gaussian function.

        mv_LogLisoGivenLogDurz%invStdSqrt2     = 1._RK / (mv_LogLisoGivenLogDurz%std * SQRT2)
        mv_LogLisoGivenLogDurz%invStdSqrt2pi   = 1._RK / (mv_LogLisoGivenLogDurz%std * SQRT2PI)

        mv_LogEpkzGivenLogDurz%invStdSqrt2     = 1._RK / (mv_LogEpkzGivenLogDurz%std * SQRT2)
        mv_LogEpkzGivenLogDurz%invStdSqrt2pi   = 1._RK / (mv_LogEpkzGivenLogDurz%std * SQRT2PI)

        mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2     = 1._RK / (mv_LogEpkzGivenLogDurzLogLiso%std * SQRT2)
        mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2pi   = 1._RK / (mv_LogEpkzGivenLogDurzLogLiso%std * SQRT2PI)

        !write(output_unit,"(*(g20.13))") ((mv_CholeskyLowerLogNormModel(i,j),j=1,NVAR),new_line("A"),i=1,NVAR)
        call getCholeskyFactor(NVAR,mv_CholeskyLowerLogNormModel,mv_DiagonalLogNormModel)
        if (mv_DiagonalLogNormModel(1)<0._RK) then
            !write(output_unit,"(*(g0))") "covariance matrix not positive definite..cycling.."
            !write(output_unit,"(*(g20.13))") ((mv_CholeskyLowerLogNormModel(i,j),j=1,NVAR),new_line("A"),i=1,NVAR)
            logNormCoefLogPostProb = NEGINF_RK
            return
        end if
 
        ! (2*pi)^(NVAR/2)(=39.478417604357434)*sqrt(determinant)
        normFac = SQRT_TWOPI_POW_NVAR * product(mv_DiagonalLogNormModel)
        if (normFac<=0) then
            write(output_unit,"(*(g0))") "sqrt of covariance determinant is <=0: ", normFac
            write(output_unit,"(*(g0))") "Cholesky mv_DiagonalLogNormModel: "
            write(output_unit,"(*(g0))") mv_DiagonalLogNormModel
            write(output_unit,"(*(g0))") "mv_CholeskyLowerLogNormModel/CovarianceMatrix: "
            write(output_unit,"(*(g20.13))") ((mv_CholeskyLowerLogNormModel(i,j),j=1,NVAR),new_line("A"),i=1,NVAR)
            stop
        end if

        ! get the full Inverse covariance matrix

        mv_InvCovMatLogNormModel = getInvMatFromCholFac(NVAR,mv_CholeskyLowerLogNormModel,mv_DiagonalLogNormModel)

        ! compute the normalization factor of the world model by integrating over all GRB attributes, subject to BATSE threshold

#if defined quadpackDPR
        call    dqag( f             = getModelIntOverLogDurzGivenRedshift   &
                    , a             = zoneMin                               &
                    , b             = zoneMax                               &
                    , epsabs        = 0._RK                                 &
                    , epsrel        = zoneTol                               &
                    , key           = 1_IK                                  &
                    , result        = modelint                              &
                    , abserr        = relerr                                &
                    , neval         = neval                                 &
                    , ier           = ierr                                  &
                    , limit         = limit                                 &
                    , lenw          = lenw                                  &
                    , last          = last                                  &
                    , iwork         = iwork                                 &
                    , work          = work                                  &
                    )
        if (mv_ierr .or. ierr/=0_IK) then
            write(output_unit,"(*(g0))") "FATAL: @qag(): error occurred while computing model integral over redshift. ierr, neval = ", mv_ierr, neval
            error stop
        end if
#elif defined quadpackSPR
        call     qag( f             = getModelIntOverLogDurzGivenRedshift   &
                    , a             = zoneMin                               &
                    , b             = zoneMax                               &
                    , epsabs        = 0._RK                                 &
                    , epsrel        = zoneTol                               &
                    , key           = 1_IK                                  &
                    , result        = modelint                              &
                    , abserr        = relerr                                &
                    , neval         = neval                                 &
                    , ier           = ierr                                  &
                    )
        if (mv_ierr .or. ierr/=0_IK) then
            write(output_unit,"(*(g0))") "FATAL: @qag(): error occurred while computing model integral over redshift. ierr=", mv_ierr
            error stop
        end if
#else
#if defined ERR_ESTIMATION_ENABLED
        zgrb_count  = 0_IK
        liso_count  = 0_IK
        epkz_count  = 0_IK
        zgrb_neval  = 0._RK
        liso_neval  = 0._RK
        epkz_neval  = 0._RK
        zgrb_relerr = 0._RK
        liso_relerr = 0._RK
        epkz_relerr = 0._RK
#endif

#if defined CAL_STAGE_0 || CAL_STAGE_1 || CAL_STAGE_2 || CAL_STAGE_3 || CAL_STAGE_4 || CAL_STAGE_5
        call doQuadRombClosed   ( getFunc           = getModelIntOverLogDurzGivenRedshift   &
                                , lowerLim          = zoneMin                               &
                                , upperLim          = zoneMax                               &
                                , maxRelativeError  = zoneTol                               &
                                , nRefinement       = zoneRef                               &
                                , integral          = modelint                              &
                                , relativeError     = relerr                                &
                                , numFuncEval       = neval                                 &
                                , ierr              = ierr                                  &
                                )
#else
        call doQuadRombOpen ( getFunc           = getModelIntOverLogDurzGivenRedshift   &
                            , integrate         = midexp                                &
                            , lowerLim          = zoneMin                               &
                            , upperLim          = zoneMax                               &
                            , maxRelativeError  = zoneTol                               &
                            , nRefinement       = zoneRef                               &
                            , integral          = modelint                              &
                            , relativeError     = relerr                                &
                            , numFuncEval       = neval                                 &
                            , ierr              = ierr                                  &
                            )
#endif
        !write(*,*) "Zone: ", neval, relerr / modelint
        if (mv_ierr/=0_IK .or. ierr/=0_IK) then
            if (ierr/=0_IK) mv_ierr = ierr
            write(output_unit,"(*(g0))") ErrorMessage(mv_ierr)
            write(getErrFileUnit(),"(*(g0,:,','))") "getModelIntOverLogDurzGivenRedshift", zoneMin, zoneMax, modelint, relerr, neval, mv_counter
            logNormCoefLogPostProb = NEGINF_RK
            return
            !error stop
        end if
#if defined ERR_ESTIMATION_ENABLED
        zone_neval  = neval
        zone_relerr = abs(relerr) / modelint
#endif
#endif
        if (modelint<=0.0_RK) then
            write(output_unit,"(*(g0))") "model_integral (variable modelint in getLogNormCoefLogPostProb.f90) is non-positive: ", modelint
            write(getErrFileUnit(),"(*(g0,:,','))") "nonPositiveModelint", zoneMin, zoneMax, modelint, relerr, neval, mv_counter
            logNormCoefLogPostProb = NEGINF_RK
            return
            !error stop
        end if

        logNormCoefLogPostProb(1) = -log(modelint*normFac)
        if (logPostProbNotNeeded) return

        ! marginalize over all possible redshifts

        logNormCoefLogPostProb(2) = GRB%count * logNormCoefLogPostProb(1)
        loopLogPostProb: do mv_igrb = 1, GRB%count
            !probGRB = getProbGRB(2.5_RK)
#if defined quadpackDPR
            call    dqag( f             = getProbGRB    &
                        , a             = zoneMin       &
                        , b             = zoneMax       &
                        , epsabs        = 0._RK         &
                        , epsrel        = zoneTol       &
                        , key           = 1_IK          &
                        , result        = probGRB       &
                        , abserr        = relerr        &
                        , neval         = neval         &
                        , ier           = ierr          &
                        , limit         = limit         &
                        , lenw          = lenw          &
                        , last          = last          &
                        , iwork         = iwork         &
                        , work          = work          &
                        )
            if (mv_ierr/=0_IK .or. ierr/=0_IK) then
                write(output_unit,"(*(g0))") "FATAL: @qag(): error occurred while computing probGRB. ierr=", mv_ierr, ierr
                error stop
            end if
#elif defined quadpackSPR
            call     qag( f             = getProbGRB    &
                        , a             = zoneMin       &
                        , b             = zoneMax       &
                        , epsabs        = 0._RK         &
                        , epsrel        = zoneTol       &
                        , key           = 1_IK          &
                        , result        = probGRB       &
                        , abserr        = relerr        &
                        , neval         = neval         &
                        , ier           = ierr          &
                        )
            if (mv_ierr/=0_IK .or. ierr/=0_IK) then
                write(output_unit,"(*(g0))") "FATAL: @qag(): error occurred while computing probGRB. ierr=", mv_ierr, ierr
                error stop
            end if
#else

#if defined CAL_STAGE_0 || CAL_STAGE_1 || CAL_STAGE_2 || CAL_STAGE_3 || CAL_STAGE_4 || CAL_STAGE_5
            call doQuadRombClosed   ( getFunc           = getProbGRB    &
                                    , lowerLim          = zoneMin       &
                                    , upperLim          = zoneMax       &
                                    , maxRelativeError  = zoneTol       &
                                    , nRefinement       = zoneRef       &
                                    , integral          = probGRB       &
                                    , relativeError     = relerr        &
                                    , numFuncEval       = neval         &
                                    , ierr              = ierr          &
                                    )
#else
            call doQuadRombOpen ( getFunc           = getProbGRB    &
                                , integrate         = midexp        &
                                , lowerLim          = zoneMin       &
                                , upperLim          = zoneMax       &
                                , maxRelativeError  = zoneTol       &
                                , nRefinement       = zoneRef       &
                                , integral          = probGRB       &
                                , relativeError     = relerr        &
                                , numFuncEval       = neval         &
                                , ierr              = ierr          &
                                )
#endif
            !write(*,*) "Zone, ith GRB: ", mv_igrb, neval, relerr / probGRB
            if (mv_ierr/=0_IK .or. ierr/=0_IK) then
                if (ierr/=0_IK) mv_ierr = ierr
                write(output_unit,"(*(g0))") ErrorMessage(mv_ierr)
                write(getErrFileUnit(),"(*(g0,:,','))") "getProbGRB", zoneMin, zoneMax, probGRB, relerr, neval, mv_counter
                logNormCoefLogPostProb(2) = NEGINF_RK
                return
                !error stop
            end if
#if defined ERR_ESTIMATION_ENABLED
        zgrb_count  = zgrb_count + 1_IK
        zgrb_neval  = zgrb_neval + neval
        zgrb_relerr = zgrb_relerr + abs(relerr) / probGRB
#endif
#endif
            if (probGRB<=0.0_RK) then
                !write(output_unit,"(*(g0))") "WARNING: probGRB <= 0.0_RK: ", probGRB, ". Setting logNormCoefLogPostProb(2) = NEGINF_RK ..."
                write(getErrFileUnit(),"(*(g0,:,','))") "nonPositiveProbGRB", zoneMin, zoneMax, probGRB, relerr, neval, mv_counter
                logNormCoefLogPostProb(2) = NEGINF_RK
                return
                !exit loopLogPostProb
            end if
            logNormCoefLogPostProb(2) = logNormCoefLogPostProb(2) + log(probGRB)
        end do loopLogPostProb

#if defined ERR_ESTIMATION_ENABLED
        zgrb_neval = zgrb_neval / zgrb_count
        liso_neval = liso_neval / liso_count
        epkz_neval = epkz_neval / epkz_count
        zgrb_relerr = zgrb_relerr / zgrb_count
        liso_relerr = liso_relerr / liso_count
        epkz_relerr = epkz_relerr / epkz_count
#endif

    end function getLogNormCoefLogPostProb

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! integral of grb world model at given redshift z, detemined by zone.
    function getModelIntOverLogDurzGivenRedshift(zone) result(modelIntOverLogDurzGivenRedshift)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Cosmology_mod, only: LOGMPC2CMSQ4PI, getLogLumDisWicMpc
       !use IntegrationOverLiso_mod, only: doQuadRombClosed, ErrorMessage
        use Integration_mod, only: doQuadRombClosed, ErrorMessage
        !use StarFormation_mod, only: getBinaryMergerRateS15
        use Constants_mod, only: RK, SPR
        use Batse_mod, only: THRESH_ERFC_AVG
        implicit none

        character(*), parameter :: PROCEDURE_NAME = "@getModelIntOverLogDurzGivenRedshift()"
        real(RK), intent(in)    :: zone
        real(RK)                :: relerr 
        real(RK)                :: modelIntOverLogDurzGivenRedshift
        integer(IK)             :: neval, ierr

        if (mv_ierr/=0_IK) then
            modelIntOverLogDurzGivenRedshift = NEGINF_RK
            return
        end if

        mv_logZone = log(zone)
#if defined kfacOneThird
        mv_logZoneKfacCorrected = mv_logZone * TIME_DILATION_EXPO
#elif defined kfacNone
        mv_logZoneKfacCorrected = mv_logZone
#else
#error "Unknown kfactor in BatseSgrbWorldModel_mod.f90"
#endif
        mv_logLisoLogPbolDiff = LOGMPC2CMSQ4PI + 2_IK * getLogLumDisWicMpc(zone) ! This is later used in 

        ! These are used in getModelIntOverLogLisoGivenRedshiftDurz()

        mv_thresh_erfc_avg_at_z = THRESH_ERFC_AVG - mv_logZoneKfacCorrected

        ! These are used in getModelIntOverLogLisoGivenDurz()

        mv_logEpkzMin = INTEGRATION_LIMIT_LOGEPK_MIN + mv_logZone
        mv_logEpkzMax = INTEGRATION_LIMIT_LOGEPK_MAX + mv_logZone

        ! compute world model integral over logDurz in the varying BATSE efficiency range, for the given z.

        !mv_logLisoAtFullEfficiency = mv_Thresh%logPbolMax + mv_logLisoLogPbolDiff
        call doQuadRombClosed   ( getFunc           = getModelIntOverLogLisoGivenRedshiftDurz                   &
                                , lowerLim          = INTEGRATION_LIMIT_LOGDUR_MIN - mv_logZoneKfacCorrected    &
                                , upperLim          = INTEGRATION_LIMIT_LOGDUR_MAX - mv_logZoneKfacCorrected    &
                                , maxRelativeError  = durzTol                                                   &
                                , nRefinement       = durzRef                                                   &
                                , integral          = modelIntOverLogDurzGivenRedshift                          &
                                , relativeError     = relerr                                                    &
                                , numFuncEval       = neval                                                     &
                                , ierr              = ierr                                                      &
                                )
!write(*,*) modelIntOverLogDurzGivenRedshift
#if defined ERR_ESTIMATION_ENABLED
        durz_count  = durz_count + 1_IK
        durz_neval  = durz_neval + neval
        durz_relerr = durz_relerr + abs(relerr) / modelIntOverLogDurzGivenRedshift
#endif
        if (ierr/=0_IK .or. mv_ierr/=0_IK) then
            if (ierr/=0_IK) mv_ierr = ierr
            write(output_unit,"(*(g0))") PROCEDURE_NAME // ErrorMessage(mv_ierr)
            write(getErrFileUnit(),"(*(g0,:,','))"  ) "getModelIntOverLogLisoGivenRedshiftDurz" &
                                                    , INTEGRATION_LIMIT_LOGDUR_MIN - mv_logZoneKfacCorrected &
                                                    , INTEGRATION_LIMIT_LOGDUR_MAX - mv_logZoneKfacCorrected &
                                                    , modelIntOverLogDurzGivenRedshift &
                                                    , relerr, neval, mv_counter
            modelIntOverLogDurzGivenRedshift = NEGINF_RK
            return
            !error stop
        end if
        ! add the analytical integral of the logLiso range within which BATSE efficiency is 100%
        ! then, multiply the integral result by the GRB rate density at the given redshift

        modelIntOverLogDurzGivenRedshift    = &
                                            ( modelIntOverLogDurzGivenRedshift &
                                            + 0.5_RK * erfc( real( (mv_Thresh%logPbolMax+mv_logLisoLogPbolDiff-mv_Avg%logLiso)*mv_logLisoInvStdSqrt2, kind=ERFK) ) &
                                            ) &
                                            * exp(getLogBinaryMergerRate(mv_logZone))
                                            !* getBinaryMergerRateS15(zone-1._RK)

    end function getModelIntOverLogDurzGivenRedshift

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! integral of grb world model at given redshift z, determined by zone, and durz. equivalent to probatdurz
    function getModelIntOverLogLisoGivenRedshiftDurz(logDurz) result(modelIntOverLogLisoGivenRedshiftDurz)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Cosmology_mod, only: LOGMPC2CMSQ4PI, getLogLumDisWicMpc
       !use IntegrationOverLiso_mod, only: doQuadRombClosed, ErrorMessage
        use Integration_mod, only: doQuadRombClosed, ErrorMessage
        use Constants_mod, only: RK, SPR
        use Batse_mod, only: THRESH_ERFC_AMP, THRESH_ERFC_STD_INV
        implicit none

        character(*), parameter :: PROCEDURE_NAME = "@getModelIntOverLogLisoGivenRedshiftDurz()"
        real(RK), intent(in)    :: logDurz
        real(RK)                :: relerr
        real(RK)                :: minLogPbolGivenDur
        real(RK)                :: modelIntOverLogLisoGivenRedshiftDurz
        integer(IK)             :: neval, ierr

        if (mv_ierr/=0_IK) then
            modelIntOverLogLisoGivenRedshiftDurz = NEGINF_RK
            return
        end if

        mv_effectivePeakPhotonFluxCorrection = THRESH_ERFC_AMP * erfc(real((logDurz-mv_thresh_erfc_avg_at_z)*THRESH_ERFC_STD_INV,kind=ERFK)) 
        minLogPbolGivenDur = mv_Thresh%logPbolMin + mv_effectivePeakPhotonFluxCorrection

        mv_LogLisoGivenLogDurz%avg = mv_LogLisoGivenLogDurz%bias + mv_LogLisoGivenLogDurz%tilt * logDurz
        mv_LogEpkzGivenLogDurz%avg = mv_LogEpkzGivenLogDurz%bias + mv_LogEpkzGivenLogDurz%tilt * logDurz
        mv_LogEpkzGivenLogDurzLogLiso%bias = mv_LogEpkzGivenLogDurz%avg - mv_LogLisoGivenLogDurz%avg * mv_LogEpkzGivenLogDurzLogLiso%tilt

        ! compute world model integral over logLiso in the varying BATSE efficiency range, for the given z.
        call doQuadRombClosed   ( getFunc           = getModelIntOverLogEpkzGivenRedshiftDurzLiso   &
                                , lowerLim          = minLogPbolGivenDur + mv_logLisoLogPbolDiff    &
                                , upperLim          = mv_Thresh%logPbolMax + mv_logLisoLogPbolDiff  &
                                , maxRelativeError  = lisoTol                                       &
                                , nRefinement       = lisoRef                                       &
                                , integral          = modelIntOverLogLisoGivenRedshiftDurz          &
                                , relativeError     = relerr                                        &
                                , numFuncEval       = neval                                         &
                                , ierr              = ierr                                          &
                                )
!write(*,*) modelIntOverLogLisoGivenRedshiftDurz
#if defined ERR_ESTIMATION_ENABLED
        liso_count  = liso_count + 1_IK
        liso_neval  = liso_neval + neval
        liso_relerr = liso_relerr + abs(relerr) / modelIntOverLogLisoGivenRedshiftDurz
#endif
        if (ierr/=0_IK .or. mv_ierr/=0_IK) then
            if (ierr/=0_IK) mv_ierr = ierr
            write(output_unit,"(*(g0))") PROCEDURE_NAME // ErrorMessage(mv_ierr)
            write(getErrFileUnit(),"(*(g0,:,','))"  ) "getModelIntOverLogEpkzGivenRedshiftDurzLiso" &
                                                    , mv_Thresh%logPbolMin + mv_logLisoLogPbolDiff &
                                                    , mv_logLisoAtFullEfficiency &
                                                    , modelIntOverLogLisoGivenRedshiftDurz &
                                                    , relerr, neval, mv_counter
            modelIntOverLogLisoGivenRedshiftDurz = NEGINF_RK
            return
            !error stop
        end if

        ! add the analytical integral of the logLiso range within which BATSE efficiency is 100%
        modelIntOverLogLisoGivenRedshiftDurz    = modelIntOverLogLisoGivenRedshiftDurz * mv_logDurzInvStdSqrt2pi &
                                                * exp( -( (logDurz-mv_Avg%logDurz) * mv_logDurzInvStdSqrt2 )**2 )

    end function getModelIntOverLogLisoGivenRedshiftDurz

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function getModelIntOverLogEpkzGivenRedshiftDurzLiso(logLiso) result(modelIntOverLogEpkzGivenRedshiftDurzLiso)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Batse_mod, only: MIN_LOGPH53_4_LOGPBOLZERO
        use Integration_mod, only: doQuadRombClosed, ErrorMessage
       !use IntegrationOverEpkz_mod, only: doQuadRombClosed, ErrorMessage
        use Constants_mod, only: RK

        implicit none
        character(*), parameter :: PROCEDURE_NAME = "@getModelIntOverLogEpkzGivenRedshiftDurzLiso()"
        real(RK), intent(in)    :: logLiso
        real(RK)                :: modelIntOverLogEpkzGivenRedshiftDurzLiso
        real(RK)                :: relerr
        integer(IK)             :: neval

        if (mv_ierr/=0_IK) then
            modelIntOverLogEpkzGivenRedshiftDurzLiso = NEGINF_RK
            return
        end if

        mv_LogEpkzGivenLogDurzLogLiso%avg = mv_LogEpkzGivenLogDurzLogLiso%bias + mv_LogEpkzGivenLogDurzLogLiso%tilt * logLiso
        mv_logPbol = logLiso - mv_logLisoLogPbolDiff

        call doQuadRombClosed   ( getFunc           =  getProbEpkzGivenRedshiftDurzLiso         &
                                , lowerLim          =  mv_logEpkzMin                            &
                                , upperLim          =  mv_logEpkzMax                            &
                                , maxRelativeError  =  epkzTol                                  &
                                , nRefinement       =  epkzRef                                  &
                                , integral          =  modelIntOverLogEpkzGivenRedshiftDurzLiso &
                                , relativeError     =  relerr                                   &
                                , numFuncEval       =  neval                                    &
                                , ierr              =  mv_ierr                                  &
                                )
!write(*,*) modelIntOverLogEpkzGivenRedshiftDurzLiso
#if defined ERR_ESTIMATION_ENABLED
        epkz_count  = epkz_count + 1_IK
        epkz_neval  = epkz_neval + neval
        epkz_relerr = epkz_relerr + abs(relerr) / modelIntOverLogEpkzGivenRedshiftDurzLiso
#endif
        if (mv_ierr/=0_IK) then
            write(output_unit,"(*(g0))") PROCEDURE_NAME // ErrorMessage(mv_ierr)
            write(getErrFileUnit(),"(*(g0,:,','))"  ) "getProbEpkzGivenRedshiftDurzLiso" &
                                                    , mv_logEpkzMin &
                                                    , mv_logEpkzMax &
                                                    , modelIntOverLogEpkzGivenRedshiftDurzLiso &
                                                    , relerr, neval, mv_counter
            modelIntOverLogEpkzGivenRedshiftDurzLiso = NEGINF_RK
            return
            !error stop
        end if

        ! add integral of the tails of the conditional logEpkz distribution given mv_logLiso

        modelIntOverLogEpkzGivenRedshiftDurzLiso    = modelIntOverLogEpkzGivenRedshiftDurzLiso * mv_LogLisoGivenLogDurz%invStdSqrt2pi &
                                                    * exp( -( (logLiso-mv_LogLisoGivenLogDurz%avg) * mv_LogLisoGivenLogDurz%invStdSqrt2 )**2 )

    end function getModelIntOverLogEpkzGivenRedshiftDurzLiso

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    !pure function getProbEpkzGivenRedshiftDurzLiso(logEpkz) result(probEpkzGivenRedshiftDurzLiso)
    function getProbEpkzGivenRedshiftDurzLiso(logEpkz) result(probEpkzGivenRedshiftDurzLiso)
        use Constants_mod, only: RK
        use Batse_mod, only: getLogPF53
        implicit none
        real(RK), intent(in)    :: logEpkz
        real(RK)                :: probEpkzGivenRedshiftDurzLiso, efficiency
        real(ERFK)              :: normedLogPF53
        normedLogPF53 = ( getLogPF53(logEpkz-mv_logZone,mv_logPbol) - mv_effectivePeakPhotonFluxCorrection - mv_Thresh%avg) * mv_Thresh%invStdSqrt2
        efficiency = 0.5_RK + 0.5_RK * erf(normedLogPF53)
        probEpkzGivenRedshiftDurzLiso   = efficiency * mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2pi &
                                        * exp( -( (logEpkz-mv_LogEpkzGivenLogDurzLogLiso%avg)*mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2)**2 )
!if (probEpkzGivenRedshiftDurzLiso<=0._RK) then
!    write(*,*) "normedLogPF53 = ", normedLogPF53
!    write(*,*) "efficiency = ", efficiency
!    write(*,*) "(logEpkz-mv_LogEpkzGivenLogDurzLogLiso%avg)*mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2)**2 = ", ( (logEpkz-mv_LogEpkzGivenLogDurzLogLiso%avg)*mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2)**2
!    write(*,*) "exp( -( (logEpkz-mv_LogEpkzGivenLogDurzLogLiso%avg)*mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2)**2 ) = ", exp( -( (logEpkz-mv_LogEpkzGivenLogDurzLogLiso%avg)*mv_LogEpkzGivenLogDurzLogLiso%invStdSqrt2)**2 )
!    write(*,*) "probEpkzGivenRedshiftDurzLiso = ", probEpkzGivenRedshiftDurzLiso
!    error stop
!end if
    end function getProbEpkzGivenRedshiftDurzLiso

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function getProbGRB(zone) result(probGRB)

        !use StarFormation_mod, only: getBinaryMergerRateS15
        use Cosmology_mod, only: LOGMPC2CMSQ4PI, getLogLumDisWicMpc
        use Constants_mod, only: RK
        use Batse_mod, only: getLogPF53
        implicit none
        real(RK), intent(in)    :: zone
        real(RK)                :: probGRB
        real(RK)                :: MeanSubtractedVar(NVAR)
        real(RK)                :: logZone, logLisoLogPbolDiff
        real(ERFK)              :: normedLogPF53

        logZone = log(zone)
        logLisoLogPbolDiff = LOGMPC2CMSQ4PI + 2_IK * getLogLumDisWicMpc(zone) ! log(4*pi*dl^2) where dl is luminosity distance in units of mpc

        ! observed data probability

        MeanSubtractedVar(1) = GRB%Event(mv_igrb)%logPbol - mv_Avg%logLiso + logLisoLogPbolDiff
        MeanSubtractedVar(2) = GRB%Event(mv_igrb)%logEpk  - mv_Avg%logEpkz + logZone
        MeanSubtractedVar(3) = GRB%Event(mv_igrb)%logsbol - mv_Avg%logEiso - logZone + logLisoLogPbolDiff
#if defined kfacOneThird
        MeanSubtractedVar(4) = GRB%Event(mv_igrb)%logt90  - mv_Avg%logDurz - logZone * TIME_DILATION_EXPO
#elif defined kfacNone
        MeanSubtractedVar(4) = GRB%Event(mv_igrb)%logt90  - mv_Avg%logDurz - logZone
#else
#error "kfactor model requested in BatseSgrbWorldModel_mod.f90"
#endif

        normedLogPF53 = mv_Thresh%invStdSqrt2 * ( GRB%Event(mv_igrb)%logPF53 - mv_Thresh%avg )
        

        probGRB = (0.5_RK + 0.5_RK * erf(normedLogPF53))    &   ! BATSE efficiency
                * exp( getLogBinaryMergerRate(logZone) - 0.5_RK * dot_product( MeanSubtractedVar , matmul(mv_InvCovMatLogNormModel,MeanSubtractedVar) ) )

    end function getProbGRB

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getBatseEfficiency(normedLogPF53) result(batseEfficiency)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in) :: normedLogPF53
        real(RK)             :: batseEfficiency
        batseEfficiency = 0.5_RK + 0.5_RK * erf( real( normedLogPF53 , kind=ERFK ) )
    end function getBatseEfficiency

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getBatseEfficiencyApprox(logEpk,logPbol) result(batseEfficiency)
        use Constants_mod, only: RK
        use Batse_mod, only: getLogPF53
        implicit none
        real(RK), intent(in) :: logEpk,logPbol
        real(RK)             :: batseEfficiency
        real(RK)             :: normedLogPF53
        if ( logPbol < mv_Thresh%logPbolMin ) then
            batseEfficiency = 0._RK
        elseif ( logPbol < mv_Thresh%logPbolMax ) then
            normedLogPF53 = ( getLogPF53(logEpk,logPbol) - mv_Thresh%avg ) * mv_Thresh%invStdSqrt2
            batseEfficiency = 0.5_RK + 0.5_RK * erf( real( normedLogPF53 , kind=ERFK ) )
        elseif ( logPbol >= mv_Thresh%logPbolMax ) then
            batseEfficiency = 1._RK
        end if
    end function getBatseEfficiencyApprox

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! create the error-catching report file
    function getErrFileUnit() result(errFileUnit)
        use String_mod, only: num2str
        implicit none
        integer(IK) :: errFileUnit
        integer(IK), save :: imageID = 1_IK, imageCount = 1_IK, dummyFileUnit = 1_IK
        if (dummyFileUnit>0_IK) then
#if defined MPI_ENABLED
            block
                use mpi
                integer(IK) :: ierrMPI
                logical     :: isInitialized
                call mpi_initialized( isInitialized, ierrMPI )
                if (.not. isInitialized) call mpi_init(ierrMPI)
                call mpi_comm_rank(mpi_comm_world, imageID, ierrMPI)
                call mpi_comm_size(mpi_comm_world, imageCount, ierrMPI)
                imageID = imageID + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
            end block
#endif
            open(newunit=dummyFileUnit,file="divergenceErrorReport_"//num2str(imageID)//".txt",status="replace")
            write(dummyFileUnit,"(*(g0,:,','))" ) "errorLocation", "integrationLowerLimit", "integrationUpperLimit", "integrationResult", "relerr", "neval", "MCMCStep"
        end if
        errFileUnit = dummyFileUnit
    end function getErrFileUnit

!***********************************************************************************************************************************
!***********************************************************************************************************************************

        !close(errFileUnit)
!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module BatseSgrbWorldModel_mod