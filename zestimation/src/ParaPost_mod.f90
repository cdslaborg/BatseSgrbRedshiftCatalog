module ParaPost_mod

    use Constants_mod, only: IK, RK, PI

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaPost_mod"

    integer(IK) , parameter :: NSAMPLE = 5000_IK
    integer(IK) , parameter :: NVAR = 4     ! number of GRB attributes used in the world model
    integer(IK) , parameter :: NPAR = 16    ! number of world model's parameters

    ! the normalization factor of the multivariate log-normal distribution

    real(RK)    , parameter :: INV_TWOPI_POW_HALFNVAR = 1._RK / sqrt( (2._RK*PI)**NVAR )

    ! BATSE threshold parameters

    type :: Thresh_type
        real(RK) :: avg, invStdSqrt2
    end type Thresh_type

    ! Posterior sample

    type :: ParaPost_type
        ! order of variables: logLiso, logEpkz, logEiso, logT90z
        real(RK) :: logNormCoef, Avg(NVAR), InvCovMat(NVAR,NVAR)
        type(Thresh_type) :: Thresh
    end type ParaPost_type

    type (ParaPost_type), allocatable :: ParaPost(:)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! All parameters are assumed to be read in log Neper (not log10) from the input file, wherever needed.
    subroutine getParaPost(sampleFilePath)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Integration_mod, only: doQuadRombOpen, ErrorMessage, midexp
        use Constants_mod, only: IK, RK, SQRT2, SQRT2PI, CARRIAGE_RETURN, CLOCK_TICK
        use Matrix_mod, only: getInvPosDefMatSqrtDet
        use String_mod, only: String_type
        use Timer_mod, only: Timer_type

        implicit none

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@getModelIntegral()"

        character(*), intent(in)            :: sampleFilePath
        integer(IK)                         :: i, j, isample, ierr, itick, logFuncOffset, sampleFileUnit, numFuncEval
        real(RK), parameter                 :: INV_SQRT2 = 1._RK / SQRT2
        real(RK)                            :: rhoLisoEpkz, modelIntegral, logFunc, relativeError, dummy(2)
        real(RK), allocatable               :: LogFuncState(:)
        type(String_type)                   :: Record

        type(Timer_type)                    :: Timer

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        !::::                                           get the sample file contents                                            ::::
        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        open(newunit=sampleFileUnit,file=sampleFilePath,status="old")

        read(sampleFileUnit,*) ! read the header

        allocate(ParaPost(NSAMPLE))
        allocate(LogFuncState(0:NPAR))
        allocate( character(600) :: Record%value )

        write(output_unit,"(*(g0))")

        itick = 0
        call Timer%tic()

        do isample = 1, NSAMPLE

            ! parse the sample line

            read(sampleFileUnit, "(A)" ) Record%value
            Record%Parts = Record%SplitStr(trim(adjustl(Record%value)),",",Record%nPart)
            do j = 0, NPAR
                read(Record%Parts(j+1)%record,*) LogFuncState(j)
            end do

            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            !::::                              translate the sample to the model parameters                                 ::::
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            ! mean vector of the LGRB variables' distribution

            ParaPost(isample)%Avg(1:NVAR) = LogFuncState(1:NVAR)

            ! compute the standard deviations

            do i = 1,NVAR
                ParaPost(isample)%InvCovMat(i,i) = exp( LogFuncState(i+NVAR) )
            end do

            ! convert inverse Fisher-transform to correlation coefficient: LisoEpkz, LisoEiso, LisoT90z

            ParaPost(isample)%InvCovMat(2:4,1) = tanh(LogFuncState(1+2*NVAR:3+2*NVAR))

            ! convert inverse Fisher-transform to correlation coefficient: EpkzEiso, EpkzT90z

            ParaPost(isample)%InvCovMat(3:4,2) = tanh(LogFuncState(4+2*NVAR:5+2*NVAR))

            ! convert inverse Fisher-transform to correlation coefficient: EisoT90z

            ParaPost(isample)%InvCovMat(4:4,3) = tanh(LogFuncState(6+2*NVAR:6+2*NVAR))

            ! log of threshold mean

            ParaPost(isample)%Thresh%avg = LogFuncState(15)

            ! the threshold standard deviation

            ParaPost(isample)%Thresh%invStdSqrt2 = INV_SQRT2 / exp(LogFuncState(16))

            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            !::::                                      begin transformation of parapost                                     ::::
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            ! compute the covariance matrix's upper-triangle elements (excluding diagonal variance elements)

            do i = 1, NVAR
                do j = i+1,NVAR
                    ParaPost(isample)%InvCovMat(i,j)    = ParaPost(isample)%InvCovMat(j,i) &
                                                        * ParaPost(isample)%InvCovMat(i,i) &
                                                        * ParaPost(isample)%InvCovMat(j,j)
                end do
            end do

            ! now convert the standard deviations to variances

            do i = 1, NVAR
                ParaPost(isample)%InvCovMat(i,i) = ParaPost(isample)%InvCovMat(i,i)**2  ! This is now squared standard deviation
            end do

            ! get the inverse covariance matrix and determinant

            call getInvPosDefMatSqrtDet ( nd = NVAR &
                                        , MatInvMat = ParaPost(isample)%InvCovMat &
                                        , sqrtDetInvPosDefMat = dummy(1) )
            if (dummy(1)<0._RK) then
                write(*,*) "FATAL: Covariance Matrix to corresponding to the parameter sample ", isample, "not positive-definite."
                error stop
            end if

            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            !::::                                           compute model integral                                          ::::
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            block
                use Constants_mod, only: NEGINF_RK
                use BatseSgrbWorldModel_mod, only: getLogNormCoefLogPostProb
                dummy = getLogNormCoefLogPostProb   ( npar = NPAR &
                                                    , Param  = LogFuncState(1:NPAR) &
                                                    , logPostProbNotNeeded = .true. &
                                                    )
                ParaPost(isample)%logNormCoef = dummy(1)
                if (ParaPost(isample)%logNormCoef==NEGINF_RK) then
                    write(*,*) "FATAL: Covariance Matrix to corresponding to the parameter sample ", isample, "not positive-definite."
                    error stop
                end if
            end block

            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
            !::::                                     End of model integral computation                                     ::::
            !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            ! report the timing

            call Timer%toc()
            itick = itick + 1
            write(output_unit,"(*(' ',g0))", advance="no"   ) CARRIAGE_RETURN, CLOCK_TICK(itick) &
                                                        , isample, " out of  ", NSAMPLE &
                                                        , "parameter posterior samples processed in" &
                                                        , Timer%Time%total, "seconds."
            flush(output_unit)
            if (itick==4) itick = 0

        end do

        write(output_unit,"(*(g0))")

    end subroutine getParaPost

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module ParaPost_mod