! Author: Amir Shahmoradi, last updated on Sunday 10:51 PM, Dec 22, 2013, IFS/ICMB, UT Austin
! Correlation matrix is now sampled from its Fisher transformation.

module BatseSgrbWorldModel_mod

#ifdef H06
    use StarFormation_mod, only: getLogSFR => getLogSFRH06
#elif defined L08
    use StarFormation_mod, only: getLogSFR => getLogSFRL08
#elif defined B10
    use StarFormation_mod, only: getLogSFR => getLogSFRB10
#elif defined M14
    use StarFormation_mod, only: getLogSFR => getLogSFRM14
#else
#error "Unknown SFR model in BatseSgrbWorldModel_mod.f90"
#endif
    use Constants_mod, only: IK, RK, PI, LN10
    use Batse_mod, only: GRB

    implicit none

    character(*), parameter :: MODULE_NAME = "@BatseSgrbWorldModel_mod"

    ! ******************************************************************************************************************************
    ! world model parameters
    ! ******************************************************************************************************************************

    integer(IK), parameter :: NVAR = 4     ! number of GRB attributes used in the world model
    integer(IK), parameter :: NPAR = 16    ! number of world model's parameters

    ! the normalization factor of the multivariate log-normal distribution
    real(RK), parameter :: SQRT_TWOPI_POW_NVAR = sqrt((2._RK*PI)**NVAR)

    ! the exponent of zone in time-dilation translation of T90 to T90z
#ifdef KFAC_ONETHIRD_ENABLED
    real(RK), parameter :: TIME_DILATION_EXPO = 0.666666666666667_RK
#endif

    ! the half-width of efficiency curve from 0 to 1, in units of threshold standard deviation
    real(RK), parameter :: THRESH_HALF_WIDTH = 4._RK

    real(RK), parameter :: INTEGRATION_LIMIT_LOGEPK_MIN = -6.712165960423344_RK
    real(RK), parameter :: INTEGRATION_LIMIT_LOGEPK_MAX = 12.455573549219071_RK

    ! observer-frame durations above and below which the average threshold at the given duration becomes practically independent of the duration.
    real(RK), parameter :: INTEGRATION_LIMIT_LOGDUR_MIN = LN10 * -3.2_RK
    real(RK), parameter :: INTEGRATION_LIMIT_LOGDUR_MAX = LN10 * +2.5_RK

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function getLogPostProb(PARAM)    ! Stands for Gaussian Mixture Log(Likelihood)
    use MODELparameters
    use Zparameters, only: zmin,zmax    !,nz
    use constants, only: sqrt2,sqrt2Pi    ! ,sqrtPiO2
    use GRBworld, only: npar
    use OBSGRBDATA
    use detection
    implicit none
    INTEGER :: i,j,posdef
    real(RK) :: getLogPostProb,determinant
    real(RK) :: lisosigma,epkzsigma,normfac,t90zsigma,dumvar
    real(RK), DIMENSION(npar) :: PARAM
    !real(RK), EXTERNAL :: probatz,modelintz
    !real(RK) :: z,modelintz,probatz        ! Here it is assumed that z=0, therefore, modelintz=modelint
    real(RK) :: probDetectionGRB 
    
    real(RK) :: time_integration_end,time_integration_start 
    
    ! LGRB model parameters:
        loglisomean=PARAM(1)
        logepkzmean=PARAM(2)
        logeisomean=PARAM(3)
        logt90zmean=PARAM(4)
        lisosigma=10.d0**PARAM(5)
        epkzsigma=10.d0**PARAM(6)
        ! The conditional meand & stdev of Epkz distribution on the conditional_distribution_of_Liso_on_duration will be calculated in function ProbatLisoGivenDurz.f90
        t90zsigma=10.d0**PARAM(8)
        ! The correlation Coefficients array (apply inverse Fisher transformation):
        do i=9,14
            dumvar=2.d0*PARAM(i)
            RHO(i)=(dexp(dumvar)-1.d0)/(dexp(dumvar)+1.d0)
        end do
    ! Other parameters needed for the purpose of model integrations:
        meanthresh=PARAM(15)
        stdevthresh=10.d0**PARAM(16)
        eff_min_pph=meanthresh-signif*stdevthresh    ! The EFFective Log(P1024ph) limit below which no trigger happens.
                                                    ! It is equivalent to minimum Log(P1024ph) at very long durations.
        eff_max_pph=meanthresh+signif*stdevthresh    ! The EFFective Log(P1024ph) limit above which trigger efficiency is 100%.
                                                    ! It is equivalent to maximum Log(P1024ph) at very long durations.
        eff_min_lpb=eff_min_pph-logphminmaxdiff-logp1024min    ! Effective LogPbol limit below which no trigger happens, for any Log(Epk).
                                                            ! It is equivalent to minimum Log(Pbol) at very long durations.
        !eff_max_lpb=eff_max_pph+logphminmaxdiff-logp1024min    ! Effective LogPbol limit above which trigger efficiency is 100%, for any Log(Epk).
                                                            ! It is equivalent to maximum Log(Pbol) at very long durations.
        glb_max_lpb=eff_max_pph+lpb_correction
        ! The parameter rhoLE_given_Durz is the partial correlation of Liso & Epkz conditional on rest-frame duration:
        ! It is used in the integration of the model where the joint distribution of Lis-Epkz given rest-frame duration
        ! has to be integrated. The relation can be found at http://en.wikipedia.org/wiki/Partial_correlation
        ! LGRBs:
            rhoLE_given_Durz=(RHO(9)-RHO(11)*RHO(13))/&
            dsqrt((1.d0-RHO(11)**2)*(1.d0-RHO(13)**2))
        ! SGRBs:
            ! ATT: Change to RHO required: rhoLE_given_Durz(2)=(PARAM(24)-PARAM(26)*PARAM(28))/&
            !dsqrt((1.d0-PARAM(26)**2)*(1.d0-PARAM(28)**2))
        conlisosigmaD=dsqrt(1.d0-RHO(11)**2)*lisosigma    ! conditional stdev of LGRB Liso distribution on rest-frame Duration variable. needed for model integration
        !conlisosigmaD(2)=dsqrt(1.d0-PARAM(26)**2)*lisosigma(2)    ! conditional stdev of LGRB Liso distribution on rest-frame Duration variable. needed for model integration
        ! The two parameters below are needed for the calculation of bELD parameters:
            conepkzsigmaD=dsqrt(1.d0-RHO(13)**2)*epkzsigma    ! conditional stdev of LGRB Epkz distribution on rest-frame Duration variable.
            !conepkzsigmaD(2)=dsqrt(1.d0-PARAM(28)**2)*epkzsigma(2)    ! conditional stdev of LGRB Epkz distribution on rest-frame Duration variable.
        ! The two parameters below is the conditional stdev of Epkz distribution on the conditional dist. of Liso distribution on Durationz.
            conepkzsigmaLD=dsqrt(1.d0-rhoLE_given_Durz**2)*conepkzsigmaD
            !conepkzsigmaLD(2)=dsqrt(1.d0-rhoLE_given_Durz(2)**2)*conepkzsigmaD(2)
        ! The means of the conditional distributions of Liso given durz will be calculated in function probatdurz.f90
        bLD=RHO(11)*lisosigma/t90zsigma    ! term used in the coditional distribution of Liso given durationZ.
        !bLD(2)=PARAM(26)*lisosigma(2)/t90zsigma(2)    ! term used in the coditional distribution of Liso given durationZ.
        aLD=loglisomean-logt90zmean*bLD    ! term used in the coditional distribution of Liso given durationz.
        !aLD(2)=loglisomean(2)-logt90zmean(2)*bLD(2)    ! term used in the coditional distribution of Liso given durationz.
        bED=RHO(13)*epkzsigma/t90zsigma    ! term used in the coditional distribution of Epkz given durationZ.
        !bED(2)=PARAM(28)*epkzsigma(2)/t90zsigma(2)    ! term used in the coditional distribution of Epkz given durationZ.
        aED=logepkzmean-logt90zmean*bED    ! term used in the coditional distribution of Epkz given durationz.
        !aED(2)=logepkzmean(2)-logt90zmean(2)*bED(2)    ! term used in the coditional distribution of Epkz given durationz.
        bELD=rhoLE_given_Durz*conepkzsigmaD/conlisosigmaD    ! term used in the coditional distribution of Epkz given Luminosity.
        !bELD(2)=rhoLE_given_Durz(2)*conepkzsigmaD(2)/conlisosigmaD(2)    ! term used in the coditional distribution of Epkz given Luminosity.
        ! aELD parameters will be calculated in Lisointergation.f90 since it requires the conditional mean of Liso distribution on durationz.
        sqrt2lisosigma=sqrt2*lisosigma                ! defined for code efficiency, used in model integration.
        sqrt2t90zsigma=sqrt2*t90zsigma                ! defined for code efficiency, used in model integration.
        sqrt2Pit90zsigma=sqrt2Pi*t90zsigma                ! defined for code efficiency, used in model integration.
        sqrt2conlisosigmaD=sqrt2*conlisosigmaD        ! This is the square root of the scale factor in the exponent of the Gaussian dist. function for LGRBs.
        sqrt2conepkzsigmaLD=sqrt2*conepkzsigmaLD        ! This is the square root of the scale factor in the exponent of the Gaussian dist. function for LGRBs.
        sqrt2PiconlisosigmaD=sqrt2Pi*conlisosigmaD    ! This is the normalization constant of the univariate Gaussian function for LGRBs.
        sqrt2PiconepkzsigmaLD=sqrt2Pi*conepkzsigmaLD    ! This is the normalization constant of the univariate Gaussian function for LGRBs.
        sqrt2stdevthresh=sqrt2*stdevthresh
    ! SGRBs: Construct the covariance matrix of the SGRBs Gaussian model:
        do i=1,nvar
            ! expression below corresponds to lisosigma,epkzsigma,eisosigma for LGRB model respectively.
                MVNCOV(i,i)=10.d0**PARAM(i+4)
        end do
        ! Expressions below correspond to corr(Liso,Epkz) = correlation of Log(Liso) and Log(Epkz)
            MVNCOV(1,2)=MVNCOV(1,1)*MVNCOV(2,2)*RHO(9)
            MVNCOV(2,1)=MVNCOV(1,2)
        ! Expressions below correspond to corr(Liso,Eiso) = correlation of Log(Liso) and Log(Eiso).
            MVNCOV(1,3)=MVNCOV(1,1)*MVNCOV(3,3)*RHO(10)
            MVNCOV(3,1)=MVNCOV(1,3)
        ! Expressions below correspond to corr(Liso,T90z) = correlation of Log(Liso) and Log(Epkz)
            MVNCOV(1,4)=MVNCOV(1,1)*MVNCOV(4,4)*RHO(11)
            MVNCOV(4,1)=MVNCOV(1,4)
        ! Expressions below correspond to corr(Epkz,Eiso) = correlation of Log(Epkz) and Log(Eiso)
            MVNCOV(2,3)=MVNCOV(2,2)*MVNCOV(3,3)*RHO(12)
            MVNCOV(3,2)=MVNCOV(2,3)
        ! Expressions below correspond to corr(Epkz,T90z) = correlation of Log(Epkz) and Log(T90z)
            MVNCOV(2,4)=MVNCOV(2,2)*MVNCOV(4,4)*RHO(13)
            MVNCOV(4,2)=MVNCOV(2,4)
        ! Expressions below correspond to corr(Eiso,T90z) = correlation of Log(Eiso) and Log(T90z)
            MVNCOV(3,4)=MVNCOV(3,3)*MVNCOV(4,4)*RHO(14)
            MVNCOV(4,3)=MVNCOV(3,4)
        do i=1,nvar
            MVNCOV(i,i)=MVNCOV(i,i)**2
            ! This corresponds to lisosigma^2,epkzsigma^2,eisosigma^2,t90zsigma^2 respectively.
        end do
        INVMVNCOV(1:nvar,1:nvar)=MVNCOV(1:nvar,1:nvar)
        if (posdef(INVMVNCOV(1:nvar,1:nvar),nvar)==0) then
            write(*,*) 'Covariance matrix of LGRBs model not positive definite..cycling..'
            getLogPostProb=-huge(getLogPostProb)
            RETURN
        end if
        normfac=determinant(nvar,nvar,MVNCOV(1:nvar,1:nvar))
        if (normfac<=0) then
            write(*,*) 'Covariance determinant of SGRBs model is <=0',normfac,posdef(MVNCOV(1:nvar,1:nvar),nvar)
            do i=1,nvar
                write(*,*) (MVNCOV(i,j),j=1,nvar)
            end do
            STOP
        end if
        normfac=39.478417604357434*dsqrt(normfac)    ! (2*Pi)^(nvar/2)=39.478417604357434
        call inversematrix(nvar,nvar,MVNCOV(1:nvar,1:nvar),INVMVNCOV(1:nvar,1:nvar))
        !WRITE(*,*) 'HELLO!!!!'
    ! Do model integrations:
        !z=0.d0    ! This is only for the observer frame fitting
        
        CALL CPU_TIME(time_integration_start)
        
        !modelint=modelintz(z)
        
        call qromb(modelintz,zmin,zmax,modelint)
        
        !do i=1,nz
        !    call midexp(modelintz,zmin,zmax,modelint,i)
        !    if (modelint<0.d0) then
        !        !write(*,*) 'model_integral (variable modelint in getLogPostProb.f90) is <zero:',modelint
        !        write(*,*) 'Recycling...'
        !        !write(9834,'(1I30,16E30.5)') i,(PARAM(j),j=1,npar)  
        !        getLogPostProb=-huge(zmin)
        !        return
        !    end if
        !end do

        CALL CPU_TIME(time_integration_end)
        WRITE(888,'(2F25.9)') modelint,time_integration_end-time_integration_start

        !!! write(*,*) 'model integral: ',modelint
        if (modelint<=0.0d0) then
            write(*,*) 'model_integral (variable modelint in getLogPostProb.f90) is <=zero:',modelint
            write(*,*) 'Press Enter to continue...'
            read(*,*)
        end if
    ! End of model integration.
    ! Now calculate the likelihood of data:
        getLogPostProb = 0._RK
        GRBPROBINTEG: do idata=1,ndata
            ! Do redshift integration:
                !call qromo(probatz,zmin,zmax,probDetectionGRB,midexp)
                call qromb(probatz,zmin,zmax,probDetectionGRB)
                !z=0.d0
                !do i=1,nz
                !    call midexp(probatz,zmin,zmax,probDetectionGRB,i)
                !end do
                
                !probDetectionGRB=probatz(z)
                
                if (probDetectionGRB<=0.0d0) then
                    write(*,*) 'logprob is negative or zero:', probDetectionGRB
                    write(*,*) 'GRB number in the input file:', idata
                    !write(*,*) 'Press Enter to continue...'
                    !read(*,*)
                end if
                probDetectionGRB=probDetectionGRB/(modelint*normfac)
            ! End of redshift integration
            if (probDetectionGRB==0.d0) then
                getLogPostProb = -huge(probDetectionGRB)
                write(*,*) 'Log(Likelihood) is set to -huge'
                return
            else
                getLogPostProb = getLogPostProb + dlog(probDetectionGRB)
            end if
        end do GRBPROBINTEG
    END function getLogPostProb

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function probatz(z)
    use Batse_mod, only: ERFC_HEIGHT
    use MODELparameters
    use COSMOparameters
    use Zparameters    !, only: zplus1,logzplus1,lumdisterm,dvdzOzplus1
    use constants, only: sqrt2,sqrt2Pi,log10Mpc2cmSq4pi
    !? use GRBworld, only: npar
    use OBSGRBDATA
    use detection
    implicit none
    real(RK) :: probatz,erfcc,PbolEpk2P1024ph,ldisWickram,delayed_rate_Belz_Li    ! function
    real(RK) :: lumdisMPc,lumdisMPcSq
    real(RK) :: expterm,efficiency
    real(RK), intent(in) :: z
    zplus1      = z+1.0d0
    logzplus1   = dlog10(zplus1)
    lumdisMPc   = ldisWickram(z)
    lumdisMPcSq = lumdisMPc*lumdisMPc
    lumdisterm  = log10Mpc2cmSq4pi+dlog10(lumdisMPcSq)    ! This is Log10(4*pi*DL^2) where DL is luminosity distance in units of MPc
    ! Observed Data Probability
        X(1)=GRB%Event(idata)%logpbol+lumdisterm-loglisomean                ! logliso-loglisomean
        X(2)=GRB%Event(idata)%logepk+logzplus1-logepkzmean                ! logepkz-logepkzmean
        X(3)=GRB%Event(idata)%logsbol+lumdisterm-logzplus1-logeisomean    ! logeiso-logeisomean
#ifdef KFAC_ONETHIRD_ENABLED
        X(4)=GRB%Event(idata)%logt90-logzplus1*TIME_DILATION_EXPO-logt90zmean                ! logt90z-logt90zmean
#else
        X(4)=GRB%Event(idata)%logt90-logzplus1-logt90zmean                ! logt90z-logt90zmean
#endif
        expterm=0.5d0*dot_product(X,matmul(INVMVNCOV(1:nvar,1:nvar),X))
        !logp1024ph_eff=GRB%Event(idata)%logp1024ph+&
        !ERFC_HEIGHT*erfcc((GRB%Event(idata)%logt90-ERFC_AVG)/ERFC_STD)
        !if (logp1024ph_eff<eff_min_pph) then
        !    efficiency=0.d0
        if (GRB%Event(idata)%logPF53<eff_max_pph) then
            efficiency=0.5d0+0.5d0*(1.d0-erfcc((GRB%Event(idata)%logPF53-meanthresh)/sqrt2stdevthresh))
        else ! if (GRB%Event(idata)%logPF53>=eff_max_pph) then
            efficiency=1.d0
        !else
        !    write(*,*) 'Wrong efficiency limit in getLogPostProb.f90, eff_min_pph,eff_max_pph:',eff_min_pph,eff_max_pph
        !    STOP
        end if
        
        probatz=delayed_rate_Belz_Li(z)*efficiency/dexp(expterm)
        
        ! Calculate the LGRB rate at redshift z:
        ! ATTN: Note that a factor of CoverH*4.*pi is intentionally dropped in the calculation of dVdZ below for a higher efficiency of the code.
        ! Also the luminosity distance used here must be in units of MPc, if the Hubble parameter is also in units of MPc (which is, as in CoverH).
        ! Also here dvdz is infact divided by an extra zplus1 (which gives zplus1**3 in the denominator), again for efficiency purposes.
        !! not needed for delayed_rate: dvdzOzplus1=lumdisMPcSq/(zplus1**3*dsqrt(omega_DM*zplus1**3+omega_DE))    !*CoverH*4.d0*pi*MPc2cm*MPc2cm*MPc2cm
        ! note that dvdzOzplus1 is infact dvdz/(z+1) a factor that appears later in calculations, however is included here for code efficiency.
        ! note that dvdz lacks some constant factors, that I droped for higher code efficiency. dvdz is in units of MPc^3
        !! not needed for delayed rate:  if (z<z0) then
        !! not needed for delayed rate:      probatz=dvdzOzplus1*zplus1**g0*efficiency/dexp(expterm)            ! rate=dvdzOzplus1*zplus1**g0
        !! not needed for delayed rate:  else if (z>=z0 .and. z<z1) then
        !! not needed for delayed rate:      probatz=gamma1*dvdzOzplus1*zplus1**g1*efficiency/dexp(expterm)    ! rate=gamma1*dvdzOzplus1*zplus1**g1
        !! not needed for delayed rate:  else if (z>=z1) then
        !! not needed for delayed rate:      probatz=gamma2*dvdzOzplus1*zplus1**g2*efficiency/dexp(expterm)    ! rate=gamma2*dvdzOzplus1*zplus1**g2
        !! not needed for delayed rate:  else
        !! not needed for delayed rate:      write(*,*) 'invalid redshift input in Likelihood integration, z:', z
        !! not needed for delayed rate:      write(*,*) 'Press Enter to continue...'
        !! not needed for delayed rate:      read(*,*)
        !! not needed for delayed rate:  end if
        !probatz=efficiency/dexp(expterm)
        !write(*,*) 'expterm in parobatz.f90 = ',expterm
        !write(*,*) 'efficiency in parobatz.f90 = ',efficiency
        !write(*,*) 'Press Enter to continue...'
        !read(*,*)
    END function probatz

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function modelintz(z)    ! integral of GRB world model at given redshift z.
    use Batse_mod, only: ERFC_AVG
    use MODELparameters, only: sqrt2lisosigma,loglisomean
    use COSMOparameters
    use Zparameters    !, only: zplus1,lumdisterm,logzplus1,dvdz,probz
    use constants, only: log10Mpc2cmSq4pi
    !? use GRBworld, only: npar
    use OBSGRBDATA
    use detection
    implicit none
    real(RK), intent(in) :: z
    real(RK) :: lumdisMPc,lumdisMPcSq,erfcc
    real(RK) :: modelintz,ldisWickram,delayed_rate_Belz_Li        ! function
    !real(RK), EXTERNAL :: probatdurz    ! function
    zplus1      = z+1.0d0
    logzplus1   = dlog10(zplus1)
#ifdef KFAC_ONETHIRD_ENABLED
    logzplus1_dur = logzplus1*TIME_DILATION_EXPO
#else
    logzplus1_dur = logzplus1
#endif
    lumdisMPc   = ldisWickram(z)
    lumdisMPcSq = lumdisMPc*lumdisMPc
    lumdisterm  = log10Mpc2cmSq4pi+dlog10(lumdisMPcSq) ! This is Log10(4*pi*DL^2) where DL is luminosity distance in units of MPc
    ! INTEGRATION_LIMIT_LOGDUR_MIN & INTEGRATION_LIMIT_LOGDUR_MAX, is zero, such that at the given redshift z, that correspond to the rest-frame durations:
    meandurz = ERFC_AVG-logzplus1_dur
    logdurzmin = INTEGRATION_LIMIT_LOGDUR_MIN-logzplus1_dur
    logdurzmax = INTEGRATION_LIMIT_LOGDUR_MAX-logzplus1_dur
    ! This means 3 times above and below ERFC_AVG for which erf3==Erf(3.d0)=0.9999779095030014, which I assum to be 1.
    ! The two parameters below will be used for integration over Epk:
    logepkzmin=logepkmin+logzplus1
    logepkzmax=logepkmax+logzplus1
    ! Beginning of model integration.
    ! Here I divide the integration over duration into two parts: The duration distribution is a Gaussian with mean logdurzmean.
    ! Since the integration is over an infinite range of duration, the method used here is midexp algorithm.
    ! Beacuse midexp is only accurate for an ever_decreasing function, I use it only for the durz>logdurzmean+4sigma=uppert90zlim
    ! I use my own revised code midexp_mirror for the durz<logdurzmean-4sigma=lowert90zlim where the Gaussian is ever_increasing.
        ! Here is the integration between the two ranges:
            call qrombDur(probatdurz,logdurzmin,logdurzmax,modelintz)
        ! Now add the integration part for which the trigger efficiency is 100%, for any GRB observed duration:
            modelintz=modelintz+0.5d0*erfcc((glb_max_lpb+lumdisterm-loglisomean)/sqrt2lisosigma)    ! Integral of model at redshift z.
            
            modelintz=modelintz*delayed_rate_Belz_Li(z)
            
            !! not needed for delayed rate:  dvdzOzplus1=lumdisMPcSq/(zplus1**3*dsqrt(omega_DM*zplus1**3+omega_DE))    !*CoverH*4.d0*pi*MPc2cm*MPc2cm*MPc2cm
        ! note that dvdzOzplus1 is infact dvdz/(z+1) a factor that appears later in calculations, however is included here for code efficiency.
        ! note that dvdz lacks some constant factors, that I droped for higher code efficiency. dvdz is in units of MPc^3
            !! not needed for delayed rate:  if (z<z0) then
            !! not needed for delayed rate:      modelintz=modelintz*dvdzOzplus1*zplus1**g0            ! rate=dvdzOzplus1*zplus1**g0
            !! not needed for delayed rate:  else if (z>=z0 .and. z<z1) then
            !! not needed for delayed rate:      modelintz=modelintz*gamma1*dvdzOzplus1*zplus1**g1    ! rate=gamma1*dvdzOzplus1*zplus1**g1
            !! not needed for delayed rate:  else if (z>=z1) then
            !! not needed for delayed rate:      modelintz=modelintz*gamma2*dvdzOzplus1*zplus1**g2    ! rate=gamma2*dvdzOzplus1*zplus1**g2
            !! not needed for delayed rate:  else
            !! not needed for delayed rate:      write(*,*) 'invalid redshift input in Likelihood integration, z:', z
            !! not needed for delayed rate:      write(*,*) 'Press Enter to continue...'
            !! not needed for delayed rate:      read(*,*)
            !! not needed for delayed rate:  end if
    ! End of model integration.
    END function modelintz
!    ---------------------------------
    function probatdurz(logdurz)    ! This function refers to Integrates over Liso and Epk, by the use of the function ProbatLisoGivenDurz.
    use Batse_mod, only: ERFC_HEIGHT, ERFC_STD
    use MODELparameters    !, only: loglisomean,logt90zmean
    use COSMOparameters
    use Zparameters
    use detection
    implicit none
    real(RK) :: probatdurz,logdurz
    real(RK) :: erfcc    ! function
    real(RK), EXTERNAL :: ProbatLisoGivenDurz    ! function
    pph_correction=ERFC_HEIGHT*erfcc((logdurz-meandurz)/ERFC_STD)
    min_lpb_at_dur=eff_min_lpb+pph_correction
    conlisomeanD=aLD+bLD*logdurz
    conepkzmeanD=aED+bED*logdurz
    aELD=conepkzmeanD-conlisomeanD*bELD    ! term used in the coditional distribution of Epk given Luminosity.
    ! Integrate Epk-Liso bivariate distribution at the given logdurz:
        call qrombPbol(ProbatLisoGivenDurz,min_lpb_at_dur+lumdisterm,glb_max_lpb+lumdisterm,probatdurz)
        probatdurz=probatdurz/&
        (sqrt2Pit90zsigma*dexp(((logdurz-logt90zmean)/sqrt2t90zsigma)**2))
    END function probatdurz
!    ---------------------------------
    function ProbatLisoGivenDurz(logliso)
    use MODELparameters    !, only: sqrt2conlisosigmaD,sqrt2conepkzsigmaLD,sqrt2PiconlisosigmaD,aEL,bEL,conepkzmeanLD,loglisomean
    use Zparameters, only: lumdisterm
    use detection
    implicit none
    real(RK) :: logliso,ProbatLisoGivenDurz,efficiency    ! erfcc
    real(RK), EXTERNAL :: EpkzProbGivenLisoDurz
    conepkzmeanLD=aELD+bELD*logliso ! ATT: This is in fact "lumdisterm" minus the mean of the conditional dist. of Epkz given Liso.
    log10pbol=logliso-lumdisterm
    call qrombEpk(EpkzProbGivenLisoDurz,logepkzmin,logepkzmax,ProbatLisoGivenDurz)
    ! Now add the integration of the tails of the conditional Epk distribution given Liso (i.e. logliso) and logdurz.
        !XXX efficiency=0.5d0+0.5d0*(1.d0-erfcc((logp1024min+log10pbol-meanthresh_at_dur)/sqrt2stdevthresh))
        !XXX ProbatLisoGivenDurz=ProbatLisoGivenDurz+efficiency*&
        !XXX (0.5d0*erfcc((logepkzmax-conepkzmeanLD)/sqrt2conepkzsigmaLD)+&
        !XXX 1.0d0-0.5d0*erfcc((logepkzmin-conepkzmeanLD)/sqrt2conepkzsigmaLD))
    ProbatLisoGivenDurz=ProbatLisoGivenDurz/&
    (sqrt2PiconlisosigmaD*dexp(((logliso-conlisomeanD)/sqrt2conlisosigmaD)**2))
    END function ProbatLisoGivenDurz
!    ---------------------------------
    function EpkzProbGivenLisoDurz(logepkz)
    use MODELparameters !, only: sqrt2conepkzsigmaLD,sqrt2PiconepkzsigmaLD,conepkzmean
    use Zparameters, only: logzplus1
    use detection
    implicit none
    real(RK) :: EpkzProbGivenLisoDurz,erfcc,logepkz,efficiency,PbolEpk2P1024ph
    efficiency=0.5d0+0.5d0*&
    (1.d0-erfcc((PbolEpk2P1024ph(logepkz-logzplus1,log10pbol)-pph_correction-meanthresh)/sqrt2stdevthresh))
    EpkzProbGivenLisoDurz=efficiency/&
    (sqrt2PiconepkzsigmaLD*dexp(((logepkz-conepkzmeanLD)/sqrt2conepkzsigmaLD)**2))
    END function EpkzProbGivenLisoDurz

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module BatseSgrbWorldModel_mod