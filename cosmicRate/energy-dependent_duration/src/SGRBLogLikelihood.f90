! Author: Amir Shahmoradi, last updated on Sunday 10:51 PM, Dec 22, 2013, IFS/ICMB, UT Austin

FUNCTION SGRBLogLikelihood(PARAM)
  USE MODELparameters
  USE Zparameters, ONLY: zmin,zmax
  USE constants, ONLY: sqrt2,sqrt2Pi
  USE GRBworld, ONLY: npar
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  INTEGER :: i,j,posdef
  DOUBLE PRECISION :: SGRBLogLikelihood,determinant
  DOUBLE PRECISION :: lisosigma,epkzsigma,normfac,t90zsigma,dumvar
  DOUBLE PRECISION, DIMENSION(npar) :: PARAM
  DOUBLE PRECISION, EXTERNAL :: probatz,modelintz
  
  DOUBLE PRECISION :: time_integration_end,time_integration_start
  
  loglisomean=PARAM(1)
  logepkzmean=PARAM(2)
  logeisomean=PARAM(3)
  logt90zmean=PARAM(4)
  lisosigma=10.d0**PARAM(5)
  epkzsigma=10.d0**PARAM(6)
  t90zsigma=10.d0**PARAM(8)
  do i=9,14
    dumvar=2.d0*PARAM(i)
    RHO(i)=(dexp(dumvar)-1.d0)/(dexp(dumvar)+1.d0)
  end do
  meanthresh=PARAM(15)
  stdevthresh=10.d0**PARAM(16)
  eff_min_pph=meanthresh-signif*stdevthresh
  eff_max_pph=meanthresh+signif*stdevthresh
  eff_min_lpb=eff_min_pph-logphminmaxdiff-logp1024min
  glb_max_lpb=eff_max_pph+lpb_correction
  rhoLE_given_Durz=(RHO(9)-RHO(11)*RHO(13))/dsqrt((1.d0-RHO(11)**2)*(1.d0-RHO(13)**2))
  conlisosigmaD=dsqrt(1.d0-RHO(11)**2)*lisosigma
  conepkzsigmaD=dsqrt(1.d0-RHO(13)**2)*epkzsigma
  conepkzsigmaLD=dsqrt(1.d0-rhoLE_given_Durz**2)*conepkzsigmaD
  bLD=RHO(11)*lisosigma/t90zsigma
  aLD=loglisomean-logt90zmean*bLD
  bED=RHO(13)*epkzsigma/t90zsigma
  aED=logepkzmean-logt90zmean*bED
  bELD=rhoLE_given_Durz*conepkzsigmaD/conlisosigmaD
  sqrt2lisosigma=sqrt2*lisosigma
  sqrt2t90zsigma=sqrt2*t90zsigma
  sqrt2Pit90zsigma=sqrt2Pi*t90zsigma
  sqrt2conlisosigmaD=sqrt2*conlisosigmaD
  sqrt2conepkzsigmaLD=sqrt2*conepkzsigmaLD
  sqrt2PiconlisosigmaD=sqrt2Pi*conlisosigmaD
  sqrt2PiconepkzsigmaLD=sqrt2Pi*conepkzsigmaLD
  sqrt2stdevthresh=sqrt2*stdevthresh
  do i=1,nvar
    MVNCOV(i,i)=10.d0**PARAM(i+4)
  end do
  MVNCOV(1,2)=MVNCOV(1,1)*MVNCOV(2,2)*RHO(9)
  MVNCOV(2,1)=MVNCOV(1,2)
  MVNCOV(1,3)=MVNCOV(1,1)*MVNCOV(3,3)*RHO(10)
  MVNCOV(3,1)=MVNCOV(1,3)
  MVNCOV(1,4)=MVNCOV(1,1)*MVNCOV(4,4)*RHO(11)
  MVNCOV(4,1)=MVNCOV(1,4)
  MVNCOV(2,3)=MVNCOV(2,2)*MVNCOV(3,3)*RHO(12)
  MVNCOV(3,2)=MVNCOV(2,3)
  MVNCOV(2,4)=MVNCOV(2,2)*MVNCOV(4,4)*RHO(13)
  MVNCOV(4,2)=MVNCOV(2,4)
  MVNCOV(3,4)=MVNCOV(3,3)*MVNCOV(4,4)*RHO(14)
  MVNCOV(4,3)=MVNCOV(3,4)
  do i=1,nvar
    MVNCOV(i,i)=MVNCOV(i,i)**2
  end do
  INVMVNCOV(1:nvar,1:nvar)=MVNCOV(1:nvar,1:nvar)
  if (posdef(INVMVNCOV(1:nvar,1:nvar),nvar)==0) then
    write(*,*) 'Covariance matrix of LGRBs model not positive definite..cycling..'
    SGRBLogLikelihood=-huge(SGRBLogLikelihood)
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
  normfac=39.478417604357434*dsqrt(normfac)
  call inversematrix(nvar,nvar,MVNCOV(1:nvar,1:nvar),INVMVNCOV(1:nvar,1:nvar))
  call CPU_TIME(time_integration_start)
  call qromb(modelintz,zmin,zmax,modelint)
  call CPU_TIME(time_integration_end)
  WRITE(888,'(2F25.9)') modelint,time_integration_end-time_integration_start
  if (modelint<=0.0d0) then
    write(*,*) 'model_integral (variable modelint in SGRBLogLikelihood.f90) is <=zero:',modelint
    write(*,*) 'Press Enter to continue...'
    read(*,*)
  end if
  GRBPROBINTEG: do idata=1,ndata
    call qromb(probatz,zmin,zmax,GRB(idata)%logprob)
    if (GRB(idata)%logprob<=0.0d0) then
      write(*,*) 'logprob is negative or zero:', GRB(idata)%logprob
      write(*,*) 'GRB number in the input file:', idata
    end if
    GRB(idata)%logprob=GRB(idata)%logprob/(modelint*normfac)
    if (GRB(idata)%logprob==0.d0) then
      GRB(idata)%logprob=-huge(zmin)
    else
      GRB(idata)%logprob=dlog(GRB(idata)%logprob)
    end if
  end do GRBPROBINTEG
  SGRBLogLikelihood=sum(GRB(1:ndata)%logprob)
  if (SGRBLogLikelihood<-huge(zmin)) then
    SGRBLogLikelihood=-huge(zmin)
    write(*,*) 'Log(Likelihood) is set to -huge'
  end if

END FUNCTION SGRBLogLikelihood

FUNCTION probatz(z)
  USE MODELparameters
  USE COSMOparameters
  USE Zparameters
  USE constants, ONLY: sqrt2,sqrt2Pi,log10Mpc2cmSq4pi
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: probatz,erfcc,PbolEpk2P1024ph,ldisWickram,delayed_rate_Belz_Li  ! function
  DOUBLE PRECISION :: lumdisMPc,lumdisMPcSq
  DOUBLE PRECISION :: expterm,efficiency
  DOUBLE PRECISION, INTENT(IN) :: z
  zplus1      = z+1.0d0
  logzplus1   = dlog10(zplus1)
  lumdisMPc   = ldisWickram(z)
  lumdisMPcSq = lumdisMPc*lumdisMPc
  lumdisterm  = log10Mpc2cmSq4pi+dlog10(lumdisMPcSq)
  X(1)=GRB(idata)%logpbol+lumdisterm-loglisomean    
  X(2)=GRB(idata)%logepk+logzplus1-logepkzmean      
  X(3)=GRB(idata)%logsbol+lumdisterm-logzplus1-logeisomean
  X(4)=GRB(idata)%logt90-logzplus1*0.66-logt90zmean       
  expterm=0.5d0*dot_product(X,matmul(INVMVNCOV(1:nvar,1:nvar),X))
  if (GRB(idata)%efflogp1024ph<eff_max_pph) then
    efficiency=0.5d0+0.5d0*(1.d0-erfcc((GRB(idata)%efflogp1024ph-meanthresh)/sqrt2stdevthresh))
  else 
    efficiency=1.d0
  end if
  probatz=delayed_rate_Belz_Li(z)*efficiency/dexp(expterm)
END FUNCTION probatz