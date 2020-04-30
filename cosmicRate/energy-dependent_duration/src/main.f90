! Amir Shahmoradi, Sunday 9:37 PM, Dec 10, 2013, IFS/ICMB, The University of Texas at Austin

include 'modules.f90'

program GRBWORLDMODEL

  use GRBworld
  use OBSGRBDATA
  use detection
  use Zparameters, ONLY: z0,z1,g0,g1,g2,gamma1,gamma2
  use constants, ONLY: log10Mpc2cmSq4pi
  use Bandmodel
  
  implicit none

  integer :: i,j
  integer, parameter :: iter=10**5,ireport=iter/10**4,idelay=10**5,idum_in=-4735
  integer, parameter :: mcmc_unit=22,mcmc_timing_unit=23
  double precision, parameter :: scalefactor=2.38**2/dble(npar)  ! Haar (2001) scale factor.
  double precision, external :: SGRBLogLikelihood  ! function
  double precision :: zprob,PbolEpk2P1024ph_default,BandFluxErgs,BandFluxPh  ! function
  double precision :: dummy,normfact,erfcc
  character(LEN=72) :: labels
  character(LEN=30), dimension(npar) :: LABEL
  character(LEN=300) :: file_iniparam,file_data_in,file_data_out,file_mcmc_out,file_mcmc_timing_out,file_integration_test
  character(LEN=300) :: path='..\'
  
  interface
    subroutine AMH(npar,iter,ireport,idelay,X,COV,SR,func,mcmc_unit,mcmc_timing_unit,LABEL,idum_in)
      inTEGER, intent(in) :: npar,iter,ireport,idelay,mcmc_unit,mcmc_timing_unit
      inTEGER, OPTIONAL :: idum_in
      double precision, EXTERNAL :: func
      double precision, dimension(npar), intent(inout) :: X
      double precision, dimension(npar,npar), intent(inout) :: COV
      double precision, dimension(npar,2), intent(in) :: SR
      CHARACTER(LEN=30), intent(in), dimension(npar) :: LABEL
    end subroutine AMH
  end interface
  
  gamma1=(1.+z0)**(g0-g1)
  gamma2=gamma1*(1.+z1)**(g1-g2)

  ! Read the command line arguments.
  if (command_argument_count()/=6) then
      if (command_argument_count()==0) then
        write(*,*)
    write(*,*) "No input arguments found on the command line."
        write(*,*)
        file_iniparam = trim(adjustl(path))//'in\iniparam_565_DRBL.in'
        file_data_in  = trim(adjustl(path))//'in\data_565.in'
        file_data_out = trim(adjustl(path))//'out\data_565.out'
        file_mcmc_out = trim(adjustl(path))//'out\AMHMCMC_DRBL.out'
        file_mcmc_timing_out = trim(adjustl(path))//'out\AMHMCMC_timing_DRBL.out'
        file_integration_test = trim(adjustl(path))//'out\integration_test_timing.out'
      else
    write(*,*)
    write(*,*) "Invalid number of input arguments on the command line."
    write(*,*) "Correct use:"
    write(*,'(1A120)') "./a.out <iniparam.txt> <data.in> <data.out> <AMHMCMC.out> <AMHMCMC_TIMinG.out> <file_integration_test>"
    write(*,*)
    stop
    end if
    else
    call get_command_argument(1,file_iniparam)
    call get_command_argument(2,file_data_in)
    call get_command_argument(3,file_data_out)
    call get_command_argument(4,file_mcmc_out)
    call get_command_argument(5,file_mcmc_timing_out)
    call get_command_argument(6,file_integration_test)
  end if
  
  ! Read input GRB data.
  open(unit=11, file=file_data_in, status = 'old')
  open(unit=21, file=file_data_out, status = 'replace')
  open(unit=mcmc_unit,file=file_mcmc_out,status='replace')
  open(unit=mcmc_timing_unit,file=file_mcmc_timing_out,status='replace')
  open(unit=888,file=file_integration_test,status='replace')
  write(888,'(2A25)') 'Integral','time_taken'
  
  write(21,*) 'comment: BATSE_565_SGRB_bol64(0.001,20000)_data_6904'
  write(21,'(23A30)') 'trigger',&
            'Pbol64[0.001,20000]keV','Sbol[0.001,20000]keV','Epk','T90',&
            'Log(Pbol64[0.001,20000]keV)','Log(Sbol[0.001,20000]keV)','Log(Epk)','Log(T90)','Log(P64ph)','Log(Eff_Trig_PF)',&
            'Log(EPR)','Log(EFR)','Log(ETR)','Log(FPR)','Log(TPR)','Log(TFR)',&
            'Log(EPM)','Log(EFM)','Log(ETM)','Log(FPM)','Log(TPM)','Log(TFM)'
  
  read(11,*) labels; read(11,*) labels
  do i=1,ndata
    read(11,*) GRB(i)%trigger,GRB(i)%logp1024ph,GRB(i)%logsbol,GRB(i)%logepk,GRB(i)%logt90
    GRB(i)%efflogp1024ph=GRB(i)%logp1024ph-Serfc*erfcc((GRB(i)%logt90-meandur)/scaledur)
    GRB(i)%logpbol=GRB(i)%logp1024ph-PbolEpk2P1024ph_default(GRB(i)%logepk)
    normfact=GRB(i)%logsbol-dlog10(BandFluxErgs(eminbol,emaxbol,alpha,beta,10.d0**GRB(i)%logepk))
    GRB(i)%logsbol=normfact+dlog10(BandFluxPh(emindet,emaxdet,alpha,beta,10.d0**GRB(i)%logepk))
    GRB(i)%logsbol=GRB(i)%logsbol-PbolEpk2P1024ph_default(GRB(i)%logepk)
    
    write(21,'(I30,22E30.6)') GRB(i)%trigger,&
    10.d0**GRB(i)%logpbol,10.d0**GRB(i)%logsbol,10.d0**GRB(i)%logepk,10.d0**GRB(i)%logt90,&
    GRB(i)%logpbol,GRB(i)%logsbol,GRB(i)%logepk,GRB(i)%logt90,GRB(i)%logp1024ph,GRB(i)%efflogp1024ph,&
    GRB(i)%logepk-GRB(i)%logpbol,GRB(i)%logepk-GRB(i)%logsbol,GRB(i)%logepk-GRB(i)%logt90,&
    GRB(i)%logsbol-GRB(i)%logpbol,GRB(i)%logt90-GRB(i)%logpbol,GRB(i)%logt90-GRB(i)%logsbol,&
    GRB(i)%logepk+GRB(i)%logpbol,GRB(i)%logepk+GRB(i)%logsbol,GRB(i)%logepk+GRB(i)%logt90,&
    GRB(i)%logsbol+GRB(i)%logpbol,GRB(i)%logt90+GRB(i)%logpbol,GRB(i)%logt90+GRB(i)%logsbol
  end do
  
  close(11); close(21)
  
  open(unit=13,file=file_iniparam,status='old')
  read(13,*) labels
  read(13,*) labels,(STCV(j),j=1,npar)
  write(*,*) 'MEAN of the parameters:'
  write(*,*) labels,(STCV(j),j=1,npar)
  read(13,*) labels,(STDEV(j),j=1,npar)
  write(*,*) 'Standard Deviations of the parameters:'
  write(*,*) labels,(STDEV(j),j=1,npar)
  read(13,*) labels
  write(*,*) 'Correlation Matrix:'
  
  do i=1,npar
    read(13,*) LABEL(i),(COR(i,j),j=1,npar)
    write(*,*) LABEL(i),(COR(i,j),j=1,npar)
    !read(*,*)
  end do
  
  write(*,*) 'Covariance Matrix:'
  
  do i=1,npar
    do j=1,npar
      COV(i,j)=STDEV(i)*STDEV(j)*COR(i,j)*scalefactor
    end do
    write(*,*) (COV(i,j),j=1,npar)
  end do
  close(13)
  
  SR(1,1)  = -huge(scalefactor)  
  SR(2,1)  = -huge(scalefactor)  
  SR(3,1)  = -huge(scalefactor)  
  SR(4,1)  = -huge(scalefactor)  
  SR(5,1)  = -huge(scalefactor)  
  SR(6,1)  = -huge(scalefactor)  
  SR(7,1)  = -huge(scalefactor)  
  SR(8,1)  = -huge(scalefactor)  
  SR(9,1)  = -huge(scalefactor)  
  SR(10,1) = -huge(scalefactor)  
  SR(11,1) = -huge(scalefactor)  
  SR(12,1) = -huge(scalefactor)  
  SR(13,1) = -huge(scalefactor)  
  SR(14,1) = -huge(scalefactor)  
  SR(15,1) = -huge(scalefactor)  
  SR(16,1) = -huge(scalefactor)  
  SR(1,2)  =  huge(scalefactor)  
  SR(2,2)  =  huge(scalefactor)  
  SR(3,2)  =  huge(scalefactor)  
  SR(4,2)  =  huge(scalefactor)  
  SR(5,2)  =  huge(scalefactor)  
  SR(6,2)  =  huge(scalefactor)  
  SR(7,2)  =  huge(scalefactor)  
  SR(8,2)  =  huge(scalefactor)  
  SR(9,2)  =  huge(scalefactor)  
  SR(10,2) =  huge(scalefactor)  
  SR(11,2) =  huge(scalefactor)  
  SR(12,2) =  huge(scalefactor)  
  SR(13,2) =  huge(scalefactor)  
  SR(14,2) =  huge(scalefactor)  
  SR(15,2) =  huge(scalefactor)  
  SR(16,2) =  huge(scalefactor)  

  ! Start the MCMC sampling
  call AMH(npar,iter,ireport,idelay,STCV,COV,SR,SGRBLogLikelihood,mcmc_unit,mcmc_timing_unit,LABEL,idum_in)

end program GRBWORLDMODEL

include 'SGRBLogLikelihood.f90'
include 'modelint.f90'
include 'AMH_TOF.f90'
include 'choldc.f90'
include 'gasdevran2.f90'
include 'MVNRND.f90'
include 'ran2.f90'
include 'samcovmat.f90'
include 'posdef.f90'
include 'erfcc.f90'
include 'inversematrix.f90'
include 'lubksb.f90'
include 'ludcmp.f90'
include 'determinant.f90'
include 'qrombPbol.f90'
include 'qrombEpk.f90'
include 'PbolEpk2P1024ph_default.f90'
include 'PbolEpk2P1024ph.f90'
include 'qromo.f90'
include 'qromb.f90'
include 'trapzd.f90'
include 'midexp.f90'
include 'polint.f90'
include 'ldisWickram.f90'
include 'BandFluxErgs.f90'
include 'BandFluxPh.f90'
include 'qrombDur.f90'
include 'delayed_rate_Belz_Li.f90'