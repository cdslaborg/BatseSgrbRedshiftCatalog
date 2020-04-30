! This version of the program now takes in the name of the input and output files as command line arguments.
! INPUT:    1.  Data file, typically named: BATSE_565_SGRBs_P64ph_Sbol_Epk_T90_6904.in.
!                The relevant parameters and structure of data are given in the module GRB OBSGRBDATA in modules.f90
!            2.    The MCMC initialization parameters are read in from the file: named iniparam.in
! OUTPUT: 
!            1.    Revised input data is dumped into the file BATSE_565_SGRB_bol64(0.001,20000)_data_6904
!            2.    The MCMC chain of parameters is dumped into AMHMCMC_SAMPLE.txt
!            3.    Relevant information about the simulation progress and the MCMC chain is dumped into the file MCMC_Time.txt.

! Amir Shahmoradi, Sunday 9:37 PM, Dec 10, 2013, IFS/ICMB, The University of Texas at Austin

!include 'modules.f90'
!include 'BatseSgrbWorldModel_mod.f90'


program BatseWorldModelSimualtion

use, intrinsic :: iso_fortran_env, only: output_unit
use BatseSgrbWorldModel_mod, only: getLogPostProb, NPAR
use Constants_mod, only: IK, RK
use System_mod, only: CmdArg_type
use ParaDRAM_mod, only: runParaDRAM

use GRBworld
use OBSGRBDATA
use detection
use Zparameters, ONLY: z0,z1,g0,g1,g2,gamma1,gamma2
use constants, ONLY: log10Mpc2cmSq4pi
use Bandmodel

implicit none

integer(IK)                 :: i, j, inFileUnit, ngrb
character(:), allocatable   :: inputBatseDataFile, outputBatseDataFile
type(CmdArg_type)           :: CmdArg

namelist /InputData/ ngrb, inputBatseDataFile, outputBatseDataFile
namelist /InputData/ zoneMin, zoneMax
namelist /InputData/ zoneTol, lisoTol, epkzTol
namelist /InputData/ zoneRef, lisoRef, epkzRef

! query input data file name from the command line
call CmdArg%query()
if (CmdArg%count/=2) then
    write(output_unit,"(*(g0))")
    write(output_unit,"(*(g0))") "FATAL: Invalid number of command-line arguments: ", CmdArg%count
    write(output_unit,"(*(g0))") "       Use the following example syntax to invoke the program: "
    write(output_unit,"(*(g0))") "       a.exe <input file path ending with a slash separator: ../in/> <input file name: WorldModelSimualtionSGRB.nml>"
end if

! read simulation input data
open( newunit = inFileUnit, file = CmdArg%Arg(1)%record//CmdArg%Arg(2)%record, status="old" )
    allocate(character(1000) :: inputBatseDataFile, outputBatseDataFile)
    read(inFileUnit,nml=InputData)
    inputBatseDataFile  = CmdArg%Arg(1)%record//trim(adjustl(inputBatseDataFile))
    outputBatseDataFile = CmdArg%Arg(1)%record//trim(adjustl(outputBatseDataFile))
    write(*,*) "zoneTol: ", zoneTol
    write(*,*) "lisoTol: ", lisoTol
    write(*,*) "epkzTol: ", epkzTol
close(inFileUnit)

!***********************************************************************************************************************************
!***********************************************************************************************************************************

! read observed grb input data
call readDataGRB( inputBatseDataFile    &
                , outputBatseDataFile   &
                , isLGRB = .false.      &
                )

! start sampling the WoldModel's parameters
call runParaDRAM( ndim          = NPAR &
                , getLogFunc    = getLogPostProb &
                , inputFile     = CmdArg%Arg(1)%record//CmdArg%Arg(2)%record &
                )

end program BatseWorldModelSimualtion

!include 'SGRBLogLikelihood.f90'
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
include 'qrombDur.f90'
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
include 'delayed_rate_Belz_Li.f90'