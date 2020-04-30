program BatseWorldModelSimualtion

use, intrinsic :: iso_fortran_env, only: output_unit
use System_mod, only: CmdArg_type
use Batse_mod, only: readDataGRB
use WorldModelForBatseSGRB_mod, only: NPAR, getLogPostProb
use WorldModelForBatseSGRB_mod, only: zoneMin, zoneMax
use WorldModelForBatseSGRB_mod, only: zoneTol, lisoTol, epkzTol
use WorldModelForBatseSGRB_mod, only: zoneRef, lisoRef, epkzRef
use Constants_mod, only: IK, RK
use ParaDRAM_mod, only: runParaDRAM

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
    write(output_unit,"(*(g0))") "       a.exe <input file path: ../in/> <input file name: WorldModelSimualtionSGRB.nml>"
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
call readDataGRB(ngrb,inputBatseDataFile,outputBatseDataFile)

! start sampling the WoldModel's parameters
call runParaDRAM( ndim          = NPAR &
                , getLogFunc    = getLogPostProb &
                , inputFile     = CmdArg%Arg(1)%record//CmdArg%Arg(2)%record &
                )

end program BatseWorldModelSimualtion
