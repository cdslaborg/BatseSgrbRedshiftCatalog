:: NOTE: This windows batch script builds the MergerDelayDist object files and the executable

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: build ParaMonte library and MergerDelayDist object files and executables
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set ERRORLEVEL=0
cd %~dp0

set BUILD_NAME=MergerDelayDist

echo.
echo. :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                       ::::
echo.                                                MergerDelayDist Build
echo. ::::                                                                                                                       ::::
echo. :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

echo. 
echo. Configuring build...
echo. 

REM set ParaMonte_ROOT_RELATIVE_PATH=..\..\..\20180101_ParaMonte\git

call configMergerDelayDist.bat
if %ERRORLEVEL%==1 (
    echo. 
    echo. -- !BUILD_NAME! - Fatal Error: Unable to configure and build flags. exiting...
    echo. 
    cd %~dp0
    exit /B 1
)
cd %~dp0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library type
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

echo !ParaMonte_LIB_NAME!|find "dynamic" >nul
if errorlevel 1 (
    set LTYPE=static
) else (
    echo.
    echo. -- !BUILD_NAME! - Fatal Error: ParaMonte library type could not be recognized.
    echo. -- !BUILD_NAME! - build failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

echo. -- !BUILD_NAME! - ParaMonte library name: !ParaMonte_LIB_NAME!
echo. -- !BUILD_NAME! - ParaMonte library type: !LTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library build
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set BTYPE=
echo !ParaMonte_LIB_NAME!|find "release" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "testing" >nul
    if errorlevel 1 (
        echo !ParaMonte_LIB_NAME!|find "debug" >nul
        if errorlevel 1 (
            echo.
            echo. -- !BUILD_NAME! - Fatal Error: ParaMonte library build could not be recognized.
            echo. -- !BUILD_NAME! - build failed. exiting...
            echo.
            cd %~dp0
            set ERRORLEVEL=1
            exit /B 1
        ) else (
            set BTYPE=debug
        )
    ) else (
        set BTYPE=testing
    )
) else (
    set BTYPE=release
)

echo. -- !BUILD_NAME! - ParaMonte library build: !BTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library parallelism
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set PTYPE=
set CAF_ENABLED=false
echo !ParaMonte_LIB_NAME!|find "mpi" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "cafshared" >nul
    if errorlevel 1 (
        echo !ParaMonte_LIB_NAME!|find "cafsingle" >nul
        if errorlevel 1 (
            set PTYPE=serial
            set CAFTYPE=none
        ) else (
            set CAFTYPE=single
            set PTYPE=cafsingle
            set CAF_ENABLED=true
        )
    ) else (
        set CAFTYPE=shared
        set PTYPE=cafshared
            set CAF_ENABLED=true
    )
) else (
    set PTYPE=mpi
)

echo. -- !BUILD_NAME! - ParaMonte library parallelism: !PTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: determine library target language
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set TARGET_LANG=
set CFI_ENABLED=false
echo !ParaMonte_LIB_NAME!|find "_fortran_" >nul
if errorlevel 1 (
    echo !ParaMonte_LIB_NAME!|find "_c_" >nul
    if errorlevel 1 (
        echo.
        echo. -- ParaMonte - Fatal Error: ParaMonte library target lanugage could not be recognized.
        echo. -- ParaMonte - build failed. exiting...
        echo.
        cd %~dp0
        set ERRORLEVEL=1
        exit /B 1
    ) else (
        set TARGET_LANG=C
        set CFI_ENABLED=true
    )
) else (
    set TARGET_LANG=Fortran
)

echo. -- !BUILD_NAME! - target language: !TARGET_LANG!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set default C/CPP/Fortran compilers/linkers
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set EXE_NAME=main.exe
set COMPILER_SUITE=intel
if !COMPILER_SUITE!==intel (

    set FCL=ifort
    set CCL=icl
    
    if !PTYPE!==mpi (
        set FCL=mpiifort -fc=ifort
    )

    REM set FCL_FLAGS=/threads /libs:static
    set FCL_FLAGS=/threads
    if !BTYPE!==debug   set FCL_FLAGS=!FCL_FLAGS! !INTEL_FORTRAN_DEBUG_FLAGS!
    if !BTYPE!==release set FCL_FLAGS=!FCL_FLAGS! !INTEL_FORTRAN_RELEASE_FLAGS!
    if !BTYPE!==testing set FCL_FLAGS=!FCL_FLAGS! !INTEL_FORTRAN_TESTING_FLAGS!
    REM if !LTYPE!==dynamic set FCL_FLAGS=!FCL_FLAGS! /fpp /DIS_COMPATIBLE_COMPILER

    set FL_FLAGS=

) else (

    echo. 
    echo. -- !BUILD_NAME! - Fatal Error: unsupported compiler suite: !COMPILER_SUITE!
    echo. -- !BUILD_NAME! - build failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B

)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set Fortran compiler version
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not defined COMPILER_VERSION (
    echo. -- !BUILD_NAME! - Detecting Fortran compiler version...
    cd %~dp0
    cd ..\lib\
    call getCompilerVersion.bat
    cd %~dp0
)

echo. -- !BUILD_NAME! - COMPILER_VERSION: !COMPILER_VERSION!
echo. -- !BUILD_NAME! - COMPILER_NAME: !FCL!
echo. -- !BUILD_NAME! - CFI_ENABLED: !CFI_ENABLED!
echo. -- !BUILD_NAME! - parallelism: !PTYPE!
echo. -- !BUILD_NAME! - build type: !BTYPE!
echo. -- !BUILD_NAME! - link type: !LTYPE!

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: initialize preprocessor flags
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set FPP_FLAGS=/fpp
if !CFI_ENABLED!==true set FPP_FLAGS=!FPP_FLAGS! /define:CFI_ENABLED
if !CAF_ENABLED!==true set FPP_FLAGS=!FPP_FLAGS! /define:CAF_ENABLED
if !PTYPE!==mpi set FPP_FLAGS=!FPP_FLAGS! /define:MPI_ENABLED

:: add Kfactor correction if needed

set FPP_FLAGS=!FPP_FLAGS! /define:!RATE_DENSITY_MODEL!

if !INTEGRATION_METHOD!==quadpackDPR set FPP_FLAGS=!FPP_FLAGS! /define:quadpackDPR
if !INTEGRATION_METHOD!==quadpackSPR set FPP_FLAGS=!FPP_FLAGS! /define:quadpackSPR
if !INTEGRATION_METHOD!==romberg set FPP_FLAGS=!FPP_FLAGS! /define:romberg

REM if !KFAC_CORRECTION!==kfacOneThird (
REM     set FPP_FLAGS=!FPP_FLAGS! /define:kfacOneThird
REM ) else (
REM     if not !KFAC_CORRECTION!==kfacNone (
REM         echo.
REM         echo. -- !BUILD_NAME! - Fatal Error occurred: KFAC_CORRECTION=!KFAC_CORRECTION! is not recognized as an option.
REM         echo. -- !BUILD_NAME! - exiting...
REM         echo.
REM         cd %~dp0
REM         set ERRORLEVEL=1
REM         exit /B 1
REM     )
REM )
set "KFAC_CORRECTION="
set MTYPE=!KFAC_CORRECTION!!RATE_DENSITY_MODEL!

:: set the executable's name

echo.
echo. -- !BUILD_NAME! - Fortran preprocessor macros: !FPP_FLAGS!
echo. -- !BUILD_NAME! - executable name: !EXE_NAME!
echo. -- !BUILD_NAME! - model: !MTYPE!
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set and make MergerDelayDist directories
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set ParaMonte_LIB_DIR=!ParaMonte_LIB_ROOT!\lib
set ParaMonte_MOD_DIR=!ParaMonte_LIB_ROOT!\mod

set MergerDelayDist_ROOT_PATH=%~dp0
set MergerDelayDist_BLD_DIR=!MergerDelayDist_ROOT_PATH!build\win!PLATFORM!\!COMPILER_SUITE!\!COMPILER_VERSION!\!BTYPE!\!LTYPE!\!PTYPE!\!MTYPE!\!INTEGRATION_METHOD!
set MergerDelayDist_SRC_DIR=!MergerDelayDist_ROOT_PATH!src
set MergerDelayDist_BIN_DIR=!MergerDelayDist_BLD_DIR!\bin
set MergerDelayDist_MOD_DIR=!MergerDelayDist_BLD_DIR!\mod
set MergerDelayDist_OBJ_DIR=!MergerDelayDist_BLD_DIR!\obj

:: loop over MergerDelayDist directories and generate them

echo.
for %%A in (
    !MergerDelayDist_BLD_DIR!
    !MergerDelayDist_BIN_DIR!
    !MergerDelayDist_MOD_DIR!
    !MergerDelayDist_OBJ_DIR!
    ) do (  if exist %%A (
                echo. -- %%A already exists. skipping...
            ) else (
                echo. -- !BUILD_NAME! - generating directory: %%A
                mkdir %%A
            )
)
echo.

if not !MergerDelayDist_OBJ_ENABLED!==true (
    echo.
    echo. -- !BUILD_NAME! - Warning: skipping object files build...
    echo.
    goto LABEL_MergerDelayDist_EXE_ENABLED
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: set and make MergerDelayDist directories
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: Read the name of each file from the ordered list of filenames in filelist.txt to compile

cd !MergerDelayDist_OBJ_DIR!
echo.
echo. -- !BUILD_NAME! - building for the rate model of !RATE_DENSITY_MODEL!...

:: First verify the source filelist exists

set MergerDelayDist_FILE_LIST=!MergerDelayDist_SRC_DIR!\filelist.txt
if not exist !MergerDelayDist_FILE_LIST! (
    echo.
    echo. -- !BUILD_NAME! - Fatal Error: The filelist.txt containing the source filenames does not exist. Path: !MergerDelayDist_FILE_LIST!
    echo. -- !BUILD_NAME! - build failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: generate object files

for /F "eol=! tokens=*" %%A in (!MergerDelayDist_FILE_LIST!) do (

    echo. 
    echo. -- !BUILD_NAME! - generating object file for %%A
    echo. 

    if !PTYPE!==mpi (
        call !FCL! !FCL_FLAGS! !FPP_FLAGS! ^
        /module:!MergerDelayDist_MOD_DIR!               %=path to output MergerDelayDist module files=% ^
        /I:!MergerDelayDist_MOD_DIR!                    %=path to output MergerDelayDist module files, needed 4 dependencies=%  ^
        /I:!ParaMonte_MOD_DIR!                     %=path to input ParaMonte module files=%  ^
        /c !MergerDelayDist_SRC_DIR!\%%A                %=path to input MergerDelayDist source file=%  ^
        || (
            echo. 
            echo. -- !BUILD_NAME! - Fatal Error: compilation of the object file for %%A failed.
            echo. -- !BUILD_NAME! - build failed. exiting...
            echo. 
            cd %~dp0
            set ERRORLEVEL=1
            exit /B
        )
    ) else (
        !FCL! !FCL_FLAGS! !FPP_FLAGS! ^
        /module:!MergerDelayDist_MOD_DIR!               %=path to output MergerDelayDist module files=% ^
        /I:!MergerDelayDist_MOD_DIR!                    %=path to output MergerDelayDist module files, needed 4 dependencies=%  ^
        /I:!ParaMonte_MOD_DIR!                     %=path to input ParaMonte module files=%  ^
        /c !MergerDelayDist_SRC_DIR!\%%A                %=path to input MergerDelayDist source file=%  ^
        || (
            echo. 
            echo. -- !BUILD_NAME! - Fatal Error: compilation of the object file for %%A failed.
            echo. -- !BUILD_NAME! - build failed. exiting...
            echo. 
            cd %~dp0
            set ERRORLEVEL=1
            exit /B
        )
    )
)
echo.

:LABEL_MergerDelayDist_EXE_ENABLED

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: generate MergerDelayDist executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if not !MergerDelayDist_EXE_ENABLED!==true (
    echo.
    echo. -- !BUILD_NAME! - Warning: skipping exectuable build...
    echo.
    goto LABEL_MergerDelayDist_RUN_ENABLED
)

echo.

if !LTYPE!==dynamic (

    echo.
    echo. -- !BUILD_NAME! - Fatal: dynamically-linked executable not implemented. This requires significant changes in the library interfaces.
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1

) else (    %= static linking requested =%

    echo. -- !BUILD_NAME! - generating statically-linked executable at: !MergerDelayDist_BIN_DIR!

    REM set EXE_NAME=!RATE_DENSITY_MODEL!.exe
    set REQUIRED_OBJECT_FILES=!MergerDelayDist_OBJ_DIR!\*.obj !ParaMonte_LIB_DIR!\!ParaMonte_LIB_NAME!

)

:: delete the old executable first
echo. -- !BUILD_NAME! - deleting old executable (if any) at: !MergerDelayDist_BIN_DIR!\!EXE_NAME!

cd !MergerDelayDist_BIN_DIR!
del !EXE_NAME!
if !ERRORLEVEL!==1 (
    echo. 
    echo. -- !BUILD_NAME! - Fatal Error: deletion of the old executable at !MergerDelayDist_BIN_DIR!\!EXE_NAME! failed. exiting...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)

:: build the executable

echo. 
echo. -- !BUILD_NAME! - Link command: !FCL! !FCL_FLAGS! !FL_FLAGS! !REQUIRED_OBJECT_FILES! /link /out:!MergerDelayDist_BIN_DIR!\!EXE_NAME!
echo. 

if !PTYPE!==mpi (
    cd !MergerDelayDist_BIN_DIR!
    call !FCL! !FCL_FLAGS! !FL_FLAGS! !REQUIRED_OBJECT_FILES! ^
    /module:!MergerDelayDist_MOD_DIR!               %=path to output MergerDelayDist module files=% ^
    /I:!MergerDelayDist_MOD_DIR!                    %=path to output MergerDelayDist module files, needed 4 dependencies=%  ^
    /I:!ParaMonte_MOD_DIR!                     %=path to input ParaMonte module files=%  ^
    /link /out:!EXE_NAME!
    REM call !FCL! !FCL_FLAGS! !FL_FLAGS! !REQUIRED_OBJECT_FILES! !MergerDelayDist_SRC_DIR!\ /link /out:!MergerDelayDist_BIN_DIR!\!EXE_NAME!
    cd %~dp0
) else (
    !FCL! !FCL_FLAGS! !FL_FLAGS! !REQUIRED_OBJECT_FILES! /link /out:!MergerDelayDist_BIN_DIR!\!EXE_NAME!
)

if !ERRORLEVEL!==1 ( 
    echo. 
    echo. -- !BUILD_NAME! - Fatal Error: linking of the object files might have failed.
    echo. -- !BUILD_NAME! - build might have failed. continuing...
    echo. 
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)
echo.

:: copy necessary input files in the executable's directory

echo. -- !BUILD_NAME! - copying input files to the executable's directory
echo. -- !BUILD_NAME! - from: !MergerDelayDist_ROOT_PATH!\in\!MTYPE!.nml   %= no need for final slash here =%
echo. -- !BUILD_NAME! -   to: !MergerDelayDist_BIN_DIR!\in\  %= final slash tells this is folder =%
xcopy /s /Y "!MergerDelayDist_ROOT_PATH!\in\!MTYPE!.nml" "!MergerDelayDist_BIN_DIR!\in\"
echo.
echo. -- !BUILD_NAME! - copying BATSE data to the executable's directory
echo. -- !BUILD_NAME! - from: !MergerDelayDist_ROOT_PATH!\in\!BATSE_DATA_FILE_NAME!   %= no need for final slash here =%
echo. -- !BUILD_NAME! -   to: !MergerDelayDist_BIN_DIR!\in\  %= final slash tells this is folder =%
xcopy /s /Y "!MergerDelayDist_ROOT_PATH!\in\!BATSE_DATA_FILE_NAME!" "!MergerDelayDist_BIN_DIR!\in\"
echo.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: run MergerDelayDist executable
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:LABEL_MergerDelayDist_RUN_ENABLED

:: run MergerDelayDist
:: if not !MergerDelayDist_RUN_ENABLED!==true goto LABEL_EXAMPLE_ENABLED
if not !MergerDelayDist_RUN_ENABLED!==true (
    echo.
    echo. -- !BUILD_NAME! - Warning: skipping the executable run...
    echo.
    goto :eof
)

echo. 
echo. :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo. ::::                                                                                                                       ::::
echo.                                             Running MergerDelayDist executable
echo. ::::                                                                                                                       ::::
echo. :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
echo.

editbin /STACK:99999999 !EXE_NAME!
REM if !LTYPE!==static (
REM     editbin /STACK:99999999 !EXE_NAME!
REM )

cd !MergerDelayDist_BIN_DIR!
!EXE_NAME! && ( 
    echo.
    echo.
    echo. -- !BUILD_NAME! - executable run successful. 
    echo.
) || ( 
    echo.
    echo.
    echo. -- !BUILD_NAME! - executable run failed. exiting...
    echo.
    cd %~dp0
    set ERRORLEVEL=1
    exit /B 1
)


cd %~dp0

set ERRORLEVEL=0
exit /B 0
