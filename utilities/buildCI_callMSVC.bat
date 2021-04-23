rem ====================================================================================================   
rem                                                                                                        
rem   SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
rem            An implementation of the generalised motion modelling and image registration framework      
rem                                                                                                        
rem   Copyright (c) University College London (UCL). All rights reserved.                                  
rem                                                                                                        
rem   This software is distributed WITHOUT ANY WARRANTY; without even                                      
rem   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
rem   PURPOSE.                                                                                             
rem                                                                                                        
rem   See LICENSE.txt in the top level directory for details.                                              
rem                                                                                                        
rem ====================================================================================================   

@echo off

set PROJECT_NAME=%1
set TARGET_NAME=%2
set MSVC_VERSION=%3

rem The script is invoked from the gitlab-ci.yml job 
rem and will invoke the build of the specified project
rem using Microsoft Visual Studio 12 (2013).
rem Two parameters are expected
rem 1st --> the name of the project that will be opened by Visual Studio
rem 2nd --> the build target (such as ALL_BUILD, INSTALL, ...)
rem Currently only the 64-bit version is tested.

rem Enable command extensions to allow checking exit status of commands.
setlocal EnableExtensions

rem Allowed values of BTYPE are "x64" and "Win32"
if "%BTYPE%" == "" (
  set BTYPE=x64
)
if "%MSVC_VERSION%" == "12" (
	set "VS_DIR=c:/Program Files (x86)/Microsoft Visual Studio 12.0"

	if "%BTYPE%" == "x64" (
	  call "%VS_DIR%/VC/bin/amd64/vcvars64.bat"
	) 
	else (
	  call "%VS_DIR%/VC/bin/vcvars32.bat"
	)
)
if "%MSVC_VERSION%" == "15" (
	set "VS_DIR=D:/development/VisualStudio2017"

	if "%BTYPE%" == "x64" (
	  call "%VS_DIR%/VC/Auxiallry/Build/vcvars64.bat"
	) 
	else (
	  call "%VS_DIR%/VC/Auxiallry/Build/vcvars32.bat"
	)
)

rem set GIT_SSL_NO_VERIFY=1
set "VS_COMMAND=devenv.com"

echo Visual Studio folder:   %VS_DIR%
echo Visual Studio command:  %VS_COMMAND%
echo Bit depth:              %BTYPE%
echo Current work dir:       %cd%

rem stop visual studio recycling already running instances of msbuild.exe. we want clean ones.
rem http://stackoverflow.com/questions/12174877/visual-studio-2012-rtm-has-msbuild-exe-in-memory-after-close
set MSBUILDDISABLENODEREUSE=1

set BCONF=Release
set "VS_CONF=%BCONF%^|%BTYPE%"
echo Visual Studio config:   %VS_CONF%

rem The git usr/bin directory is needed for the 'tee' command.
set "PATH=%CMAKE_DIR%/bin;c:/Program Files/Git/bin;%VS_DIR%/Common7/IDE;%PATH%"

@echo on
%VS_COMMAND% /build %BCONF% /project %TARGET_NAME% /projectconfig %VS_CONF% %cd%/%PROJECT_NAME%.sln | tee %cd%/build.log 2>&1
rem %VS_COMMAND% /build %BCONF% /project %TARGET_NAME% /projectconfig %VS_CONF% %cd%/%PROJECT_NAME%.sln 
if %ERRORLEVEL% NEQ 0 exit /B 2
