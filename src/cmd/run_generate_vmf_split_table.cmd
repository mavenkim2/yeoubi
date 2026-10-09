@echo off
setlocal
pushd "%~dp0..\.."
if defined VCINSTALLDIR goto compiler_ready
set "VMF_VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VMF_VSWHERE%" goto failed
for /f "usebackq tokens=*" %%i in (`"%VMF_VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VMF_VS_PATH=%%i"
if not defined VMF_VS_PATH goto failed
call "%VMF_VS_PATH%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 goto failed
:compiler_ready
if not exist build mkdir build
cl /nologo /O2 /EHsc /std:c++17 /W4 /I src src\cmd\generate_vmf_split_table.cpp /Fo:build\generate_vmf_split_table.obj /Fe:build\generate_vmf_split_table.exe >nul
if errorlevel 1 goto failed
build\generate_vmf_split_table.exe %*
set "VMF_EXIT_CODE=%errorlevel%"
popd
exit /b %VMF_EXIT_CODE%
:failed
echo The vMF split-table build failed. Visual Studio C++ tools are required. 1>&2
popd
exit /b 1
