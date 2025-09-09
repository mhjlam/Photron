@echo off
REM Batch script to run Photron in release mode with proper DLL paths
REM Add vcpkg release bin directory to PATH
set PATH=%VCPKG_ROOT%\installed\x64-windows\bin;%PATH%

REM Run the executable with the first argument as config file
if "%~1"=="" (
    echo Usage: %0 ^<config-file^>
    echo Example: %0 config\config1.in
    exit /b 1
)

bin\Photron.exe %1
