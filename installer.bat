@echo off
setlocal enabledelayedexpansion

echo         Installing Apsis Interpreter

set "INSTALL_DIR=%LOCALAPPDATA%\Apsis"
set "EXE_PATH=%INSTALL_DIR%\apsis.exe"

if not exist "%INSTALL_DIR%" (
    mkdir "%INSTALL_DIR%"
)

if exist "build\Release\apsis.exe" (
    copy /Y "build\Release\apsis.exe" "%EXE_PATH%"
) else if exist "apsis.exe" (
    copy /Y "apsis.exe" "%EXE_PATH%"
) else (
    echo Building Apsis using CMake...
    cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release
    copy /Y "build\Release\apsis.exe" "%EXE_PATH%"
)

echo Checking System PATH...
powershell -NoProfile -Command ^
    "$oldPath = [Environment]::GetEnvironmentVariable('Path', 'User');" ^
    "if ($oldPath -notlike '*%INSTALL_DIR%*') {" ^
    "    [Environment]::SetEnvironmentVariable('Path', $oldPath + ';%INSTALL_DIR%', 'User');" ^
    "    Write-Host 'Added Apsis to User PATH successfully.';" ^
    "} else {" ^
    "    Write-Host 'Apsis is already in your PATH.';" ^
    "}"

echo.
echo  Installation Complete!
echo  Please restart your terminal/PowerShell.
echo  Usage: apsis script.jdx
pause