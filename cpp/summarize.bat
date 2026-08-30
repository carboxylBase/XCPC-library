@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

for /f %%i in ('powershell -command "Get-Date -Format \"yyyy-MM-dd HH:mm:ss\""') do set datetime=%%i

(
    echo # Summary
    echo Run summarize.bat to update this file.
    echo ## Last Updated: %datetime%
    echo ## Author: HuangZy
    echo [TOC]
    echo ## DataStructure
    for /r "DataStructure" %%f in (*.cpp) do (
        echo.
        echo ## %%~nxf
        echo ```cpp
        type "%%f"
        echo "\n"
        echo ```
    )

    for /r "Geometry" %%f in (*.cpp) do (
        echo.
        echo ## %%~nxf
        echo ```cpp
        type "%%f"
        echo "\n"
        echo ```
    )

    for /r "Graph" %%f in (*.cpp) do (
        echo.
        echo ## %%~nxf
        echo ```cpp
        type "%%f"
        echo "\n"
        echo ```
    )

    for /r "Math" %%f in (*.cpp) do (
        echo.
        echo ## %%~nxf
        echo ```cpp
        type "%%f"
        echo "\n"
        echo ```
    )
) > README.md

echo README.md has been written.
