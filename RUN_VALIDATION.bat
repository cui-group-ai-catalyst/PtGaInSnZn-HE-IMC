@echo off
setlocal
cd /d "%~dp0"

where python >nul 2>nul
if errorlevel 1 (
  echo Python was not found on PATH.
  echo Create or activate the project environment, then run this file again.
  pause
  exit /b 1
)

python -c "import numpy, pandas, matplotlib" >nul 2>nul
if errorlevel 1 (
  echo The minimal validation dependencies are not available.
  echo Run: python -m pip install -r requirements_validation.txt
  pause
  exit /b 1
)

echo Running bounded Reviewer 2 Comment 1 validation...
python experimental_extensions\run_validation.py
if errorlevel 1 (
  echo Validation failed. Review the message above.
  pause
  exit /b 1
)

echo.
echo Validation passed. Opening the static evidence report...
start "" "experimental_extensions\outputs\system_validation\evidence_report.html"
pause
