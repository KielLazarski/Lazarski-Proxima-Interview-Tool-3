@echo off
rem GlueForge. Double-click this; a page opens in your browser and everything
rem else happens there, including installing what is missing.
rem
rem Leave this window open. It *is* the server -- the page in your browser
rem talks to it, and closing this window is what makes the page start saying
rem "the GlueForge server is not answering".
setlocal
cd /d "%~dp0"

if not exist "glueforge.py" goto :inzip

rem Finding a Python that actually runs.
rem
rem `where python` is not enough, and this is the trap: Windows ships a stub
rem at WindowsApps\python.exe that exists, answers `where`, and then refuses
rem to run, printing "Python was not found; run without arguments to install
rem from the Microsoft Store". So each candidate is tested by running it.
set "PY="
for %%C in (py python python3) do (
  if not defined PY (
    %%C -c "import sys; sys.exit(0)" >nul 2>nul && set "PY=%%C"
  )
)
if not defined PY goto :nopython

echo.
echo   Starting GlueForge with %PY%.
echo   Leave this window open while you use it.
echo.
%PY% glueforge.py %*
set "CODE=%ERRORLEVEL%"

rem A server that stops on its own has printed why, and without this the
rem window closes over the explanation and the whole thing just looks like
rem nothing happened.
if not "%CODE%"=="0" (
  echo.
  echo   GlueForge stopped with error %CODE%. The reason should be above.
  echo.
  pause
)
goto :end

:nopython
echo.
echo   No working Python was found.
echo.
echo   If Windows just offered to install Python from the Microsoft Store,
echo   that is a stub rather than Python itself. Close it and install from
echo   https://www.python.org/downloads/ instead, ticking
echo   "Add python.exe to PATH" on the first screen of the installer.
echo.
echo   Everything else GlueForge needs installs from inside the app.
echo.
pause
goto :end

:inzip
echo.
echo   glueforge.py is not next to this file.
echo.
echo   This usually means you are running it from inside the zip. Windows
echo   lets you open a zip and double-click what is in it, but the files
echo   land in a temporary folder that gets wiped, so nothing works.
echo.
echo   Right-click the zip, choose "Extract All...", pick a real folder,
echo   and run GlueForge.bat from there instead.
echo.
pause

:end
endlocal
