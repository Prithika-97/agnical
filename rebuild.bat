@echo off
rem Re-run the full analysis after adding new FIRMS files to data\raw.
cd /d "%~dp0"
".venv\Scripts\python.exe" pipeline.py
pause
