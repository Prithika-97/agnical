@echo off
rem Open the AgniCal website locally. Double-click this file after setup.
cd /d "%~dp0"
echo Opening AgniCal at http://localhost:8000
echo Keep this window open while you use the site. Close it to stop.
start "" http://localhost:8000
".venv\Scripts\python.exe" -m http.server 8000 -d docs
pause
