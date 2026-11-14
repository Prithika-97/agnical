@echo off
rem Start the AgniCal app. Double-click this file after setup.
cd /d "%~dp0"
echo Starting AgniCal... your browser will open at http://localhost:8501
echo Keep this window open while you use the app. Close it to stop the app.
start "" http://localhost:8501
".venv\Scripts\python.exe" -m streamlit run app.py --server.headless true --browser.gatherUsageStats false
pause
