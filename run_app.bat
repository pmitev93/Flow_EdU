@echo off
cd /d "%~dp0"
Rscript -e "shiny::runApp('Flow_GUI.r', launch.browser=TRUE)"
