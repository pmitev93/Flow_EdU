#!/bin/bash
# Double-click this file to launch the Flow_EdU app
cd "$(dirname "$0")"
Rscript -e "shiny::runApp('Flow_GUI.r', launch.browser=TRUE)"
