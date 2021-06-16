# SIR-MODEL
This repository contains all files related to SIR model developed on Covid data

## How to run this project?

* Step 1: Move the SIR-MODEL directory.
* Step 2: Execute main.sh or main.bat depending on the OS.

## Notes:

* The main.bat or main.sh may need permission. Grant it execute permission if needed.
* python3 needs to be installed and environment variable must be set to python3 for Linux/Mac and python for Windows.

## Installing required modules for the code:

    pip install -r requirements.txt

## If in case you want to run without script:

    python -m src.python.SecondWave_All_SIR_MGDM

## To add cities and state names:

* Add names of the cities to data/cityNames.txt and add names of the states to data/stateNames.txt.
* Then run data/CityAndStateJsonMaker.py to update the JSON file.
* The program reads data from the cityAndStateNames.json file.
