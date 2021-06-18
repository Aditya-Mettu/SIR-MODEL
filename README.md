# SIR-MODEL
This repository contains all files related to SIR model developed on Covid data
We have made code in two formats.

1. Python-module
2. Standalone Jupyter notebook

## How to run this Python module?

* Step 1: Move the SIR-MODEL directory.
* Step 2: Execute main.sh or main.bat depending on the OS.

## Notes:

* The main.bat or main.sh may need permission. Grant it execute permission if needed.
* python3 needs to be installed and environment variable must be set to python3 for Linux/Mac and python for Windows.

## Installing required modules for the code:

    pip install -r requirements.txt

## If in case you want to run without script:

    git clone https://github.com/Aditya-Mettu/SIR-MODEL
    cd SIR-MODEL
    python -m src.python.SecondWave_All_SIR_MGDM

## To add cities and state names:

* Add names of the cities to data/cityNames.txt and add names of the states to data/stateNames.txt.
* Then run data/CityAndStateJsonMaker.py to update the JSON file.
* The program reads data from the cityAndStateNames.json file.

## Other details:

* The main python file is [SecondWave_All_SIR_MGDM.py]
* It uses a toImport module which has a functions module in it.
* toImport mainly defines every module and other functions needed for the main file to run.
* toImport has some constants.py with dummyFileName and clear.py which has clear function which is matlab clc equivalent.
* functions sub-module has getDistrictNames and getStateNames functions in GetDistrictAndStateNames.py and urlwrite in URLWrite.py
