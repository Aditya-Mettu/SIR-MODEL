import os
from sys import platform

def isWindows():
    return platform.lower().startswith("win")

def clear():    
    if isWindows(): # For Windows.
        os.system("cls")
    else: # For Linux and Mac.
        os.system("clear")
