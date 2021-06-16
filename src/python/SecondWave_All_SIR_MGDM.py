from src.python.toImport import *
from datetime import datetime

ifplotPR = True
dists = True
debug = False
ifprint = True 
span = 7

if __name__ == "__main__":
    clear()

    fullName = None
    sb_number = None 
    name = None 

    if dists:
        fullName = 'https://api.covid19india.org/csv/latest/districts.csv'
        name = getDistrictNames()

        tt1 = datetime(2020, 4, 26)
        sb_number = 4
    else:
        fullName = 'https://api.covid19india.org/csv/latest/states.csv'
        name = getStateNames()

        tt1 = datetime(2020, 3, 14)
        sb_number = 5

    urlwrite(fullName, dummyFileName)

    tableAll = pd.read_csv(dummyFileName)

    span1 = 14

    mean_si, min_mean_si, max_mean_si = 4.7, 3.7, 6.0
    std_si, min_std_si, max_std_si = 2.9, 1.9, 4.9 

    for n in range(len(name)):

        Location = name[n]

        if dists:
            pass
        else:
            pass


        print("Works Fine.")