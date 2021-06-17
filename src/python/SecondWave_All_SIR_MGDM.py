from src.python.toImport import *

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

    tableAll = pd.read_csv(dummyFileName, na_values={'Tested':['']})
    tableAll.Date = pd.to_datetime(tableAll.Date, format="%Y-%m-%d")

    span1 = 14

    mean_si, min_mean_si, max_mean_si = 4.7, 3.7, 6.0
    std_si, min_std_si, max_std_si = 2.9, 1.9, 4.9 

    for n in range(len(name)):

        Location = name[n]
        tableLocation = None

        if dists:
            tableLocation = tableAll[(tableAll.District == Location) & (tableAll.Date >= pd.to_datetime('1-3-2020', format='%d-%m-%Y'))]

            R = tableLocation['Recovered']
            D = tableLocation['Deceased']
            C = tableLocation['Confirmed']
            O = tableLocation['Other']
            T = tableLocation['Tested']

            print(tableLocation)

        else:
            tableLocation = tableAll[(tableAll.State == Location) & (tableAll.Date >= pd.to_datetime('1-3-2020', format='%d-%m-%Y'))]

            R = tableLocation['Recovered']
            D = tableLocation['Deceased']
            C = tableLocation['Confirmed']
            O = tableLocation['Other']
            T = tableLocation['Tested']

        print(tableLocation)