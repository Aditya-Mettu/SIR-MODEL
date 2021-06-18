from src.python.toImport import *

ifplotPR = True
dists = True
debug = False
ifprint = True 
span = 7

if __name__ == "__main__":
    clear()

    # To avoid scope related problems
    fullName = None
    sb_number = None 
    name = None 

    if dists:
        fullName = 'https://api.covid19india.org/csv/latest/districts.csv'
        name = getDistrictNames()

        tt1 = pd.to_datetime('2020-04-26', format = DateFormat)
        sb_number = 4
    else:
        fullName = 'https://api.covid19india.org/csv/latest/states.csv'
        name = getStateNames()

        tt1 = pd.to_datetime('2020-03-14', format = DateFormat)
        sb_number = 5

    urlwrite(fullName, dummyFileName)

    tableAll = pd.read_csv(dummyFileName, na_values={'Tested':['']})
    tableAll.Date = pd.to_datetime(tableAll.Date, format="%Y-%m-%d")

    span1 = 14

    mean_si, min_mean_si, max_mean_si = 4.7, 3.7, 6.0
    std_si, min_std_si, max_std_si = 2.9, 1.9, 4.9 

    for n in range(len(name)):

        Location = name[n]

        # To avoid scope related problems.
        tableLocation = None
        R, D, C, O, T = None, None, None, None, None

        dateInit = pd.to_datetime('2020-03-01', format = DateFormat)

        if dists:
            tableLocation = tableAll[(tableAll.District == Location) & (tableAll.Date >= dateInit)]

            R = tableLocation['Recovered'].to_numpy()
            D = tableLocation['Deceased'].to_numpy()
            C = tableLocation['Confirmed'].to_numpy()
            O = tableLocation['Other'].to_numpy()
            T = tableLocation['Tested'].to_numpy()
        else:
            tableLocation = tableAll[(tableAll.State == Location) & (tableAll.Date >= dateInit)]

            R = tableLocation['Recovered'].to_numpy()
            D = tableLocation['Deceased'].to_numpy()
            C = tableLocation['Confirmed'].to_numpy()
            O = tableLocation['Other'].to_numpy()
            T = tableLocation['Tested'].to_numpy()
 
        date = tableLocation['Date'].to_numpy()
        A = C - D - O - R # Active cases

        dailyC = np.diff(C)

        if debug:
            fig = plt.figure()

            plt.plot(date, C)
            plt.xlabel('Date')
            plt.ylabel('Confirmed Cases')
            
            plt.show()