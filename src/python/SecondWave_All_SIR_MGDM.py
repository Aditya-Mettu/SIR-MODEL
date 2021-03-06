from src.python.toImport import *
import warnings

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

    tableAll = pd.read_csv(dummyFileName, na_values={'Tested':['', 0]})
    tableAll.Date = pd.to_datetime(tableAll.Date, format="%Y-%m-%d")
    tableAll.fillna(method = 'bfill', inplace=True)

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
            fig.add_subplot(2, 1, (1))

            plt.stem(date, C)
            plt.xlabel('Date')
            plt.ylabel('Confirmed Cases')
            
            fig.add_subplot(2, 1, (2))

            plt.stem(date[1:], dailyC)
            plt.xlabel('Date')
            plt.ylabel('Infection rate')

            plt.show()

        PR, CFR, Testing = None, None, None
        if not dists:    
            # PR is not used for districts, it is just for the states.
            PR = (np.diff(movmean(C,span)) / np.diff(movmean(T,span))) *100
            Testing = np.diff(T)
        CFR = (np.diff(movmean(D,span)) / np.diff(movmean(C,span))) *100

        date0 = np.copy(date)
        tt1 = None

        
        if dists:
            tt1 = str(date0[0])
            sb_number = 4
        else:
            tt1 = str(date0[0])
            sb_number = 5
        
        date1 = pd.to_datetime('2021-4-1', format = DateFormat)

        dateEnd = date[-1]

        tableLocation1 = tableLocation[(tableLocation.Date > date1) & (tableLocation.Date < dateEnd)]

        date = tableLocation1['Date'].to_numpy()
        C = tableLocation1['Confirmed'].to_numpy()
        D = tableLocation1['Deceased'].to_numpy()
        R = tableLocation1['Recovered'].to_numpy()
        Npop = int('80e6', 16)

        init = Init(Location, date, C, D, R, Npop)
        # Init is a simple class to save data like a struct object.

        diffC = np.diff(init.C)

        q1 = diffC[0]
        q2 = diffC[-1]
        t = diffC.shape[0]

        ## Using exponential model

        Cexp = np.copy(init.C)
        # print(Cexp)
        texp = np.array(list(range(1, Cexp.shape[0] + 1)))

        f = fit(texp, Cexp, 'exp1')

        coeff = list(f)
        td_2 = log(2)/coeff[1]

        tnew = texp + 120

        EXP = ModelClass(date, funcExp1(tnew, *f))

        if debug == 0:
            fig = plt.figure()
            plt.plot(texp, EXP.C)
            plt.show()

        ## Using Logistic model

        flog = fit(texp, Cexp, 'logistic')

        ## Using For coefficients
        
        logmdl = fit(texp, Cexp, 'logmdl')
        LOG = ModelClass(date, funcLogmdl(tnew, *logmdl))

        print("OK")
