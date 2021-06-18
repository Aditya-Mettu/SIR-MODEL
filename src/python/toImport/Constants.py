from src.python.toImport.Clear import isWindows

dummyFileName = "data\dummy.csv" if isWindows() else "data/dummy.csv"
DateFormat = '%Y-%m-%d'