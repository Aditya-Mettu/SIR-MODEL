from json import dumps

fileNames = ["cityNames.txt", "stateNames.txt"]
toApply = ["cityNames", "stateNames"]
writeFileName = "cityAndStateNames.json"
toWrite = dict()

if __name__ == '__main__':
    for i in range(len(fileNames)):
        with open(fileNames[i]) as fptr:
            currentNames = [x.strip() for x in fptr.readlines()]

            toWrite[toApply[i]] = currentNames 

        with open(writeFileName, 'w') as f:
            f.write(dumps(toWrite))
