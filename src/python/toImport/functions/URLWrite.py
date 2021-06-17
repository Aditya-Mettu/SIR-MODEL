import requests
import csv

def urlwrite(url, fileName):
    response = requests.get(url)

    for i in range(5):
        if response.status_code == 200:
            break
        response = requests.get(url)
    else:
        raise Exception("Internet not working")
    
    with open(fileName, 'w') as fptr:
        fptr.write(response.text)
