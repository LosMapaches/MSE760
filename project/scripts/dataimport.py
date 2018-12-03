from itertools import islice
import pandas as pd
import os


def importdata(skip):

    directory = '..'
    currentdir = os.getcwd()
    systems = os.listdir(directory)
    systems = [i for i in systems if 'scripts' != i]
    systems = [i for i in systems if 'figures' != i]

    for run in systems:
        subruns = [i[0] for i in os.walk(directory+'/'+run)]
        subruns = subruns[1:]

    files = [i+'/data.txt' for i in subruns]

    runs = {}
    for item in files:
        with open(item) as file:
            for line in islice(file, 1, 2):
                headers = line.strip().split(' ')[1:]

        data = pd.read_csv(
                           item,
                           delimiter=' ',
                           comment='#',
                           names=headers,
                           skiprows=skip+2
                           )

        item = item.split('.')[2]
        item = item.split('/')
        item = item[1]+'_'+item[2]

        runs[item] = data

    return runs
