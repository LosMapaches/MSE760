from itertools import islice
import pandas as pd
import os

def importdata(skip):
    directory = '../'
    currentdir = os.getcwd()
    folders = os.listdir(directory)
    folders = [i for i in folders if 'percent' in i]
    files = [currentdir+'/../'+i+'/data.txt' for i in folders]
    inputfile = [currentdir+'/../'+i+'/AlSm.in' for i in folders]

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
                           skiprows=skip
                           )

        runs[item.split('/')[-2]] = data

    param = {}
    for item in inputfile:
        with open(item) as file:
            for line in file:
                values = line.strip(' ').split(' ')
                if 'side' in values:
                    side = int(values[-1])

                if 'decimal' in values:
                    decimal = float(values[-1])

        atoms = 4*side**3
        param[item.split('/')[-2]] = {
                                      'atoms': atoms,
                                      'decimal': decimal
                                      }

    return runs, param
