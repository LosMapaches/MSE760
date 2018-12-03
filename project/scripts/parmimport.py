from itertools import islice
import pandas as pd
import os


def importparm():

    directory = '..'
    currentdir = os.getcwd()
    systems = os.listdir(directory)
    systems = [i for i in systems if 'scripts' != i]
    systems = [i for i in systems if 'figures' != i]

    for run in systems:
        subruns = [i[0] for i in os.walk(directory+'/'+run)]
        subruns = subruns[1:]

    inputfile = [i+'/system.in' for i in subruns]

    param = {}
    for item in inputfile:
        with open(item) as file:
            for line in file:
                values = line.strip(' ').split(' ')
                if 'side' in values:
                    side = int(values[-1])

                if 'decimal' in values:
                    decimal = float(values[-1])

                if 'holdsteps' in values:
                    hold = int(values[-1])

                if 'mytimestep' in values:
                    timestep = float(values[-1])

                if 'mydumprate' in values:
                    rate = int(values[-1])

        atoms = 4*side**3

        item = item.split('.')[2]
        item = item.split('/')
        item = item[1]+'_'+item[2]

        param[item] = {
                       'atoms': atoms,
                       'decimal': decimal,
                       'hold': hold,
                       'timestep': timestep,
                       'rate': rate
                       }

    return param
