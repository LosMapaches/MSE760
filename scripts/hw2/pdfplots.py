from matplotlib import pyplot as pl
from itertools import islice
import os

path = './pdf/'
files = os.listdir(path)

count = 0
for item in files:
    files[count] = path+files[count]
    count += 1

data = {}
for item in files:
    data[item] = {'gr': [], 'dist': []}

    with open(item) as file:
        for line in islice(file, 1, None):
            values = line.strip().split(' ')
            values = [float(i) for i in values]
            data[item]['gr'].append(values[0])
            data[item]['dist'].append(values[1])

for key in data:
    pl.plot(data[key]['dist'], data[key]['gr'])
    pl.xlabel('Distance [m]')
    pl.ylabel('g(r)')
    pl.tight_layout()
    pl.grid()
    name = key.split('/')[-1]
    pl.savefig('./figures/'+name)
    pl.clf()
