from matplotlib import pyplot as pl
from itertools import islice

step = []
energy = []
with open('./energies.txt') as file:
    for line in islice(file, 0, 1):
        headers = line.strip().split(' ')

    for line in islice(file, 1, None):
        values = line.strip().split(' ')
        values = [float(i) for i in values]

        step.append(values[0])
        energy.append(values[2])

every = 1000
pl.plot(step[::every], energy[::every], 'b.')
pl.show()
