from matplotlib import pyplot as pl
from itertools import islice

step = []
energy = []
with open('./0p84rho_energies.txt') as file:
    for line in islice(file, 0, 1):
        headers = line.strip().split(' ')

    for line in islice(file, 1, None):
        values = line.strip().split(' ')
        values = [float(i) for i in values]

        step.append(values[0])
        energy.append(values[2])

every = 1000
index1 = 400000
index2 = 1000000
average = sum(energy[index1:index2])/len(energy[index1:index2])

pl.axvline(x=index1, color='r', label='Start of Average')
pl.axvline(x=index2, color='r', label='End of Average')
pl.axhline(y=average, color='y', label='Settled Average: '+str(average)[:8]+' [eV/atom]')

pl.plot(step[::every], energy[::every], 'b.', label='Data Points')
pl.ylabel('Potential Energy [eV/atom]')
pl.legend(loc='best')
pl.xlabel('Step [-]')
pl.tight_layout()
pl.grid()
pl.savefig('part1')
pl.clf()
