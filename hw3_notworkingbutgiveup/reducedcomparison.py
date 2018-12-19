'''
Tred = 2.0
These are the values for the reduced, uncut pressures
'''

from matplotlib import pyplot as pl
from itertools import islice
import numpy as np
import os

files = os.listdir('.')
pressurefiles = [i for i in files if 'pressures' in i]

rhostate = np.arange(0.1, 1.0, 0.1)

pstate = [
          0.1776,
          0.329,
          0.489,
          0.700,
          1.071,
          1.75,
          3.028,
          5.285,
          9.12
          ]

fig, ax = pl.subplots()

ax.plot(
        rhostate,
        pstate,
        marker='.',
        color='b',
        linestyle='none',
        label='Equation of State'
        )

pvir = []
pideal = []
rhos = []
for item in pressurefiles:
    name = item.split('.')[0]
    print(name)

    rho = name.split('rho')[0]
    rho = rho.split('p')[-1]
    rho = '0.'+rho
    rho = float(rho)
    rhos.append(rho)

    pressures = []
    with open(item) as file:
        for line in islice(file, 1, None):
            pressures.append(float(line))

    pvir.append(pressures[0])
    pideal.append(pressures[1])

ax.plot(
        rhos,
        pvir,
        marker='.',
        color='r',
        linestyle='none',
        label='Virial'
        )

ax.plot(
        rhos,
        pideal,
        marker='.',
        color='g',
        linestyle='none',
        label='Ideal'
        )

ax.set_xlabel('LJ-Reduced Density')
ax.set_ylabel('LJ-Reduced Pressure')
ax.grid()
ax.legend(loc='best')
fig.tight_layout()
fig.savefig('reducedrho')
