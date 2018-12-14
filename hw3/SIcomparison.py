'''
Tred = 2.0
These are the values for the reduced, uncut pressures
'''

from matplotlib import pyplot as pl
from itertools import islice
import numpy as np
import os

epsilon = 0.0104;          # Energy [eV]
sigma = 3.4e-10;           # Length [m]


def convertp(x):
    return x*epsilon/(sigma**3)


def convertrho(x):
    return x/(sigma**3)


files = os.listdir('.')
pressurefiles = [i for i in files if 'pressures' in i]

rhostate = np.arange(0.1, 1.0, 0.1)
rhostate = [convertrho(i) for i in rhostate]

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

pstate = [convertp(i) for i in pstate]

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
    rho = convertrho(float(rho))
    rhos.append(rho)

    pressures = []
    with open(item) as file:
        for line in islice(file, 1, None):
            pressures.append(float(line))

    pvir.append(convertp(pressures[0]))
    pideal.append(convertp(pressures[1]))

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

ax.set_xlabel('Density [atoms/m^3]')
ax.set_ylabel('Pressure [eV/m^3]')
ax.grid()
ax.legend(loc='best')
fig.tight_layout()
fig.savefig('SIrho')
