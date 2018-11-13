from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as pl
import pandas as pd
import numpy as np

from dataimport import importdata as data


def function(E, T):
    kb = 8.6173303e-5  # [eV/K]
    return E-3.0*kb*T

runs, param = data(6)

for run in runs:
    runs[run] = runs[run].sort_values(by='c_mytemp', ascending=True)


    energy = runs[run]['c_PE']+runs[run]['c_KE']
    energy = [i/param[run]['atoms'] for i in energy]
    term = [function(i, j) for i, j in zip(energy, runs[run]['c_mytemp'])]

    percentage = run.split('percent')[0]

    percentage = percentage.split('p')
    percentage = percentage[0]+'.'+percentage[1]+' % Sm'

    fig, sub = pl.subplots(1, 1)
    xdata = list(runs[run]['c_mytemp'])
    sub.plot(xdata, term, '.', label='Data Points')

    try:
        npoints = 100
        kterm = 4
        xfit = np.linspace(min(xdata), max(xdata), npoints)
        spl = UnivariateSpline(xdata, term, k=kterm)
        spl.set_smoothing_factor(1e4)
        sub.plot(xfit, spl(xfit), label='Spline Fit (k='+str(kterm)+')')

    except Exception:
        pass

    sub.set_xlabel('Temperature [K] for '+percentage)
    sub.set_ylabel('E-3*k_b*T [eV/atom]')
    sub.legend(loc='best')
    sub.grid()

pl.show()
