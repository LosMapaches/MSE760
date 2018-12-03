from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as pl
import pandas as pd
import numpy as np
import math

from dataimport import importdata as data
from parmimport import importparm as parm

param = parm()
keys = list(param.keys())
runs = data(int(param[keys[0]]['hold']/param[keys[0]]['rate']))

for run in runs:
    runs[run] = runs[run].sort_values(by='c_mytemp', ascending=True)

    percentage = run.split('-')[-1]

    fig, sub = pl.subplots(1, 1)

    xdata = list(runs[run]['TimeStep'])
    xdata = [i/0.001 for i in xdata]

    ydata = list(runs[run]['c_mytemp'])
    sub.plot(xdata, ydata, 'b.', label='Data Points')

    sub.set_ylabel('Temperature [K] for '+percentage)
    sub.set_xlabel('Time [ps]')
    sub.legend(loc='upper left')
    sub.grid()
    fig.savefig('../figures/temperature_'+run)
