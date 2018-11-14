from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as pl
import pandas as pd
import numpy as np
import math

from dataimport import importdata as data


def function(E, T):
    kb = 8.6173303e-5  # [eV/K]
    return E-3.0*kb*T

runs, param = data(400)

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
    sub.plot(xdata, term, 'b.', label='Data Points')

    # Find bottom percent of data to fit
    length = len(xdata)
    start0 = 0
    end0 = math.floor(0.15*length)

    x0 = xdata[start0:end0]
    y0 = term[start0:end0]
    fit0 = np.polyfit(x0, y0, 1)
    fitf0 = np.poly1d(fit0)
    yfit0 = fitf0(xdata)

    fitrange0 = str([math.floor(xdata[start0]), math.floor(xdata[end0])])
    sub.plot(xdata, yfit0, 'r', label='Fit Range of '+fitrange0+' [K]')

    # Find top percent of data to fit
    start1 = math.ceil(0.15*length)
    end1 = -1
    x1 = xdata[start1:end1]
    y1 = term[start1:end1]
    fit1 = np.polyfit(x1, y1, 1)
    fitf1 = np.poly1d(fit1)
    yfit1 = fitf1(xdata)

    fitrange1 = str([math.floor(xdata[start1]), math.floor(xdata[end1])])
    sub.plot(xdata, yfit1, 'c', label='Fit Range of '+fitrange1+' [K]')

    try:
        npoints = 2
        point = [0]*(npoints+1)  # A point larger than 1 to start with
        tolerance = abs(max(term)-min(term))  # The initial tolerance

        # Iterate until point of intersection is found between fits
        while len(point) > npoints:
            point = np.argwhere(np.isclose(yfit0, yfit1, atol=tolerance)).reshape(-1)
            tolerance -= tolerance/10

        temps = np.array(xdata)[point]
        temp = math.floor(sum(temps)/len(temps))
        sub.axvline(x=temp, linestyle='--', label='Fit Intersection at '+str(temp)+' [K]')

    except Exception:
        pass

    sub.set_xlabel('Temperature [K] for '+percentage)
    sub.set_ylabel('E-3*k_b*T [eV/atom]')
    sub.legend(loc='upper left')
    sub.grid()
    fig.savefig('../figures/Energy_AlSm_slow_'+run)
