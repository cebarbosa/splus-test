# -*- coding: utf-8 -*-
"""

Created on 06/12/2016

@Author: Carlos Eduardo Barbosa

Calculation of the extinction coefficients

"""
import os
import fileinput

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import pymc
import corner

from config import *

def run_south():
    os.chdir(south_data_dir)
    bands = ("G", "I")
    catalogs = sorted([x for x in os.listdir(".") if x.startswith(bands) and
                x.endswith(".cat")])
    cols = {"G": (4,5,12,13,4,5,6,7,14), "I" : (8,9,12,13,4,5,6,7,14)}
    kappas = {"G" : 0.18, "I" : 0.08}
    for band in bands:
        cats = [x for x in catalogs if x.startswith(band)]
        data = np.loadtxt(fileinput.input(cats), usecols=cols[band])
        bounds = ((20, 0), (25, 1))
        sol = least_squares(res, [25, 0.5], args=(data, kappas[band]), bounds=bounds)
        print sol["x"]
        def model(data):
            M, Merr, m, merr, c1, c1err, c2, c2err, X = data.T
            zp = pymc.Uniform("zp", 0, 50)
            c = pymc.Uniform("c", -2, 2)
            @pymc.deterministic()
            def calib(zp=zp, c=c):
                return m + zp + c * (c1 - c2) - kappas[band] * X
            sigma = np.sqrt(merr**2 + Merr**2)
            y = pymc.Normal("y", mu=calib, tau=1 / sigma**2, value=M,
                            observed=True)
            return locals()
        M = pymc.MCMC(model(data))
        M.sample(5000, 1, 200)
        samples = np.array([M.trace("zp")[:], M.trace("c")[:]]).T
        tmp = corner.corner(samples[:,:])
        plt.show()
        raw_input()

def res(p, data, kappa):
    M, Merr, m, merr, c1, c1err, c2, c2err, X = data.T
    return (M - (m + p[0] + p[1] * (c1 - c2) - kappa * X)) / np.sqrt(Merr**2 + merr**2)

if __name__ == "__main__":
    run_south()