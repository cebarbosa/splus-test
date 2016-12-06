# -*- coding: utf-8 -*-
"""

Created on 06/12/2016

@Author: Carlos Eduardo Barbosa

Calculation of the extinction coefficients

"""
import os
import fileinput

import numpy as np

from config import *

def run_south():
    os.chdir(south_data_dir)
    bands = ("G", "I")
    catalogs = sorted([x for x in os.listdir(".") if x.startswith(bands) and
                x.endswith(".cat")])
    cols = {"G": (4,5,12,13,4,5,6,7,14), "I" : (8,9,12,13,4,5,6,7,14)}
    for band in bands:
        cats = [x for x in catalogs if x.startswith(band)]
        data = np.loadtxt(fileinput.input(cats), usecols=cols[band])
        print data.shape


if __name__ == "__main__":
    run_south()