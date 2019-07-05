# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:55:21 2019

@author: nilscp
"""

import numpy as np
import os

path = "D:/GeologyMoon/WATTERS/"
filename = "watters.csv"

data = np.genfromtxt(path + filename, delimiter=",", skip_header=1)

im = data[:,-2]
cav = data[:,13]
depth = data[:,9]
diam = data[:,7]
dD = depth/diam


ixfresh = np.where(im <= 0.0)
ixold = np.where(im > 3.0)

import scipy.stats
scipy.stats.ks_2samp(dD[ixfresh], dD[ixold])
print (np.nanmedian(dD[ixfresh]), np.nanstd(dD[ixfresh]))
print (np.nanmedian(dD[ixold]), np.nanstd(dD[ixold]))