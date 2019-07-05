# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:43:40 2019

@author: nilscp
"""

import os
import numpy as np

path = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/double_detrending/data/"
patha = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/ascii4R/"
filename = "crater0173XY.txt"


data = np.loadtxt(path + filename, delimiter=";", comments="#", skiprows=1)


X = data[:,0] 
Y = data[:,1]

import matplotlib.pyplot as plt
plt.plot(X,Y)

def readheader(path, filename):
    
    '''
    (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = readheader(filename)
    
    OK: 18.10.2018
    '''
    
    os.chdir(path)
    
    lines =  []
    with open(filename) as f:
        ix = 0
        for line in f:
            if ix < 6:
                lines.append(line)
                ix = ix + 1
            else:
                break
    
    for ix, line in enumerate(lines):
        if ix == 0:
            tmp = line.strip('\n')
            ncols = int(tmp.split('ncols')[1])
        elif ix == 1:
            tmp = line.strip('\n')
            nrows = int(tmp.split('nrows')[1])
        elif ix == 2:
            tmp = line.strip('\n')
            xllcorner = float(tmp.split('xllcorner')[1])            
        elif ix == 3:
            tmp = line.strip('\n')
            yllcorner = float(tmp.split('yllcorner')[1])
        elif ix == 4:
            tmp = line.strip('\n')
            cellsize = float(tmp.split('cellsize')[1])
        else:
            tmp = line.strip('\n')
            NODATA_value = float(tmp.split('NODATA_value')[1])            
            
    return (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)


(ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = readheader(patha, "crater0173.asc")

XX = xllcorner + Y
YY = yllcorner + X

allv = np.column_stack((XX,YY))
paths = "D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_LARGE_CRATERS/steps_algorithm/"

np.savetxt(paths + "test2.txt", allv, delimiter=";")