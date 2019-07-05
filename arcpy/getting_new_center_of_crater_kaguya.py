# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:39:55 2019

@author: nilscp
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
import sys
import richdem as rd

#sys.path.append('M:/Nils/Python/arcpy/')
sys.path.append('/uio/kant/geo-ceed-u1/nilscp/Nils/Python/arcpy/')

import geomorphCraters as wk

xcoords = []
ycoords = []
'''
**************************************************************************
'''

def onclick(event):
    
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print ('x = %d, y = %d' % (ix, iy))

'''
**************************************************************************
'''

def readASCII(name_ascii, scaling_factor = 1.0):
    
    '''
    need to multiply by 0.5 for SLDEM2015
    '''
    
    data = pd.read_csv(name_ascii, delimiter=' ', skiprows=6).values
    data2 = np.rot90(data)
    data3 = np.rot90(data2)
    data4 = np.rot90(data3) # need to divide by 2 for LOLA + Kaguya = SLDEM2015
    
    del data
    del data2
    del data3
    data = copy.deepcopy(data4) * scaling_factor #multiplying factor of the elevation (not for kaguya, 0.5 for Kaguya + LOLA)
    del data4
    
    return data

'''
**************************************************************************
'''

def constrainASCII(pathascii, filename, data):
    
    (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value) = wk.readheader(pathascii,filename)

    # take the minimum of ncols, nrows and ncols and nrows of data
    nrowsd, ncolsd = np.shape(data)
    ncolsnrows = np.min([ncols, nrows, ncolsd, nrowsd])
    data = data[:ncolsnrows, :ncolsnrows]
    print (ncols, nrows, ncolsd, nrowsd)
    
    x = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
    y = np.linspace(0,(ncolsnrows-1)*cellsize,ncolsnrows)
    
    xe = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
    ye = np.linspace(cellsize/2.0, (cellsize/2.0) + ((ncolsnrows-1)*cellsize), ncolsnrows)
    
    # create my own matrices for the plotting with pcolor or pcolormesh
    xc = np.zeros_like(data)
    yc = np.zeros_like(data)
    xce = np.zeros_like(data)
    yce = np.zeros_like(data)
    
    for i in range(ncolsnrows):
        xc[:,i] = x
        xce[:,i] = xe
        
    for i in range(ncolsnrows):
        yc[i,:] = y
        yce[i,:] = ye
    
            
    # middle of the dtm (this has to change!)
    ncentery = int(nrows/2) # changed
    ncenterx = int(ncols/2) # changed
    
    return (xc, yc, data, ncenterx, ncentery)

'''
**************************************************************************
'''    

def loadfigures(pathascii, pathvis16R, pathxy, filename):
    
    '''
    pathascii = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii/'
    pathvis8R = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_visible8R/'
    pathvis32R = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_visible/'
    pathxy = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/data/'
    
    
    
    pathascii = '/uio/kant/geo-ceed-u1/nilscp/Desktop/astra/TMP_DOWNLOAD/ascii/'
    pathvis8R = '/uio/kant/geo-ceed-u1/nilscp/Desktop/astra/TMP_DOWNLOAD/ascii_visible8R/'
    pathvis32R = '/uio/kant/geo-ceed-u1/nilscp/Desktop/astra/TMP_DOWNLOAD/ascii_visible/'
    pathxy = 'D:/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/data/'
    filename = 'cpcrater0000'
    '''
    name_ascii = pathascii + filename + '.asc'
    name_crater_txt = pathxy + filename + 'XY.txt'
    name_vis16R = pathvis16R + filename + '_visible.asc'
    
    # should be good this way
    data = readASCII(name_ascii)
    datavis16R = readASCII(name_vis16R)
    
    # load data and constrain size of the array
    (xc, yc, data, ncenterx, ncentery) = constrainASCII(pathascii, name_ascii, data)
    (xcv1, ycv1, datav1, ncenterxv1, ncenteryv1) = constrainASCII(pathvis16R, name_vis16R, datavis16R)
    
    # load xy
    dataxy = np.loadtxt(name_crater_txt, skiprows=1, delimiter=";")
    datax = dataxy[:,0]
    datay = dataxy[:,1]
    
    # slope
    datareload = rd.rdarray(data, no_data=-9999)
    slope = rd.TerrainAttribute(datareload, attrib='slope_riserun')
    
    # profile_curvatyre
    pfc = rd.TerrainAttribute(datareload, attrib='profile_curvature')
    
    # planform curvature
    pfc2 = rd.TerrainAttribute(datareload, attrib='planform_curvature')
    
    # curvature
    pfc3 = rd.TerrainAttribute(datareload, attrib='curvature')
    # I could calculate the slope and the curvature and then plot it
    
    
    # plot figures (1) visible 32R, (1) visible 8R, (2) DTM 8R, (3) Slope 8R, (4) Curvature 8R
    # when (4) is open, engage 
    
    fig2 = plt.figure(1)
    plt.pcolormesh(xcv1, ycv1, datav1)
    plt.colorbar()
    
    fig3 = plt.figure(2)
    plt.pcolormesh(xc, yc, data)
    plt.colorbar()
    
    fig4 = plt.figure(3)
    plt.pcolormesh(xc, yc, slope)
    plt.colorbar()
    
    return fig4
        

    
    
'''
**************************************************************************
'''

pathascii = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_16R/'
pathvis16R = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/ascii_visible16R/'
pathxy = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2013_COLDSPOTS/double_detrending/data/'

ixx = ixx + 1
filename = 'cpcrater0'+ str(ixx)
fig4 = loadfigures(pathascii, pathvis16R, pathxy, filename)

cid = fig4.canvas.mpl_connect('button_press_event', onclick)
print (filename)

fig4.canvas.mpl_disconnect(cid)
for i in range(3):
    plt.close()