# -*- coding: utf-8 -*-
"""
Created on May 02 2019

@author: Chocobo
"""
import os
import sys
import pandas as pd

sys.path.append('M:/Nils/Python/arcpy/') 

import scipy.spatial as ss
import numpy as np

'''
Problem what if, the detected rim is incomplete (then it will underestimate the 
volume of the crater)
'''

'''
**************************************************************************************************
'''
def tetrahedron_volume(a, b, c, d):
    
    '''
    from: https://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy
    '''
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

'''
**************************************************************************************************
'''

def convex_hull_volume_bis(pts):
    
    '''
    from: https://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy
    '''
    
    ch = ss.ConvexHull(pts)

    simplices = np.column_stack((np.repeat(ch.vertices[0], ch.nsimplex),
                                 ch.simplices))
    tets = ch.points[simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))
    
'''
**************************************************************************************************
'''

def load3d(filenameX, filenameY, filenameZ):

    dataX = np.loadtxt(filenameX,delimiter=";")
    dataY = np.loadtxt(filenameY,delimiter=";")
    dataZ = np.loadtxt(filenameZ,delimiter=";")
    
    return (dataX, dataY, dataZ)

'''
**************************************************************************************************
'''

def deg2rad(deg):
    
    output = deg * (np.pi/180.0)
    
    return output

'''
******************************************************************************
'''

def find_nearest(array, value):
    
    '''
    __, ixx = find_nearest(x, x2[0])
    __, iyy = find_nearest(y, y2[0])
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

'''
**************************************************************************************************
'''

def crossto3D(path, number_of_cross_sections):
    
    '''
    path same as paths
    '''
    
    os.chdir(path)
    
    # define modelname
    first_part_string = path.split('/final')[0]
    modelname = first_part_string.split('data/')[1]
    
    data = np.loadtxt(path + modelname + '_XYfinalprofile.txt')
    
    dist = data[:,0]
    zi = data[:,1]
    cellsize = dist[1] - dist[0]
    

    
    filenameX = modelname + "_cross_sectionsX.txt" 
    filenameY = modelname + "_cross_sectionsY.txt"
    filenameZ = modelname + "_cross_sectionsZ.txt"
    
    i = np.arange(0,360,360/number_of_cross_sections)
    
    
        # create empty arrays
    X = np.ones((len(i), len(zi)))
    Y = np.ones((len(i), len(zi)))
    Z = np.ones((len(i), len(zi)))
    
    for idx, angle in enumerate(i):
        
        x = np.cos(deg2rad(angle)) * dist
        y = np.sin(deg2rad(angle)) * dist
        
        idxi = int(idx)
        
        X[idxi,:] = x
        Y[idxi,:] = y
        Z[idxi,:] = zi
        
    np.savetxt(path + filenameX, X, delimiter = ";",fmt='%10.5f', comments='#')
    np.savetxt(path + filenameY, Y, delimiter = ";",fmt='%10.5f', comments='#')
    np.savetxt(path + filenameZ, Z, delimiter = ";",fmt='%10.5f', comments='#')
    
               
    print ("GENERIC CROSS SECTIONS have been calculated")
        
    

'''
**************************************************************************************************
'''

def crater_volume(path):
    
    '''
    Calculate the crater volume based on all available profiles
    
    '''
    
    
    os.chdir(path)
    
    # define modelname
    first_part_string = path.split('/final')[0]
    modelname = first_part_string.split('data/')[1]
    
    data = np.loadtxt(path + modelname + '_XYfinalprofile.txt')
    
    dist = data[:,0]
    
    # generate a "fake" 3D impact crater
    number_of_cross_sections = 15.0
    crossto3D(path, number_of_cross_sections)
    
        
    # load X, Y, Z for all profiles for specific crater_id
    filenameX = modelname + "_cross_sectionsX.txt"
    filenameY = modelname + "_cross_sectionsY.txt"
    filenameZ = modelname + "_cross_sectionsZ.txt"
    
    (dataX, dataY, dataZ) = load3d(filenameX, filenameY, filenameZ)
    
    # load the median diameter (MODIFY MAYBE from somewhere else) (most important to change)
    filename_meddiam = crater_id + "_res.txt"

    tmp_data = np.loadtxt(filename_meddiam, delimiter=";", comments="#")
    med_diam  = tmp_data[0]
    
    # number of cross sections 
    n = np.shape(dataX)[0]
    
    # transform it to a dictionnary
    data_dict = {'cross_section' : [], 'X' : [], 'Y' : [], 'elevation' : []} 
    
        
    # select only distances smaller than the median diameter
    
    try:
        idx = np.where(dist < (med_diam/2.0))
        argmax = np.max(idx)
        
        # set data in the newly created dictionnary
        for i in range(n):
            data_dict['cross_section'].append(i)  
            data_dict['X'].append(dataX[i,:argmax+1])
            data_dict['Y'].append(dataY[i, :argmax+1])
            data_dict['elevation'].append(dataZ[i, :argmax+1])
        
        # load dictionnary in pandas DataFrame   
        data_pd = pd.DataFrame(data_dict)
        
        nvalues = argmax + 1
        posit = np.ones((nvalues*n, 3))
        
        # prepare the data so we have an array with X, Y, Z (positions)
        for p in range(n):
            for i in range(nvalues):
                
                #index
                ind = (nvalues * p) + i
                
                posit[ind] = np.array([data_pd.X[p][i],data_pd.Y[p][i],data_pd.elevation[p][i]])
                
        # calculate the volume of the surface (result in m3)
        vol = convex_hull_volume_bis(posit)
        
    except:
        vol = np.nan
        
    return vol

'''
**************************************************************************************************
'''

def main(path, pathdata, filenamecrater):
    
    '''
    Calculate the crater volume for all craters contained in the filenamecrater
    file (e.g., crater0000, crater0001 ....)
    
    path = "X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_ASCII/"
    pathdata = "X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_DATA1D/"
    filenamecrater = "crater_id.txt"
    
    main(path, pathdata, filenamecrater)
    '''
    
    os.chdir(path)   
    crater_id = np.genfromtxt(filenamecrater,skip_header=1,dtype=str)
    
    
    vol = np.ones(len(crater_id))
    
    for indf, filename in enumerate(crater_id):
    
        print (indf)
        
        vol[indf]  = crater_volume(pathdata, filename)
        
        
    # saveit to 
    header_txt = ('volume')
    name_crater_v = filenamecrater.split(".txt")[0] + '_vol.txt'
    np.savetxt(pathdata + name_crater_v, np.column_stack(vol), delimiter = ";", header=header_txt,fmt='%10.5f', comments='#')
               
    return vol

'''
**************************************************************************************************
'''

#path = "X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_ASCII/"
#pathdata = "X:/Moon/ANALYSIS/SIMPLECRATERS_MOON/LINNE_DATA1D/"
#filenamecrater = "crater_id.txt"
    
#main(path, pathdata, filenamecrater)
