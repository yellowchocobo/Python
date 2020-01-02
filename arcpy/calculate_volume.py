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

def crater_volume(pathdata, crater_id, resolution_DEM = 59.2252937999999):
    
    '''
    Calculate the crater volume based on all available profiles
    
    '''
    
    os.chdir(pathdata)
    
    # second try
    try:
    
        # load X, Y, Z for all profiles for specific crater_id
        filenameX = crater_id + "_cross_sectionsX.txt"
        filenameY = crater_id + "_cross_sectionsY.txt"
        filenameZ = crater_id + "_cross_sectionsZ.txt"
        
        (dataX, dataY, dataZ) = load3d(filenameX, filenameY, filenameZ)
        
        # load the median diameter
        filename_meddiam = crater_id + "_res.txt"
    
        tmp_data = np.loadtxt(filename_meddiam, delimiter=";", comments="#")
        med_diam  = tmp_data[0]
        
        # number of cross sections 
        n = np.shape(dataX)[1]
        
        # transform it to a dictionnary
        data_dict = {'cross_section' : [], 'X' : [], 'Y' : [], 'elevation' : []} 
        
        
        # calculate the distance along the profiles (same for all)
        dist_cells = np.sqrt(((dataX[:,0] - dataX[0,0])**2.) + ((dataY[:,0] - dataY[0,0])**2.))
        
        # SLDEM resolution (does this take into account the sampling at half the SLDEM resolution?)
        # the distance should only be half
        # dist = dist_cells * (resolution_DEM/2.0)
        dist = dist_cells * resolution_DEM 
        
        # select only distances smaller than the median diameter
        
        try:
            idx = np.where(dist < (med_diam/2.0))
            argmax = np.max(idx)
            
            # set data in the newly created dictionnary
            for i in range(n):
                data_dict['cross_section'].append(i)  
                data_dict['X'].append(dataX[:argmax+1,i])
                data_dict['Y'].append(dataY[:argmax+1,i])
                data_dict['elevation'].append(dataZ[:argmax+1,i])
            
            # load dictionnary in pandas DataFrame   
            data_pd = pd.DataFrame(data_dict)
            
            nvalues = argmax + 1
            nprof = np.shape(dataX)[1]
            posit = np.ones((nvalues*nprof, 3))
            
            # prepare the data so we have an array with X, Y, Z (positions)
            for p in range(nprof):
                for i in range(nvalues):
                    
                    #index
                    ind = (nvalues * p) + i
                    
                    posit[ind] = np.array([data_pd.X[p][i]*resolution_DEM,data_pd.Y[p][i]*resolution_DEM,data_pd.elevation[p][i]])
                    
            # calculate the volume of the surface (result in m3)
            vol = convex_hull_volume_bis(posit)
            
        except:
            vol = np.nan
            
    except:
        vol = np.nan
        
    return vol

'''
**************************************************************************************************
'''

def main(path, pathdata, filenamecrater, resolution_DEM):
    
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
        
        vol[indf]  = crater_volume(pathdata, filename, resolution_DEM)
        
        
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
