# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:21:38 2019

@author: nilscp
"""
import glob, os
from osgeo import gdal
import subprocess

'''
**************************************************************************************************
'''
def convert_sc360_to_sc180(ulx_scy_360, lrx_scy_360, cellsize = 7.403161724669900):
    
    ulx_degree = ulx_scy_360 / (4096*cellsize)
    
    lrx_degree = lrx_scy_360 / (4096*cellsize)
    
    if ulx_degree >= 180.0:
        ulx_degree_180 = -(360.0 - ulx_degree)
    else:
        ulx_degree_180 = ulx_degree
        
    if lrx_degree >= 180.0:
        lrx_degree_180 = -(360.0 - lrx_degree)
    else:
        lrx_degree_180 = lrx_degree
        
    ulx_sc180 = ulx_degree_180 * (4096*cellsize)
    
    lrx_sc180 = lrx_degree_180 * (4096*cellsize)
    
    return (ulx_sc180, lrx_sc180)
    

'''
**************************************************************************************************
'''

def getExtentRasters(path_raster, list_tif):
    
    
    # change directory
    os.chdir(path_raster)
    
    for raster in list_tif:
        
        print (raster)
    
        ds=gdal.Open(raster)
    
        gt=ds.GetGeoTransform()
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        ext=GetExtent(gt,cols,rows)
        
        ulx = ext[0][0]
        uly = ext[0][1]
        lrx = ext[2][0]
        lry = ext[2][1]
    
        (ulx_sc180, lrx_sc180) = convert_sc360_to_sc180(ulx, lrx, cellsize = 7.403161724669900)
        
        nname = raster.split(".tif")[0] + "_fix.tif"
        
        translate(['-a_ullr', str(ulx_sc180), str(uly), str(lrx_sc180), str(lry), path_raster + raster, path_raster + nname])
        
'''
**************************************************************************************************
'''

def translate(args):
    
    """with a def you can easily change your subprocess call"""
    # command construction with binary and options
    options = ['/work/nilscp/anaconda2/bin/gdal_translate']
    options.extend(args)
    # call gdalwarp 
    subprocess.check_call(options)

'''
**************************************************************************************************
'''


def GetExtent(gt,cols,rows):
    ''' Return list of corner coordinates from a geotransform

        @type gt:   C{tuple/list}
        @param gt: geotransform
        @type cols:   C{int}
        @param cols: number of columns in the dataset
        @type rows:   C{int}
        @param rows: number of rows in the dataset
        @rtype:    C{[float,...,float]}
        @return:   coordinates of each corner
    
    raster=r'DTM_MAP_01_S45E350S46E351SC.tif'
    ds=gdal.Open(raster)

    gt=ds.GetGeoTransform()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    ext=GetExtent(gt,cols,rows)
    
    ulx = ext[0][0]
    uly = ext[0][1]
    lrx = ext[2][0]
    lry = ext[2][1]
    
    (ulx_sc180, lrx_sc180) = convert_sc360_to_sc180(ulx, lrx, cellsize = 7.403161724669900)
    '''
    ext=[]
    xarr=[0,cols]
    yarr=[0,rows]

    for px in xarr:
        for py in yarr:
            x=gt[0]+(px*gt[1])+(py*gt[2])
            y=gt[3]+(px*gt[4])+(py*gt[5])
            ext.append([x,y])
            print (x,y)
        yarr.reverse()
    return ext

'''
**************************************************************************************************
'''
path_raster = '/run/media/nilscp/pampa/Kaguya/SLDEM2013/'

os.chdir(path_raster)

list_tif = glob.glob("*.tif")

getExtentRasters(path_raster, list_tif)