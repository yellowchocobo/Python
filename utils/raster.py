#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import rasterio as rio
import fiona
import rasterio.mask as mask

from rasterio.plot import reshape_as_raster, reshape_as_image
from rasterio.enums import Resampling

import matplotlib.pyplot as plt


'''
*****************************************************************************
'''

def read(in_raster):
    
    """
    
    Args:
        in_raster : absolute path to raster
        
    Returns:
        Numpy array
    """
    
    with rio.open(in_raster) as src:
        
        array = src.read()
        
    reshaped_array = reshape_as_image(array)
    
    
    return reshaped_array


'''
*****************************************************************************
'''

def resample(in_raster, resampling_factor, resampled_raster = ""):
    
    """
    
    Args:
        in_raster : absolute path to raster
        
    Returns:
        Numpy array
    """
    
    
    with rio.open(in_raster) as src:
        array = src.read(
            out_shape=(src.count, int(src.height * resampling_factor), 
                       int(src.width * resampling_factor)),
                       resampling=Resampling.cubic)
        
        if resampled_raster:
            profile = src.profile.copy()
            
            transform, width, height = rio.warp.calculate_default_transform(
                            src.crs, src.crs, int(src.width * resampling_factor), 
                            int(src.height * resampling_factor), *src.bounds)
            
            profile.update({
            'transform': transform,
            'width': width,
            'height': height})
    
            with rio.open(resampled_raster, 'w', **profile) as dst:
                for i in range(1, src.count + 1):
                    rio.warp.reproject(
                        source=rio.band(src, i),
                        destination=rio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=src.crs,
                        resampling=Resampling.cubic)
                
        else:
            None
            
    return reshape_as_image(array)


'''
*****************************************************************************
'''

def clip(in_raster, in_polygon, cliptype, clipped_raster = ""):
    
    """
    Args:
        in_raster: raster to be clipped
        in_polygon: as 
        cliptype : either 'bbox', 'shp', 'geojson'
        clipped_raster: raster to be saved if not it will .....
    """
    # shapefile
    if cliptype == "shp":

        with fiona.open(in_polygon, "r") as polygon:
            shapes = [feature["geometry"] for feature in polygon]
            
    # bounding box       
    elif cliptype == "bbox":
        None
    
    # if geojson polygon (e.g., in geopandas)    
    elif cliptype == "geojson":
        shapes = in_polygon
        
        
    with rio.open(in_raster) as src:
        out_image, out_transform = mask.mask(src, shapes, crop=True)
        out_meta = src.meta
            
    out_meta.update({"driver": "GTiff",
             "height": out_image.shape[1],
             "width": out_image.shape[2],
             "transform": out_transform})
    
    with rio.open(clipped_raster, "w", **out_meta) as dst:
        dst.write(out_image)
    
'''
*****************************************************************************
'''

def reproject(in_raster, dest_crs, reproj_raster):
    
    with rio.open(in_raster) as src:
        transform, width, height = rio.warp.calculate_default_transform(
            src.crs, dest_crs, src.width, src.height, *src.bounds)
        
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dest_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
    
        with rio.open(reproj_raster, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                rio.warp.reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dest_crs,
                    resampling=Resampling.cubic)    

'''
*****************************************************************************
'''

def extract_values_from_polygon(in_raster, bbox, clipped_raster = ""):

    None    

'''
*****************************************************************************
'''

#path = '/run/media/nilscp/pampa/ANALYSIS/SIMPLECRATERS_MOON/SLDEM_2015_LARGE_CRATERS/ascii'
#filename = os.path.join(path, 'crater0209.asc')
#
#data = read(filename)
#
#data2 = resample(filename, 2.0, './test.asc')
#
#in_raster = './test.asc'
#in_polygon = './layer.shp'
#clipped_raster = './clip.asc'
#clip_from_polygon(in_raster, in_polygon, clipped_raster)