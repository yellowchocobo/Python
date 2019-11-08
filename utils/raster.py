#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import rasterio as rio
import fiona
import rasterio.mask as mask
import geopandas as gpd
import json
import numpy as np

from rasterio.plot import reshape_as_raster, reshape_as_image
from rasterio.enums import Resampling
from shapely.geometry import Polygon, box
from itertools import product

import matplotlib.pyplot as plt

'''
*****************************************************************************
'''

def crs_eqc(crs_wkt_src, lat):
        
    # standard parallel should be replaced by the latitude
    crs_wkt_dst = crs_wkt_src.replace('["standard_parallel_1",0]', '["standard_parallel_1",' + str(int(lat)) + ']')
    
    # central meridian should be replaced by the longitude
    #crs_wkt = crs_wkt.replace('["central_meridian",0]', '["central_meridian",' + str(int(lon)) + ']')
    
    return crs_wkt_dst
'''
*****************************************************************************
'''

def parse_srid(src):    
    """
    Parse the SRID (EPSG code) from a raster open with rasterio    
    :param src: geodataframe with polygons    
    :returns: coordinate systems as EPSG code (integer)
    """
    return src.crs.to_epsg()

'''
*****************************************************************************
'''


def srid_to_wkt(crs_epsg):
    """Convert the SRID (EPSG) to GDAL-compatible projection metadata.
    The return value can be set directly on a data set with ds.SetProjetion()
    :param crs_epsg: coordinate systems as EPSG code (integer)  
    :returns: coordinate systems as wkt string
    """
    srs = rio.crs.CRS.from_epsg(crs_epsg)
    return srs.to_wkt()


'''
*****************************************************************************
'''

def boundary(src):
    """Get boundary of gdal tif.
    :param ds: rio dataset
    :returns: shapely polygon
    """
    ulx, xres, _, uly, _, yres = src.transform.to_gdal()
    lrx = ulx + (src.width * xres)
    lry = uly + (src.height * yres)
    
    bbox = [ulx, lry, lrx, uly]
    
    return bbox

'''
*****************************************************************************
'''

def boundary_to_polygon(bbox):
    
    [ulx, lry, lrx, uly] = bbox
    
    x_coords = [ulx, lrx, lrx, ulx, ulx] # will this be equivalent?
    y_coords = [uly, uly, lry, lry, uly]
    pol = Polygon(zip(x_coords, y_coords))
    
    return pol

'''
*****************************************************************************
'''

def get_extent(src, bbox):
    
    row_ul, col_ul = src.index(bbox[0], bbox[3])
    
    row_lr, col_lr = src.index(bbox[2], bbox[1])
    
    extent = (col_ul, 
              row_ul, 
              col_lr - col_ul,
              row_lr - row_ul)
    
    return extent

'''
*****************************************************************************
'''

def read(in_raster, pixelextent=None):
    """Read a raster.

    Args:
        filename: File to read
        pixelextent: (xoff, yoff, nx, ny)

    Returns:
        numpy array of shape (bands, ny, nx)
    """
    
    # open in_raster image
    with rio.open(in_raster) as src:
        
        # if pixel extent is provided, get raster data only for pixels
        if pixelextent:
            
            # if a rio.Windows.Window is provided
            if type(pixelextent) == rio.windows.Window:
                array = src.read(window=pixelextent)
                
            # if a list provided
            else:
                array = src.read(window=rio.windows.Window(*pixelextent))
        else:
            array = src.read()
      
    # reshape to (rows, columns, bands) from (bands, rows, columns)
    return reshape_as_image(array)




'''
*****************************************************************************
'''

def clip(in_raster, bbox, clipped_raster = ""):
    
    with rio.open(in_raster) as src:
        out_meta = src.meta
        
        # if we get a rio.windows.Windows directly
        if type(bbox) == rio.windows.Window:
            
            # read array for window
            array = src.read(window=bbox)
            
        else:
            # get window (can be done with src.index and indexes too)
            extent = get_extent(src, bbox)
                
            win = rio.windows.Window(*extent)
            
            # read array for window
            array = src.read(window=win)
        
        # shape of array
        dst_channel, dst_height, dst_width = np.shape(array)
        
        # get new transform
        win_transform = src.window_transform(win)
        
        # update meta information
        if clipped_raster:          
            out_meta.update({"driver": "GTiff",
                     "height": dst_height,
                     "width": dst_width,
                     "transform": win_transform})
            
            with rio.open(clipped_raster, "w", **out_meta) as dst:
                dst.write(array)
                
        else:
            None
            
        return array

'''
*****************************************************************************
'''

def clip_advanced(in_raster, in_polygon, cliptype, clipped_raster = ""):
    
    """
    Args:
        in_raster: raster to be clipped
        in_polygon: as 
        cliptype : either 'bbox', 'shp', 'geojson'
        clipped_raster: raster to be saved if not it will .....
    """
    
    with rio.open(in_raster) as src:
        out_meta = src.meta
    
    
        # get clip shape from polygon shape
        if cliptype == "shp":
            
            with fiona.open(in_polygon, "r") as polygon:
                shapes = [feature["geometry"] for feature in polygon]
                
        # get clip shape from bounding box       
        elif cliptype == "bbox":
            
            # convert it to a polygon bbox
            bbox_pol = box(in_polygon[0], in_polygon[1], in_polygon[2], in_polygon[3]) 
            
            # use geopandas to transform the bounding box into a json polygon
            geo = gpd.GeoDataFrame({'geometry': bbox_pol}, index=[0], crs=src.crs)
                
            # convert it to a json polygon
            shapes = [json.loads(geo.to_json())['features'][0]['geometry']]
        
        # or getting directly a geojson string (as above)   
        elif cliptype == "geojson":
            shapes = in_polygon
                    
        else:
            None
            # print error 
                        
        # clipping of raster
        out_image, out_transform = mask.mask(src, shapes, crop=True)
        
    # if clipped raster, we save the data to a tif (otherwise just return out_image)    
    if clipped_raster:          
        out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})
        
        with rio.open(clipped_raster, "w", **out_meta) as dst:
            dst.write(out_image)
            
    else:
        None
        
    return reshape_as_image(out_image)

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

def reproject(in_raster, dest_crs_wkt, reproj_raster):
    
    
    with rio.open(in_raster) as src:
        
        # get the new coordinate system
        dest_crs = src.crs.from_wkt(dest_crs_wkt)
        
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

def tile_windows(src, width = 512, height = 512):
    
    # number of columns and rows
    ncols, nrows = src.meta['width'], src.meta['height']
    
    offsets = product(range(0, ncols, width), range(0, nrows, height))
    
    big_window = rio.windows.Window(col_off=0, row_off=0, width=ncols, height=nrows)
    
    tile_window = []
    tile_transform = []
    
    for col_off, row_off in  offsets:
        window =rio.windows.Window(col_off=col_off, row_off=row_off, width=width, height=height).intersection(big_window)
        transform = rio.windows.transform(window, src.transform)
        
        tile_window.append(window)
        tile_transform.append(transform)
        
    return tile_window, tile_transform

'''
*****************************************************************************
'''

def tile_bounds(src, width = 512, height = 512):
    
    # number of columns and rows
    ncols, nrows = src.meta['width'], src.meta['height']
    
    offsets = product(range(0, ncols, width), range(0, nrows, height))
    
    big_window = rio.windows.Window(col_off=0, row_off=0, width=ncols, height=nrows)
    
    tile_bounds = []
    
    for col_off, row_off in  offsets:
        window =rio.windows.Window(col_off=col_off, row_off=row_off, width=width, height=height).intersection(big_window)
        win_transform = src.window_transform(window)
        
        # some problem here
        tile_bounds.append(rio.windows.bounds(window, win_transform))
        
    return tile_bounds

'''
*****************************************************************************
'''

def extract_values_from_polygon(in_raster, bbox, clipped_raster = ""):
    
    # TODO
    None    
    
'''
*****************************************************************************
'''

def vrt(in_raster, bbox, clipped_raster = ""):
    
    # TODO
    None   

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