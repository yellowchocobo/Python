# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:25:01 2019

@author: nilscp

for SLDEM (1x1):
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N00E000S01E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N01E000N00E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N02E000N01E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N90E000N89E001SC.lbl
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_S01E000S02E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_S02E000S03E001SC.img
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon359/data/DTM_MAP_01_N00E359S01E000SC.img
    

"""

'''
**************************************************************************************************
'''

import os
import numpy as np

'''
**************************************************************************************************
'''

def download(path_to_save, default_link, default_link2, min_latitude, max_latitude, min_longitude, max_longitude, degrees_per_latlon = 1):
    
    '''
    generate links to download SLDEM2013 (1x1 degree)
    
    default_link = 'http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-ortho-map-v2.0/' # (images)
    default_link = 'http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/' # (dtm)
    
    default_link2 = '/data/TCO_MAP_02_'
    default_link2 = '/data/DTM_MAP_01_'
    '''    
    img_to_download = []
    lbl_to_download = []
    
    #default_link = 'http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/'
    
    for i in range(min_longitude, max_longitude, degrees_per_latlon):
        
        lon1 = 'lon' + str(int(i)).zfill(3)
        lon2 = 'E' + str(int(i)).zfill(3)
        
        if i + degrees_per_latlon == 360:
            #lon3 = 'E' + str(int(0)).zfill(3) # for DTM_MAP
            lon3 = 'E' + str(int(i+degrees_per_latlon)).zfill(3)
        else:    
            lon3 = 'E' + str(int(i+degrees_per_latlon)).zfill(3)
        
        '''
        Must modify if we want data all the way up to 90 degrees
        '''
        for j in range(min_latitude, max_latitude, degrees_per_latlon):
            
            if j  < 0:
                lat1 = 'S' + str(int(abs(j))).zfill(2)
                lat2 = 'S' + str(int(abs(j-degrees_per_latlon))).zfill(2)
                
            elif j == 0:
                lat1 = 'N' + str(int(abs(j))).zfill(2)
                lat2 = 'S' + str(int(abs(j-degrees_per_latlon))).zfill(2)
                
            else:
                lat1 = 'N' + str(int(abs(j))).zfill(2)
                lat2 = 'N' + str(int(abs(j-degrees_per_latlon))).zfill(2)

            # do
            img_to_download.append(default_link + lon1 + default_link2 + 
                                   lat1 + lon2 + lat2 + lon3 + 'SC.img')            
            
            lbl_to_download.append(default_link + lon1 + default_link2 + 
                                   lat1 + lon2 + lat2 + lon3 + 'SC.lbl')
            
            
    np.savetxt(path_to_save + 'download.txt', np.array(img_to_download + lbl_to_download), fmt="%s")
    print ('DONE')
         
'''
**************************************************************************************************
'''

def select_DTMs(bottom, top, left, right):
    
    '''
    equivalent to bottom, top, left, right
    
    
    
    '''    
    dtm_selection = []
        
    min_latitude = np.int(np.floor(bottom))
    max_latitude = np.int(np.ceil(top))
    min_longitude = np.int(np.floor(left))
    max_longitude = np.int(np.ceil(right))
    
    for i in range(min_longitude, max_longitude):
        
        lon2 = 'E' + str(int(i)).zfill(3)
        
        if i + 1 == 360:
            lon3 = 'E' + str(int(0)).zfill(3)
        else:    
            lon3 = 'E' + str(int(i+1)).zfill(3)
        
        '''
        Must modify if we want data all the way up to 90 degrees
        '''
        for j in range(min_latitude, max_latitude):
            
            if j  < 0:
                lat1 = 'S' + str(int(abs(j+1))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            elif j == 0:
                lat1 = 'N' + str(int(abs(j+1))).zfill(2)
                lat2 = 'S' + str(int(abs(j))).zfill(2)
                
            else:
                lat1 = 'N' + str(int(abs(j+1))).zfill(2)
                lat2 = 'N' + str(int(abs(j))).zfill(2)

            # do
            dtm_selection.append('DTM_MAP_01_' + 
                                   lat1 + lon2 + lat2 + lon3 + 'SC.tif')            
            
    return (dtm_selection)
         
'''
**************************************************************************************************
'''

path_to_save = 'D:/Kaguya/ORTHOIMAGES/'
default_link = 'http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-ortho-map-v2.0/'
default_link2 = '/data/TCO_MAP_02_'
min_latitude = -60 # 60
max_latitude = 62 # 60
min_longitude = 0
max_longitude = 360 # actually 360
 
download(path_to_save, default_link, default_link2, min_latitude, max_latitude, min_longitude, max_longitude, 3)    