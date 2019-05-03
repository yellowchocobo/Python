# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:27:24 2019

@author: nilscp
"""

'''
where to download:
    https://darts.isas.jaxa.jp/planet/pdap/selene/product_search.html

product:
    SLN-L-TC-5-DTM-MAP-SEAMLESS-V2.0 (3x3 degrees)
    SLN-L-TC-5-DTM-MAP-V2.0 (3x3 degrees)
    SLDEM (1x1 degrees)
    
    what are the  differences between them?
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-v2.0/ (SLN-L-TC-5-DTM-MAP-V2.0)
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/ (SLDEM2013)

ftp:
    /lon000/
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-seamless-v2.0/lon000/data/DTM_MAPs02_N00E000S03E003SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-seamless-v2.0/lon000/data/DTM_MAPs02_N03E000N00E003SC.lbl
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-seamless-v2.0/lon000/data/DTM_MAPs02_N06E000N03E003SC.img
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-seamless-v2.0/lon000/data/DTM_MAPs02_S03E000S06E003SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-dtm-map-seamless-v2.0/lon000/data/DTM_MAPs02_S06E000S09E003SC.img
    
    /lon003/ to /lon357/
    
    /lon357/
    
for SLDEM (1x1):
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N00E000S01E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N01E000N00E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N02E000N01E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_N90E000N89E001SC.lbl
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_S01E000S02E001SC.img
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon000/data/DTM_MAP_01_S02E000S03E001SC.img
    
    http://darts.isas.jaxa.jp/pub/pds3/sln-l-tc-5-sldem2013-v1.0/lon359/data/DTM_MAP_01_N00E359S01E000SC.img
    

to transfrom from .img to .tif:
    gdal_translate DTM_MAP_01_N03E060N02E061SC.lbl DTM_MAP_01_N03E060N02E061SC.tif
'''

import matplotlib.pyplot as plt
import numpy as np

# get from the label
bands = 1
lines = 12288
line_samples = 12288
data = np.fromfile(filename, dtype='>i2')
img = data.reshape(bands, lines, line_samples)
plt.imshow(img[0], cmap='gray')
plt.show()