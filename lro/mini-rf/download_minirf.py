# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:25:01 2019

@author: nilscp

for Mini-RF Bistatic Radar Data:
    
    http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-2_3_5-bistatic-v2/lromrf_2xxx/data/rdr/
    
for Mini-RF Global Mosaics:
    
    http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-5-global-mosaic-v1/lromrf_1001/data/128ppd/
    
for all other datasets:
    
    http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/
    
    Volume	Data acquisition dates	Orbit range
    LROMRF_0001	2007-07-13 to 2010-01-19	200-2599
    LROMRF_0002	2010-01-19 to 2010-06-17	2600-4499
    LROMRF_0003	2010-06-17 to 2010-07-23	4500-4999
    LROMRF_0004	2010-07-23 to 2010-12-05	5000-6699
    LROMRF_0005	2010-12-05 to 2011-01-23	6700-7302
    
https://lunar.gsfc.nasa.gov/images/DataUsersWorkshop/Mini-RF.pdf

for 15m data?:
    
    grouped by 100 orbits:
    http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/06700_06799/level2/
    
    
https://isis.astrogeology.usgs.gov/fixit/projects/isis/wiki/Working_with_Lunar_Reconnaissance_Orbiter_MiniRF_Data

but need to get that from there



LROMRF_0004	2010-07-23 to 2010-12-05	5000-6699
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0004/data/sar/05000_05099/
-------------------------------
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0004/data/sar/06600_06699/

LROMRF_0005 2010-12-05 to 2011-01-23	6700-7302
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/06700_06799/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/06800_06899/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/06900_06999/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/07000_07099/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/07100_07199/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/07200_07299/level1/
http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_0005/data/sar/07300_07399/level1/


Need to search for data



In order to get some data for 
    
"""

'''
**************************************************************************************************
'''

import os
import numpy as np
import re
import urllib
import urllib2

'''
**************************************************************************************************
'''

def generateString():
    
    '''
    default_link = 'http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_'
    
    default_link_part2 = ['_0001/data/sar/', _0002/data/sar/, _0003/data/sar/, _0004/data/sar/, _0005/data/sar/]
    start_end_values = np.array(([200, 2599], [2600, 4499], [4500, 4999], [5000, 6699], [6700, 7399]))
    
    #need to create '07300_07399' type of strings
    
    '''

    default_link = 'http://pds-geosciences.wustl.edu/lro/lro-l-mrflro-4-cdr-v1/lromrf_'
    
    default_link_part2 = ['0001/data/sar/', '0002/data/sar/', '0003/data/sar/', 
                          '0004/data/sar/', '0005/data/sar/']
    
    start_end_values = np.array(([200, 2599], [2600, 4499], [4500, 4999], [5000, 6699], [6700, 7399]))
    
    
    default_link_part3 = []
    
    for i in range(len(default_link_part2)):
        
        # get start and end values
        start_value = start_end_values[i][0]
        end_value = start_end_values[i][1]
        
        # generate array from them
        start_values_array = np.arange(start_value, end_value, 100)
        end_values_array = start_values_array + 99
        
        default_link_tmp = []
        
        # create string array
        for j, v in enumerate(start_values_array):
        
            default_link_tmp.append(str(v).zfill(5) + '_' + str(end_values_array[j]).zfill(5))
        
        
        default_link_part3.append(default_link_tmp)
        
    # generate absolute path

    absolute_path = []
    
    for i in range(len(default_link_part2)):
        
        for j in range(len(default_link_part3[i])):
            absolute_path.append(default_link + default_link_part2[i] + default_link_part3[i][j] + '/level1/')
        
    return (absolute_path)

'''
**************************************************************************************************
''' 

def search(min_latitude, max_latitude):
    
    '''
    return
    '''
    
    # get all the absolute paths    
    (absolute_path) = generateString()
    
    # create empty list of array
    download_list = []
    
    # loop through folders that containt data and get all raw products
    for path in absolute_path:
        
        try:
            response = urllib2.urlopen(path)
            html = response.read()
        
            data = html.split('</URL>')
                    
            #
            content_webpage = data[0]
            lines = content_webpage.split('</A><br>')
            
            # split
            for line in lines:
                if line.endswith('.img'):
                    download_list.append(path + line.split('">')[1])
                    
        except:
            None
            
    
    # take only values that are contains or equal to min and max latitude
    array_latitudes = range(min_latitude, max_latitude + 1)
    array_latitudes_string = []
    
    for latitude in array_latitudes:
        if latitude >= 0:
            array_latitudes_string.append(str(latitude) + 'n')
        else:
            array_latitudes_string.append(str(abs(latitude)) + 's')
    
    download_selected = []
    
    for string in array_latitudes_string:
    
        pattern = re.compile(string)
        
        for product in download_list:
            if pattern.search(product):
                download_selected.append(product)
    
    

        
    # download_list should contain all data
        
    txt = []
    lbl = []
        
    for product in download_selected:
        txt.append(product.split(".img")[0] + ".txt")
        lbl.append(product.split(".img")[0] + ".lbl")
        
    global_list = download_selected + txt + lbl
    
    return global_list

'''
**************************************************************************************************
'''

min_latitude = 27
max_latitude = 30
global_list = search(min_latitude, max_latitude)
paths = "/uio/kant/geo-ceed-u1/nilscp/Desktop/astra/TMP_DOWNLOAD/MINI-RF/N30N27/"
np.savetxt(paths + "download_minirf.txt", global_list, fmt='%s')
                
'''
**************************************************************************************************
'''    

 