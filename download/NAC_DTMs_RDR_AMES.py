# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 11:14:49 2019

@author: nilscp

This script create text files to download all the high-resolution digital terrain
model for the Moon (DTM generated from NAC imagesS)
"""

import os
import urllib
from bs4 import BeautifulSoup
import re
import numpy as np

'''
**************************************************************************************************
'''

# several pages of data
html_page_link1 = ("http://wms.lroc.asu.edu/lroc/rdr_product_select?filter%5Btext%5D=&" +
            "filter%5Blat%5D=&filter%5Blon%5D=&filter%5Brad%5D=&filter%5Bwest%5D" +
            "=&filter%5Beast%5D=&filter%5Bsouth%5D=&filter%5Bnorth%5D=&filter%5B" +
            "prefix%5D%5B%5D=NAC_DTM&show_thumbs=0&per_page=100&commit=Search")

html_page_link2 = ("http://wms.lroc.asu.edu/lroc/rdr_product_select?commit=Search" + 
                   "&filter%5Blat%5D=&filter%5Blon%5D=&filter%5Bnorth%5D=&filter%" +
                   "5Bprefix%5D%5B%5D=NAC_DTM&filter%5Brad%5D=&filter%5Bsouth%5D=&" +
                   "filter%5Btext%5D=&filter%5Btopographic%5D=either&filter%5Bwest%" +
                   "5D=&page=2&per_page=100&show_thumbs=0&sort=time_reverse")

html_page_link3 = html_page_link2.split('&page=2')[0] + '&page=3' + html_page_link2.split('&page=2')[1]

html_page_link4 = html_page_link2.split('&page=2')[0] + '&page=4' + html_page_link2.split('&page=2')[1]

html_page_link5 = html_page_link2.split('&page=2')[0] + '&page=5' + html_page_link2.split('&page=2')[1]

'''
**************************************************************************************************
'''

# list containing all the pages
html_pages = [html_page_link1, html_page_link2, html_page_link3, html_page_link4, html_page_link5]
links = []


for html_page in html_pages:
    
    html_page_data = urllib.request.urlopen(html_page)
    
    soup = BeautifulSoup(html_page_data)
        
    for link in soup.findAll('a', attrs={'href': re.compile("^/lroc/view_rdr/")}):
        links.append('http://wms.lroc.asu.edu' + link.get('href'))
    
    
'''
**************************************************************************************************
'''
NACDTMlinks = []
DTMfolder = 'http://lroc.sese.asu.edu/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/DATA/SDP/NAC_DTM/'

#I need to loop through each of those links
for sublink in links:
    NACDTM_name = sublink.split('/')[-1]
    region_name = NACDTM_name.split('NAC_DTM_')[1]
    NACDTMlinks.append(DTMfolder + region_name + "/" + NACDTM_name + ".TIF")
    
    
    
'''
**************************************************************************************************
'''
html_DTM_AMES = 'http://lroc.sese.asu.edu/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/AMES_DTM/LRONAC_DTMS/'
html_DTM_AMES_data = urllib.request.urlopen(html_DTM_AMES)
soup_DTM_AMES = BeautifulSoup(html_DTM_AMES_data)
DTM_AMES = []


for link in soup_DTM_AMES.findAll('a', attrs={'href': re.compile("^NAC_DTM_M")}):
    DTM_AMES.append(link.get('href'))

DTM_AMES_FINAL = []

for DTM in DTM_AMES:
    if (DTM.endswith("_DEM.LBL") | DTM.endswith("_DEM.TIF")):
        DTM_AMES_FINAL.append(html_DTM_AMES+ DTM)
        
        
np.savetxt("NAC_DTM_AMES.txt", DTM_AMES_FINAL, fmt="%s")