#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:45:44 2019

@author: nilscp
"""

import numpy as np

'''
*****************************************************************************
'''

def crater_mineral_abundance(plg, ol, opx, cpx):
    
    '''
    plg : Plagioclase abundance in %
    ol: Olivine abundance in %
    opx: O
    
    
    
    Possible problems:
        What if, for some reasons the different rasters does not have the same
        resolution (maybe that is not a problem as we take the median values)
    '''
    
    # Plagioclase
    plg_median = np.median(plg)
    plg_25 = np.percentile(plg, 25)
    plg_75 = np.percentile(plg, 75)
    
    # Olivine
    ol_median = np.median(ol)
    ol_25 = np.percentile(ol, 25)
    ol_75 = np.percentile(ol, 75)
    
    # Orthopyroxene
    opx_median = np.median(opx)
    opx_25 = np.percentile(opx, 25)
    opx_75 = np.percentile(opx, 75)
    
    # Clinopyroxene
    cpx_median = np.median(cpx)
    cpx_25 = np.percentile(cpx, 25)
    cpx_75 = np.percentile(cpx, 75)
    
    # Pyroxene
    px = opx + cpx # need to have the same resolution
    px_median = np.median(px)
    px_25 = np.percentile(px, 25)
    px_75 = np.percentile(px, 75)
    
    return [plg_median, plg_25, plg_75,
            ol_median, ol_25, ol_75,
            opx_median, opx_25, opx_75,
            cpx_median, cpx_25, cpx_75,
            px_median, px_25, px_75]
    
    
'''
*****************************************************************************
'''
    

def classification(mineral_abundance):
    
    
    [plg_median, plg_25, plg_75,
            ol_median, ol_25, ol_75,
            opx_median, opx_25, opx_75,
            cpx_median, cpx_25, cpx_75,
            px_median, px_25, px_75] = mineral_abundance
     
     
    # first classification
    
    # Anorthosite
    if plg_median >= 0.90:        
         classv = 1
         
    elif np.logical_and(plg_median >= 0.775, plg_median < 0.90):
        
        if px_median > ol_median:
            classv = 2
        else:
            classv = 3
            
    elif np.logical_and(plg_median >= 0.60, plg_median < 0.775):
        if px_median > ol_median:
            classv = 4
        else:
            classv = 5
            
    elif np.logical_and(plg_median >= 0.10, plg_median < 0.60):
        if ol_median <= 0.10:
            classv = 6
        elif px_median <= 0.10:
            classv = 8
        else:
            classv = 7
    else:
        if px_median >= 0.90:
            classv = 9
        elif ol_median >= 0.90:
            classv = 11
        else:
            classv = 10
            
    # second classification
        
    if plg_median >= 0.90:        
         classv2 = 1
         
    elif np.logical_and(plg_median >= 0.775, plg_median < 0.90):
        
        if opx_median > cpx_median:
            classv2 = 2
        else:
            classv2 = 3
            
    elif np.logical_and(plg_median >= 0.60, plg_median < 0.775):
        if opx_median > cpx_median:
            classv2 = 4
        else:
            classv2 = 5
            
    elif np.logical_and(plg_median >= 0.10, plg_median < 0.60):
        if opx_median > cpx_median:
            
            if cpx_median <= 0.10:
                classv2 = 6
            else:
                classv2 = 7
        else:
            if opx_median <= 0.10:
                classv2 = 9
            else:
                classv2 = 8
    else:
        classv2 = 10

                        
    return (classv, classv2)
         
    
'''
*****************************************************************************
'''         
        