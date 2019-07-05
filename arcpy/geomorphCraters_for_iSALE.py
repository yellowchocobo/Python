# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 15:49:09 2018

@author: nilscp
"""
import sys

pathm = ['/uio/kant/geo-ceed-u1/nilscp/Nils/Python/arcpy/']


# add directories to system path
for pat in pathm:
    sys.path.append(pat)



import re, os
import numpy as np
from scipy.optimize import curve_fit
from scipy import optimize
import calculate_volume_iSALE as cv 

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
******************************************************************************
'''

def leastsq_circle(x,y):
    
    # should only take not nan values
    
    # coordinates of the barycenter
    
    xt = x[~np.isnan(x)]
    yt = y[~np.isnan(y)]
    
    x_m = np.nanmean(xt)
    y_m = np.nanmean(yt)
    center_estimate = x_m, y_m
    
    center, ier = optimize.leastsq(f, center_estimate, args=(xt,yt))    
    xc, yc = center
    Ri       = calc_R(xt, yt, *center)
    R        = Ri.mean()
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu

'''
******************************************************************************
'''

def calc_R(x,y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

'''
******************************************************************************
''' 

def f(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()


'''
******************************************************************************
''' 

def power(x,a,b,c):
    
    return a * (x**b) + c

'''
******************************************************************************
''' 

def linear(x,a,b):
    y = a*x + b
    return y


'''
******************************************************************************
''' 

def tokenize(filename):
    '''
    Function to list filenames in correct order (see
    http://stackoverflow.com/questions/5997006/sort-a-list-of-files-using-python)

    :param filename:
    :return:
    '''
    digits = re.compile(r'(\d+)')

    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))
    
    
'''
******************************************************************************
'''

def calculation(path, paths):
    
    '''
    x : distance from the centre of the crater (obtained from crater.py and run through main.py)  
    z : elevation (obtained from crater.py and run through main.py)
    Rn : Rim-to-rim crater radius could either be imported or calculated on the fly

    path: path to where the XY profile of the final crater is located
    paths: path where the geomorphological parameters will be saved                                            
       
    
    '''    
    # Find the elevation or the distance closest to the rim-to-rim crater diameter
    # 1) can either import it
    # 2) or calculate it on the fly
    
    # 1)
    os.chdir(path)
    
    # define modelname
    first_part_string = path.split('/final')[0]
    modelname = first_part_string.split('data/')[1]
    
    # read cross section at the final crater diameter and derive data
    
    '''
    I have to update how a final crater diameter is assumed to have reached
    the final crater diameter. I should have the possibility to either fix and absolute
    time for a specific projectile diameter or specify manually or take a time normalized
    to the transient crater time (need in this case to be sure that the crater volume
    evolution is correct)
    '''
    # let's assume that all models have data for at least 4 times the transient crater time
    
    try:
        data = np.loadtxt(path + modelname + '_XYfinalprofile040.txt')
        
        dist = data[:,0]
        zi = data[:,1]
        cellsize = dist[1] - dist[0]
        
        
        '''
        This is also dependent on the fact that we get a reasonable estimate of the
        rim-to-rim crater diameter. In order to get robust images, I have to somehow
        modify crater/main.py
        '''
        # read the values at the final crater diameter
        dataf = np.loadtxt(path + modelname + '_final040.txt', comments="#")
        
        diam = dataf[-2]
        radius = diam /2.0
        
        value_nearest, idx_nearest = find_nearest(dist, radius)
        
        #distance normalized 
        dist_norm = dist/dist[idx_nearest]
                
        # find index
        A , idxA = find_nearest(dist_norm,0.0)
        B , idxB = find_nearest(dist_norm,0.1)
        C , idxC = find_nearest(dist_norm,0.7)
        D , idxD = find_nearest(dist_norm,0.8)
        E , idxE = find_nearest(dist_norm,0.9)
        F , idxF = find_nearest(dist_norm,1.0)
        G , idxG = find_nearest(dist_norm,1.2)
        H , idxH = find_nearest(dist_norm,2.0) # should be the maximum or end of the profile
        
        '''
        **************************************************************************
        '''
        # for upper cavity-wall radius of curvature
        # radius of circle fitted to the profile from D to F
        interval_upcw = zi[idxD:idxF+1]
        dist_upcw = dist[idxD:idxF+1]
        
        # we are interest in R_upcw
        #x_upcw, y_upcw, R_upcw, residu =leastsq_circle(dist_upcw,interval_upcw)
        try:
            __, __, R_upcw, __ =leastsq_circle(dist_upcw,interval_upcw)
            
        except:
            R_upcw = np.nan
        
        
        '''
        **************************************************************************
        '''
        # for upper flank radius of curvature
        interval_ufrc = zi[idxF:idxG+1]
        dist_ufrc = dist[idxF:idxG+1]
        
        #x_ufrc, y_ufrc, R_ufrc, residu =leastsq_circle(dist_ufrc,interval_ufrc)
        
        # new addition
        try:
            __, __, R_ufrc, __ =leastsq_circle(dist_ufrc,interval_ufrc)
            
        except:
            R_ufrc = np.nan
        
        '''
        **************************************************************************
        '''
        
        # cavity shape exponent
        interval_cse = zi[idxB:idxE+1]
        dist_cse = dist[idxB:idxE+1]
        
        try:
            a, b = curve_fit(power,dist_cse,interval_cse)
                
            exponent = a[1]
            cse = exponent
            
        except:
            cse = np.nan
        
        '''
        **************************************************************************
        '''
        
        # middle cavity wall slope angle
        #line fitted through
        interval_mcw = zi[idxC:idxE+1]
        dist_mcw = dist[idxC:idxE+1]
        
        try:
            a, b = curve_fit(linear, dist_mcw, interval_mcw)
            xs = np.linspace(np.min(dist_mcw),np.max(dist_mcw),100)
            ys = linear(xs,*a)
            
            #calculate the slope
            tetarad = np.arctan((ys[-1] - ys[0]) / (xs[-1] - xs[0]))
            slope_mcw = tetarad * (180./np.pi)
            
        except:
            slope_mcw = np.nan
        
        '''
        **************************************************************************
        '''
        
        # upper cavity wall slope angle
        interval_ucw = zi[idxD:idxF+1]
        dist_ucw = dist[idxD:idxF+1]
        
        try:
            a, b = curve_fit(linear, dist_ucw, interval_ucw)
            xs = np.linspace(np.min(dist_ucw),np.max(dist_ucw),100)
            ys = linear(xs,*a)
            
            #calculate the slope
            tetarad = np.arctan((ys[-1] - ys[0]) / (xs[-1] - xs[0]))
            slope_ucw = tetarad * (180./np.pi)
            
        except:
            slope_ucw = np.nan
        
        
        '''
        **************************************************************************
        '''
        
        # flank slope angle
        interval_fsa = zi[idxF:idxG+1]
        dist_fsa = dist[idxF:idxG+1]             
        
        try:
            a, b = curve_fit(linear, dist_fsa, interval_fsa)
            xs = np.linspace(np.min(dist_fsa),np.max(dist_fsa),100)
            ys = linear(xs,*a)
            
            #calculate the slope
            tetarad = np.arctan((ys[-1] - ys[0]) / (xs[-1] - xs[0]))
            slope_fsa = np.abs(tetarad * (180./np.pi))
            
            #upper and lower rim span
            slope_urs = 180. -  (slope_ucw + slope_fsa)
            slope_lrs =  180. - (slope_mcw + slope_fsa)
            
        except:
            slope_fsa = np.nan
            slope_urs = np.nan
            slope_lrs = np.nan
        
        '''
        **************************************************************************
        '''
    
        #average rim height
        h = zi[idxF]
    
        '''
        **************************************************************************
        '''        
        
        interval_hr = zi[idxF:idxH] 
        
        # minimum values beyond the rim of the crater
        
        try:
            min_h = np.nanmin(interval_hr)
            
            # height from the rim to the smallest elevation beyond the rim of the  crater 
            hr = zi[idxF] - min_h
            
        except:
            hr = np.nan
        
        '''
        **************************************************************************
        '''        
        # calculate the depth (new way where the min along each cross section is taken)
        depth_tmp = np.min(zi)
        depth = depth_tmp
        
        '''
        **************************************************************************
        '''
        
        # flank rim decay length
        interval_frdl = zi[idxF:idxH] # I think I can not use idxH+1 otherwise it is outside
        dist_frdl = dist[idxF:idxH] # The question is the rim flank going all the way up to 2 radius?
        
        try:
            
            frdlx1 = dist_frdl[:-1]
            frdlx2 = dist_frdl[1:]
            frdly1 = interval_frdl[:-1]
            frdly2 = interval_frdl[1:]
            
            dx = frdlx2 - frdlx1
            dy = frdly2 - frdly1
            
            tetarad = np.arctan(dy/dx)
            slope_frdl = np.abs(tetarad * (180./np.pi))
            
            # get the maximum slope
            slope_frdl_max = np.nanmax(slope_frdl)
            
            # where it is the closest of half the maximum
            __, idx_frdl = find_nearest(slope_frdl, slope_frdl_max/2.0)
            
            # get the distance at half the maximum
            frdl = dist_frdl[idx_frdl] - dist[idxF]
            
        except:
            frdl = np.nan
            
        '''
        **************************************************************************
        '''
            
        # cavity rim decay length
        interval_crdl = zi[idxA:idxF+1] # I think I can not use idxH+1 otherwise it is outside
        dist_crdl = dist[idxA:idxF+1]
        
        try:
            crdlx1 = dist_crdl[:-1]
            crdlx2 = dist_crdl[1:]
            crdly1 = interval_crdl[:-1]
            crdly2 = interval_crdl[1:]
            
            dx = crdlx2 - crdlx1
            dy = crdly2 - crdly1
            
            tetarad = np.arctan(dy/dx)
            slope_crdl = np.abs(tetarad * (180./np.pi))
            
            # get the maximum slope
            slope_crdl_max = np.nanmax(slope_crdl)
            
            # where it is the closest of half the maximum
            __, idx_crdl = find_nearest(slope_crdl, slope_crdl_max/2.0)
            
            # get the distance at half the maximum
            crdl = dist[idxF] - dist_crdl[idx_crdl]
            
            # volume gets calculate afterwards from calculate_volume.py
            
        except:
            crdl = np.nan
            
        # calculate the maximum depth    
        dD = (h-depth) / diam
    
            
        
        header_txt = ('mdiam;' + 
                      'mdepth;' + 
                      'mh;' +
                      'mhr;' +
                      'm_mcw;' +
                      'm_ucw;' +
                      'mcse;' +
                      'dD;' +
                      'mrupcw;' +
                      'mrufrc;' +
                      'mlrs;' +
                      'murs;' +
                      'mfsa;' +
                      'mcrdl;' +
                      'mfrdl;' +
                      'vol')
        
        # calculate the volume (NEED TO FIX)
        vol = cv.main(path) # I guess need to add 040 at the end of the script
        
        # fname
        fname = modelname + '_geomorph040.txt'
        
        # save the data to a text file (in a geomorphological parameters)
        np.savetxt(fname, np.column_stack([diam, depth, h, hr, slope_mcw, slope_ucw, cse, dD, R_upcw , R_ufrc, slope_lrs ,slope_urs ,slope_fsa ,crdl, frdl ,vol]), header=header_txt,
                           delimiter=';', fmt=['%1.6e', '%1.6e','%1.6e', '%1.6e','%1.6e', '%1.6e','%1.6e', '%1.6e',
                                               '%1.6e', '%1.6e','%1.6e', '%1.6e','%1.6e', '%1.6e','%1.6e', '%1.6e'])
        
    except:
        
        # if it does not work then set all values to nan
        R_upcw = R_ufrc = cse = slope_mcw = slope_ucw = slope_fsa = slope_lrs = slope_urs = crdl = frdl = h = hr = depth = diam = np.nan
    
    
    return (R_upcw, R_ufrc, cse, slope_mcw, slope_ucw, slope_fsa, slope_lrs, slope_urs, crdl, frdl,
            h, hr, depth, diam)



'''
******************************************************************************
'''

