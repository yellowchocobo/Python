# -*- coding: utf-8 -*-

'''
******************************************************************************

     =========================================================================
     Subroutine to read, post-process and save crater dimensions. Call subroutines
     such as ejecta.py and crater.py

     Description                                     Programmer    Date
     ------------------------------------------------------------------
     Original version (1.0).............................NCP  2016/01/12
     Modified version (2.0).............................NCP  2017/11/23
     xxx                      ..........................NCP  xxxx/xx/xx
     xxx                      ..........................NCP  xxxx/xx/xx
    
    The version 2.0 includes:
    - a better description of functions
    - changed the name of some function 
    - functions:
    
    
    ==========================================================================

******************************************************************************
'''

# import Python's module
import sys
import glob
import os


pathm = ['/uio/kant/geo-ceed-u1/nilscp/Nils/Python/arcpy/', '/uio/kant/geo-ceed-u1/nilscp/Nils/Python/craterGeomorph/'] # change to the path where you have the scripts


# add directories to system path
for pat in pathm:
    sys.path.append(pat)


import crater as cr
import plot as iplt


# path to in-house scripts and pySALEPlot
path_pySALEPlot = '/work/nilscp/iSALE/Dellen/lib' # change to the path where you have the Dellen install

# path to pySALEPlot is loaded
sys.path.append(path_pySALEPlot)
import pySALEPlot as psp


'''
***********************************************************************
'''

# INPUT (2, 3 and 4)
# related to how the depth of the crater is calculated (1 is a good value)
mode = 1
poros = 1  # if porous = 1, non porous = 0
id_mat = 0  # 0 = everything, 2 only lower material, especially practical in layering modelling
transient = 0  # 0 transient at max volume, 1: transient at maximum depth
thresholdf = 0.01
g = 3.711

'''
********************CALCULATE DIMENSION EVOLUTION******************************

It might be some problems about the evolution of the crater volume through time
depending on which model is used 
'''


# path containing all the simulations (ROCK_AVG2_10 and SAND_AVG2_10)
path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/' # change to the path where you have the folders with jdata in it


os.chdir(path)

# selection of folders
search = '*AVG*'
models = glob.glob(search)

pathdata = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/' # create data folder


###############################################################################
# 2. CALCULATE CRATER DIMENSIONS DURING CRATER EVOLUTION
# 3. CALCULATE CRATER DIMENSIONS AT THE TRANSIENT CRATER
# 4. CALCULATE CRATER DIMENSIONS AT THE FINAL CRATER (LAST FIVE TIMESTEPS)
###############################################################################

cr.main(path, models, pathdata, mode, id_mat, transient)


###############################################################################
# 5. SAVE TXT FILE WITH X AND Y FOR THE TRANSIENT and FINAL CRATER
###############################################################################
cr.craterProfiles(models, path, pathdata, 0, 1)  # XY transient craters
cr.craterProfiles(models, path, pathdata, 0, 2)  # XY final craters


'''
****************PLOTS ROCK_AVG2_10********************
''' 

path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/ROCK_AVG2_10/'
paths = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/ROCK_AVG2_10/'
normpath = ''
norm = 0
zoom_id = ['close', 'mid', 'hires', 'all', 'manual'] # 'mid'; 'hires'; 'all' #depending on the zoom
param = 'Den'
vmiin = [] # or vmiin = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e-5, 0.0] # need to take into account the field factor
vmaax = [] # or vmaax = [1.0, 1.0, 2000.0, 400.0, 1.0, 1.0, 1.0, 1.0, 1.25,1.0, 1.0]
manualx = ((0,1500)) #distance in m
manualy = ((-600,500)) # distance in m

lbl = r"$\mathit{rock, L} = 100 m"

iplt.field(path1, paths, norm, normpath, zoom_id, lbl, vmiin, vmaax, manualx, manualy, param)


'''
****************PLOTS SAND_AVG2_10********************
''' 

path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/SAND_AVG2_10/'
paths = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/SAND_AVG2_10/'
normpath = ''
norm = 0
zoom_id = ['close', 'mid', 'hires', 'all', 'manual'] # 'mid'; 'hires'; 'all' #depending on the zoom
param = 'Den'
vmiin = [] # or vmiin = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e-5, 0.0] # need to take into account the field factor
vmaax = [] # or vmaax = [1.0, 1.0, 2000.0, 400.0, 1.0, 1.0, 1.0, 1.0, 1.25,1.0, 1.0]
manualx = ((0,1500))
manualy = ((-600,500)) # distance in m

lbl = r"$\mathit{sand, L} = 100 m"

iplt.field(path1, paths, norm, normpath, zoom_id, lbl, vmiin, vmaax, manualx, manualy, param)


'''
****************MAKE VIDEOS********************
'''

path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/ROCK_AVG2_10/evolution/plots/Density/manual/'
delay_number = 10
loop_number = 1
videoname = 'density_manual'
iplt.make_video(path, delay_number, loop_number, videoname)

path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/SAND_AVG2_10/evolution/plots/Density/manual/'
delay_number = 10
loop_number = 1
videoname = 'density_manual'
iplt.make_video(path, delay_number, loop_number, videoname)

path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/ROCK_AVG2_10/evolution/plots/Density/all/'
delay_number = 10
loop_number = 1
videoname = 'density_all'
iplt.make_video(path, delay_number, loop_number, videoname)

path = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo/kristina/data/SAND_AVG2_10/evolution/plots/Density/all/'
delay_number = 10
loop_number = 1
videoname = 'density_all'
iplt.make_video(path, delay_number, loop_number, videoname)