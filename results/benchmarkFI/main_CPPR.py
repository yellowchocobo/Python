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


pathm = ['/uio/kant/geo-ceed-u1/nilscp/Nils/Python/arcpy/', '/uio/kant/geo-ceed-u1/nilscp/Nils/Python/craterGeomorph/']


# add directories to system path
for pat in pathm:
    sys.path.append(pat)


import crater as cr
import plot as iplt
import geomorphCraters_for_iSALE as geom


# path to in-house scripts and pySALEPlot
path_pySALEPlot = '/work/nilscp/iSALE/Dellen/lib'

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
g = 1.62

'''
********************CALCULATE DIMENSION EVOLUTION******************************

It might be some problems about the evolution of the crater volume through time
depending on which model is used 
'''


pathm = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/'

folders =  ['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2', 
            'COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH']


#CPPR = ['CPPR5', 'CPPR10', 'CPPR20'] # (eventually to automatize also including
# the cppr)


for f in folders:
    
    path = pathm + f + '/'

    os.chdir(path)

    # selection of folders
    search = '*'
    models = glob.glob(search)

    pathdata = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/'
    pathplots = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/plots/'


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
    # evo.XY_PROFILE(folders, path, pathdata, 19, 3) # arbitrary
    
    ###############################################################################
    # 6. PLOTTING
    ###############################################################################    
    
    


'''
****************CALCULATE EXTRA GEOMORPHOLOGICAL PARAMETERS********************
''' 
for m in models:
    
        path_tmp = pathdata + m + '/'
        geom.calculation(path_tmp, paths)

'''
*********************************GENERATE PLOTS********************************
'''

###############################################################################
# 1. CRATER EVOLUTION IN FUNCTION OF CPPR and DIFFERENT TYPES OF MODEL SETUPS
###############################################################################

path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/CDILNOA/CDILNOA_L100/'
paths = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/CDILNOA_L100/'
normpath = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/CDILNOA/CDILNOA_L100/'
norm = 0 # 0 no normalization
showtransient = False
showhline = False
timei = 'pen' #'transient' or 'norm' or 'pen' (for penetration time)

lbl = r"$\mathit{U} = 12.7 km/s, \mathit{L} = 100 m$"

iplt.morphology(path1, paths, norm, normpath, timei, showtransient, showhline, lbl)


#path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/IVANOVS/IVANOVS_L100/'
#path2 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/COLDWPO/COLDWPO_L100/'
#path3 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/COLDNOP/COLDNOP_L100/'
#path4 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/SANDCOH/SANDCOH_L100/'


path1 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR5/CDILWPO/CDILWPO_L100/'
path2 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/CDILWPO/CDILWPO_L100/'
path3 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR20/CDILWPO/CDILWPO_L100/'
#path4 = '/uio/kant/geo-ceed-u1/nilscp/Desktop/stallo_work/benchmarkFI/CPPR/UAVG/CPPR10/SANDCOH/SANDCOH_L100/'

os.chdir(path1)
model1 = psp.opendatfile('jdata.dat')
den1 = model1.readStep('Den', 68)

os.chdir(path2)
model2 = psp.opendatfile('jdata.dat')
den2 = model2.readStep('Den', 68)

os.chdir(path3)
model3 = psp.opendatfile('jdata.dat')
den3 = model3.readStep('Den', 68)

os.chdir(path4)
model4 = psp.opendatfile('jdata.dat')
den4 = model4.readStep('Den', 400)

plt.figure(1)
plt.pcolormesh(model1.xc, model1.yc, den1.cmc[0])
plt.xlim((0,2000))
plt.ylim((-1500,200))

plt.figure(2)
plt.pcolormesh(model2.xc, model2.yc, den2.cmc[0])
plt.xlim((0,2000))
plt.ylim((-1500,200))

plt.figure(3)
plt.pcolormesh(model3.xc, model3.yc, den3.cmc[0])
plt.xlim((0,2000))
plt.ylim((-1500,200))

plt.figure(4)
plt.pcolormesh(model4.xc, model4.yc, den4.cmc[0])
plt.xlim((0,1500))
plt.ylim((-500,200))



model_setup =  ['COLDNOP', 'COLDWPO', 'COLDWPO2', 'IVANOVS', 'SANDCOH'] #['CDILWPO'] #['CDILNOA', 'CDILNOP', 'CDILWPA', 'CDILWPA2', 'CDILWPO', 'CDILWPO2'] 
                


for m in model_setup:
    f1 = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/' + m + '_L100/evolution/' + m + '_L100_data.txt'
    f2 = '/run/media/nilscp/Squall/benchmarkFI/CPPR5/data/' + m + '_L100/evolution/' + m + '_L100_data.txt'
    
    data1 = np.loadtxt(f1, comments="#")
    data2 = np.loadtxt(f2, comments="#")
                       
    t1 = data1[:,0]
    t2 = data2[:,0]
    v1 = data1[:,4]
    v2 = data2[:,4]
    vr1 = data1[:,-1]
    vr2 = data2[:,-1]
    
    ix1 = np.nanargmax(v1)
    ix2 = np.nanargmax(v2)
    
    plt.plot(t1,v1,"o", label = m + ' CPPR10')
    plt.plot(t2,v2,"o", label = m + ' CPPR5')
    plt.plot(t1[ix1], v1[ix1], '+')
    plt.plot(t2[ix2], v2[ix2], '+')

#f3 = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/' + 'CDILWPO' + '_L100_new/evolution/' + m + '_L100_data.txt'
#f4 = '/run/media/nilscp/Squall/benchmarkFI/CPPR10/data/' + 'CDILWPO' + '_L100_new/evolution/' + m + '_L100_data.txt'    
plt.legend()